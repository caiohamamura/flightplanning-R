Wysokosc <- 100
output = ("mytest/lot.csv")
roi = sf::st_read("mytest/lasek.gpkg")

# if( nrow(roi) > 1) {
#   for (i in seq_len(nrow(roi))) {
#     output <- paste0("mytest/lot_", i, ".csv")
#     flightplanning::litchi.plan(roi[i,],
#                                 output,
#                                 params,
#                                 flight.lines.angle = -1,
#                                 max.waypoints.distance = 4000,
#                                 max.flight.time = 16,
#                                 grid = FALSE
#     )
#   }
# }

if(nrow(roi) > 1) {
  roi <- sf::st_union(roi) |>
    sf::st_as_sf()
}

params = flightplanning::flight.parameters(height=Wysokosc,
                                           flight.speed.kmh=24,
                                           side.overlap = 0.8,
                                           front.overlap = 0.8)

# Create the csv plan
# flightplanning::litchi.plan(roi,
#   output,
#   params,
#   flight.lines.angle = -1,
#   max.waypoints.distance = 4000,
#   max.flight.time = 16,
#   grid = FALSE
# )

flight.params <- params
gimbal.pitch.angle = -90
flight.lines.angle = -1
max.waypoints.distance = 2000
max.flight.time = 15
starting.point = 1
grid = FALSE

litchi.plan = function(roi, output,
                       flight.params, gimbal.pitch.angle = -90,
                       flight.lines.angle = -1, max.waypoints.distance = 2000,
                       max.flight.time = 15, starting.point = 1, grid = FALSE) {
  # Check parameters
  roiCRS <- sf::st_crs(roi)

  if (class(roi)[1] == "sf") {
    roi <- sf::as_Spatial(roi)
  }
  if (class(roi)[1] != "SpatialPolygonsDataFrame")
    stop("ROI is not a valid polygon layer")
  if (length(grep("units=m", as.character(roi@proj4string@projargs))) == 0)
    stop("ROI is not in a metric projection")
  if (methods::is(flight.params)[1] != "Flight Parameters")
    stop("Flight parameters is not an instance returned from flight.parameters()")
  if (!is.logical(grid)) {
    stop("parallel has to be TRUE or FALSE")
  }

  # Parameters calculated
  flight.speed.kmh = flight.params@flight.speed.kmh
  flightSpeedMs = flight.speed.kmh / 3.6
  height = flight.params@height
  groundHeight = flight.params@ground.height
  groundHeightOverlap = groundHeight * flight.params@front.overlap
  flightLineDistance = flight.params@flight.line.distance
  vertices = roi@polygons[[1]]@Polygons[[1]]@coords

  # Get bounding box parameters
  if (flight.lines.angle != -1) {
    minBbox = flightplanning:::getBBoxAngle(vertices, flight.lines.angle)
  } else {
    # if angle not specified use minimum possible bounding box
    minBbox = shotGroups::getMinBBox(vertices)
  }

  if (grid == FALSE) {
    width = minBbox$height
    height = minBbox$width
    alpha = minBbox$angle
  } else {
    width = minBbox$width
    height = minBbox$height
    alpha = minBbox$angle-90
  }

  rads = alpha*pi/180
  centroid = apply(minBbox$pts, 2, mean)

  # Calculate points offset from centroid
  # based on angle and width/height offsets
  # width offsets (between flightlines)
  nLines = ceiling(width / flightLineDistance) + 1
  xWidths = (-nLines/2):(nLines/2) * flightLineDistance
  xWidths = rep(xWidths, each=2)

  # heights offset (one for upper half
  #                 one for lower half)
  heightDistance = groundHeight-groundHeightOverlap
  heightAdjusted = height + 2*heightDistance
  # Put offset to avoid intersection issues
  heightAdjusted = heightAdjusted + heightDistance*2
  heightMHalf = -heightAdjusted/2
  heightPHalf = heightAdjusted/2
  yHeights = c(heightMHalf, heightPHalf)


  # Switch position of the first point
  if (starting.point == 2) {
    yHeights = c(heightPHalf, heightMHalf)
  } else if (starting.point == 3) {
    xWidths = rev(xWidths)
    yHeights = c(heightPHalf, heightMHalf)
  } else if (starting.point == 4) {
    xWidths = rev(xWidths)
  }

  # Interleave one upper, two bottom
  # two upper, two bottom... until end
  yHeights = c(rep(c(yHeights, rev(yHeights)), nLines/2+1))
  yHeights = yHeights[1:length(xWidths)]


  # Calculate translated x and y from
  # angles and offsets from centroid
  xys = data.frame(
    x = -xWidths  * sin(rads) +
      yHeights * cos(rads),
    y = xWidths   * cos(rads) +
      yHeights  * sin(rads))

  # Initial waypoints to intersect waypoints
  waypoints = xys + rep(centroid, each=nrow(xys))

  #################################################
  # Intersect each flight line from bounding box
  # to match actual ROI
  #################################################
  # For some reason gIntersection with MULTILINESTRING
  # will return linestrings in inconsistent order
  # though it will be done in a for loop

  #
  #
  #  wktLines = paste(apply(waypoints, 1, paste, collapse=" "), collapse=", ")
  #  wktLines = paste("LINESTRING(", wktLines,")")
  #  gLines = rgeos::readWKT(wktLines, p4s = roi@proj4string)
  #  inter = rgeos::gIntersection(rgeos::gBuffer(roi, width = flightLineDistance), gLines)
  #  nLines = length(inter@lines[[1]]@Lines)

  #  flightLines = t(sapply(inter@lines[[1]]@Lines, function(x) x@coords))

  # RSB
  glist <- vector(mode="list", length=nrow(waypoints)-1)
  for (i in seq_along(glist)) glist[[i]] <- sp::Lines(list(sp::Line(waypoints[c(i, (i+1)),])), ID=as.character(i))
  gLines <- sp::SpatialLines(glist, proj4string=slot(roi, "proj4string"))

# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
  lines <- do.call(
    sf::st_sfc,
    lapply(
      1:(nrow(waypoints)-1),
      function(i) {
        sf::st_linestring(
          matrix(
            c(as.numeric(waypoints[i, ]), as.numeric(waypoints[i + 1, ])),
            ncol = 2, byrow = TRUE
          )
        )
      }
    )
  )

  lines <- lines |>
    sf::st_as_sf(crs = roiCRS)
  lines$ID <- seq.int(nrow(lines))
#  rename_geometry(lines, "geom")
  # class(lines)
  # sp::plot(roi)
  # terra::plot(lines, add = TRUE)
  mgLines <- lines |>
    sf::as_Spatial()

sp::plot(gLines)
sp::plot(mgLines, add = TRUE, col = "yellow")
class(gLines)
# ---------------------------------------------------------------------------------------------
  inter = rgeos::gIntersection(rgeos::gBuffer(roi, width = flightLineDistance), gLines, byid=TRUE)
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
  ginter <- sf::st_intersection(
    sf::st_buffer(sf::st_as_sf(roi), flightLineDistance),
    lines)
  # myinter <- ginter
  # ginter <- ginter |>
  #   sf::as_Spatial()
  # class(ginter)
  # terra::plot(myinter, add = TRUE, col = "yellow")

# ---------------------------------------------------------------------------------------------

  nLines <- length(inter)
  flightLines <- t(sapply(slot(inter, "lines"), function(x) slot(slot(x,  "Lines")[[1]], "coords")))
  # RSB
  flightLines = flightLines[,c(1,3,2,4)]
#  flightLines

# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
  nLines <- length(ginter)
  # gflightLines <-
    a <- sf::st_coordinates(ginter)[,1:2]
    # to poniżej nie potrzebne gdyż waypoints == a
  #   gflightLines = matrix(nrow = mnLines, ncol = 4)
  #   gflightLines[,1] <- a[seq(1,length(a[,1]),2),1]
  #   gflightLines[,2] <- a[seq(1,length(a[,1]),2),2]
  #   gflightLines[,3] <- a[seq(2,length(a[,1]),2),1]
  #   gflightLines[,4] <- a[seq(2,length(a[,1]),2),2]
  # flightLines <- gflightLines
# ---------------------------------------------------------------------------------------------

  waypoints = matrix(nrow=nLines * 2, ncol=2)
  waypoints[seq(1, nLines*2, 2),] = flightLines[, 1:2]
  waypoints[seq(2, nLines*2, 2),] = flightLines[, 3:4]
waypoints
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

waypoints <- a
# ---------------------------------------------------------------------------------------------


  # Calculate curves points to allow smooth curves
  curvedPoints = flightplanning:::outerCurvePoints(waypoints = waypoints,
                                  angle = alpha,
                                  flightLineDistance = flightLineDistance)

  # Adjust curve points position to avoid acute angles
  adjustedCurves = flightplanning:::adjustAcuteAngles(xy = curvedPoints,
                                     angle = alpha,
                                     minAngle = 80)

  # Concatenate regular waypoints with curve waypoints
  wgs84 = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  wptsMatrix = as.data.frame(matrix(nrow=nrow(waypoints)+nrow(adjustedCurves),ncol=4))
  colnames(wptsMatrix) = colnames=c("x", "y", "isCurve", "takePhoto")
  mat_pos = 1
  for (i in seq_len(nrow(waypoints))) {
    curve = as.vector(adjustedCurves[as.character(i),])
    hasCurve = !anyNA(curve)
    if (hasCurve) {
      if (curve$before) {
        wptsMatrix[mat_pos,] = c(curve[1:2], TRUE, FALSE)
        mat_pos = mat_pos + 1
        wptsMatrix[mat_pos,] = c(waypoints[i, 1:2], FALSE, i %% 2 == 1)
        mat_pos = mat_pos + 1
      } else {
        wptsMatrix[mat_pos,] = c(waypoints[i, 1:2], FALSE, i %% 2 == 1)
        mat_pos = mat_pos + 1
        wptsMatrix[mat_pos,] = c(curve[1:2], TRUE, FALSE)
        mat_pos = mat_pos + 1
      }
    } else {
      wptsMatrix[mat_pos,] = c(waypoints[i, 1:2], FALSE, i %% 2 == 1)
      mat_pos = mat_pos + 1
    }
  }
  waypoints = wptsMatrix

  # Break if distance greater than the maxWaypointDistance
  waypointsXY = waypoints[, c("x", "y")]
  distances = sqrt(diff(waypoints$x)**2 + diff(waypoints$y)**2)
  breakIdx = distances > max.waypoints.distance

  newSize = nrow(waypoints) + sum(breakIdx)
  if (newSize != nrow(waypoints)) {
    midpoints = (waypointsXY[breakIdx,] + waypointsXY[-1,][breakIdx,])/2
    waypoints2 = data.frame(x = numeric(newSize),
                            y = numeric(newSize),
                            isCurve = FALSE,
                            takePhoto = TRUE)

    pos = seq_along(breakIdx)[breakIdx]
    idx = pos + order(pos)
    waypoints2[idx,1:2] = midpoints
    waypoints2[-idx,] = waypoints
    waypoints = waypoints2
  }


# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
  t <- waypoints |>
    sf::st_as_sf(coords = c("x", "y"), crs = roiCRS) |>
    sf::st_transform(crs = "EPSG:4326")
  lngs <- as.numeric(sf::st_coordinates(t)[,1])
  lats <- as.numeric(sf::st_coordinates(t)[,2])
# ---------------------------------------------------------------------------------------------


  # Transform to WGS84 latitude and longitude
  transform = rgdal::rawTransform(roi@proj4string@projargs, wgs84, n=nrow(waypoints), x=waypoints[,1], y=waypoints[,2])
  lats = transform[[2]]
  lngs = transform[[1]]
  graphics::plot(waypoints[,1:2])
  graphics::polygon(roi@polygons[[1]]@Polygons[[1]]@coords)


  # Calculate heading
  nWaypoints = nrow(waypoints)
  latDiff = lats[-1]-lats[-nWaypoints]
  lngDiff = lngs[-1]-lngs[-nWaypoints]
  headDegree = atan(latDiff/lngDiff)/pi*180
  finalHeading = 270-headDegree
  finalHeading[lngDiff > 0] = 90-headDegree[lngDiff > 0]

  # Set parameters of the flight in the CSV
  dfLitchi = read.csv(system.file("extdata/litchi.csv", package = "flightplanning"))
  dfLitchi = dfLitchi[rep(1, length(lats)),]
  dfLitchi$latitude = lats
  dfLitchi$longitude = lngs
  dfLitchi$altitude.m. = flight.params@height
  dfLitchi$speed.m.s. = flightSpeedMs
  dfLitchi$heading.deg. = c(finalHeading, 90)
  dfLitchi$curvesize.m. = 0
  dfLitchi$curvesize.m.[waypoints$isCurve==1] = flightLineDistance*0.5
  dfLitchi$photo_distinterval = flight.params@ground.height
  dfLitchi$gimbalpitchangle = gimbal.pitch.angle


  # Split the flight if is too long
  dists = sqrt(diff(waypoints[,1])**2+diff(waypoints[,2])**2)
  distAcum = c(0,cumsum(dists))
  flightTime = distAcum / (flightSpeedMs*0.75) / 60
  finalSize = nrow(dfLitchi)
  totalFlightTime = flightTime[finalSize]
  dfLitchi$split = 1
  if (totalFlightTime > max.flight.time) {
    indexes = seq_len(finalSize)
    nBreaks = ceiling(totalFlightTime/max.flight.time)
    breaks = seq(0, flightTime[finalSize], length.out = nBreaks+1)[c(-1, -nBreaks-1)]
    endWaypointsIndex = indexes[waypoints$isCurve & (seq_len(finalSize) %% 2 == 0)]
    endWaypoints = flightTime[waypoints$isCurve & (seq_len(finalSize) %% 2 == 0)]
    selected = sapply(breaks, function(x) which.min(abs(endWaypoints-x)))
    waypointsBreak = endWaypointsIndex[indexes[selected]]


    dfLitchi$split = rep(1:nBreaks, diff(c(0, waypointsBreak, finalSize)))
    splits = split.data.frame(dfLitchi, f = dfLitchi$split)
    message("Your flight was splitted in ", length(splits), " splits,
because the total time would be ", round(totalFlightTime, 2), " minutes.")
    message("They were saved as:")
    first = substr(output, 1, nchar(output)-4)
    second = substr(output, nchar(output)-3, nchar(output))
    for (dataSplit in splits) {
      i = dataSplit[1, ]$split
      output2 = paste0(first, i, second)
      write.csv(dataSplit[,-ncol(dataSplit)], output2, row.names = FALSE)
      message(output2)
    }
    output2 = paste0(first, "_entire", second)
    write.csv(dfLitchi, output2, row.names = FALSE)
    message("The entire flight plan was saved as:")
    message(output2)
  } else {
    write.csv(dfLitchi, output, row.names = FALSE)
  }

  colors = grDevices::rainbow(length(unique(dfLitchi$split)))
  for (i in unique(dfLitchi$split))
  {
    graphics::lines(waypoints[dfLitchi$split == i,1:2], lty=2, col=colors[as.integer(i)])
  }
  graphics::text(waypoints[,1], waypoints[,2], seq_along(waypoints[,1]), pos=3)


  message("#####################")
  message("## Flight settings ## ")
  message("#####################")
  message("Min shutter speed: ", appendLF = FALSE)
  message(flight.params@minimum.shutter.speed)
  message("Photo interval:    ", appendLF = FALSE)
  message(flight.params@photo.interval, appendLF = FALSE)
  message(" s")
  message("Flight speed:      ", appendLF = FALSE)
  message(round(flight.params@flight.speed.kmh, 4), appendLF = FALSE)
  message(" km/h")
  message("Flight lines angle: ", appendLF = FALSE)
  message(round(alpha, 4))
  message('Total flight time: ', appendLF = FALSE)
  message(round(totalFlightTime, 4))

  return (waypoints)
}

rename_geometry <- function(df, name){
  current = attr(df, "sf_column")
  names(df)[names(df) == current] = name
  sf::st_geometry(df) = name
  df
}

f <- list.files(path = "mytest", pattern = ".csv", full.names = TRUE)
unlink(f)
