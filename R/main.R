MIN_PHOTO_INTERVAL = 2
DIAG_35MM = sqrt(36^2 + 24^2) # Classical 35mm film diagonal
MAX_WAYPOINTS = 99

#'  Function to generate Litchi csv flight plan
#'
#'
#' @rdname litchi.plan
#'
#' @return A data frame with the waypoints calculated for the flight plan
#'
#' @param roi range of interest loaded as an OGR layer, must be in
#' a metric units projection for working properly
#' @param output output path for the csv file
#' @param flight.params Flight Parameters. parameters calculated from flight.parameters()
#' @param gimbal.pitch.angle gimbal angle for taking photos, default -90 (overriden at flight time)
#' @param flight.lines.angle angle for the flight lines, default -1 (auto set based on larger direction)
#' @param max.waypoints.distance maximum distance between waypoints in meters,
#' default 2000 (some issues have been reported with distances > 2 Km)
#' @param max.flight.time maximum flight time. If mission is greater than the estimated
#' time, it will be splitted into smaller missions.
#' @param starting.point numeric (1, 2, 3 or 4). Change position from which to start the flight, default 1
#' @param grid logical (default FALSE). Change direction of the fly lines over polygon from parallel to perpendicular
#'
#' @note this function will feed the csv flight plan with the `gimbal.pitch.angle`
#' and the `photo time interval` for each waypoint, but those are not supported
#' by Litchi yet, although they are present in the exported csv from the
#' Litchi hub platform, though it may be supported in the future; when it does
#' the function will already work with this feature.
#'
#' @examples
#' library(flightplanning)
#' library(rgdal)
#'
#' exampleBoundary = readOGR(
#'                           system.file("extdata",
#'                                       "exampleBoundary.shp",
#'                                       package="flightplanning"
#'                                      ),
#'                           "exampleBoundary")
#' outPath = tempfile(fileext=".csv")
#'
#' flight.params = flight.parameters(
#'   gsd = 4,
#'   side.overlap = 0.8,
#'   front.overlap = 0.8,
#'   flight.speed.kmh = 54
#' )
#'
#' litchi.plan(exampleBoundary,
#'             outPath,
#'             flight.params,
#'             flight.lines.angle = -1,
#'             max.waypoints.distance = 2000,
#'             max.flight.time = 15)
#'
#'
#' @export
#' @import rgeos rgdal
#' @importFrom graphics text
#' @importFrom shotGroups getMinBBox
#' @importFrom methods slot
#' @importFrom sp Line Lines SpatialLines
#' @importFrom utils data read.csv write.csv
#'
litchi.plan = function(roi,
                       output,
                       flight.params,
                       gimbal.pitch.angle = -90,
                       flight.lines.angle = -1,
                       max.waypoints.distance = 2000,
                       max.flight.time = 15,
                       starting.point = 1,
                       grid = FALSE) {
  # Check parameters
  if (class(roi)[1] == "sf") {
    roi <- sf::as_Spatial(roi)
  }

  if (class(roi)[1] != "SpatialPolygonsDataFrame")
    stop("ROI is not a valid polygon layer")
  if (length(grep("units=m", as.character(roi@proj4string@projargs))) == 0)
    stop("ROI is not in a metric projection")
  if (methods::is(flight.params)[1] != "Flight Parameters")
    stop("Flight parameters is not an instance returned from flight.parameters()")
  # TODO add a test
  if (!is.logical(grid)) {
    stop("grid has to be TRUE or FALSE")
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
    minBbox = getBBoxAngle(vertices, flight.lines.angle)
  } else {
    # if angle not specified use minimum possible bounding box
    minBbox = shotGroups::getMinBBox(vertices)
#    minBbox = getMinBBox(vertices)
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
  inter = rgeos::gIntersection(rgeos::gBuffer(roi, width = flightLineDistance), gLines, byid=TRUE)
  nLines <- length(inter)
  flightLines <- t(sapply(slot(inter, "lines"), function(x) slot(slot(x,  "Lines")[[1]], "coords")))

# RSB
  flightLines = flightLines[,c(1,3,2,4)]


  waypoints = matrix(nrow=nLines * 2, ncol=2)
  waypoints[seq(1, nLines*2, 2),] = flightLines[, 1:2]
  waypoints[seq(2, nLines*2, 2),] = flightLines[, 3:4]

  # Calculate curves points to allow smooth curves
  curvedPoints = outerCurvePoints(waypoints = waypoints,
                  angle = alpha,
                  flightLineDistance = flightLineDistance)

  # Adjust curve points position to avoid acute angles
  adjustedCurves = adjustAcuteAngles(xy = curvedPoints,
                   angle = alpha,
                   minAngle = 80)

  # Concatenate regular waypoints with curve waypoints
  wgs84 = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  wptsMatrix = as.data.frame(matrix(nrow=nrow(waypoints)+nrow(adjustedCurves),ncol=4))
  colnames(wptsMatrix) = colnames=c("x", "y", "isCurve", "takePhoto")
  mat_pos = 1
  for (i in seq_len(nrow(waypoints))) {
    curve = as.vector(adjustedCurves[as.character(i),])

    # From R 4.2.0 onwards:
    # "as.vector() gains a data.frame method which returns a simple
    # named list, also clearing a long standing ‘FIXME’ to enable
    # as.vector(<data.frame>, mode="list").  This breaks code relying
    # on as.vector(<data.frame>) to return the unchanged data frame."
    # therefore changing curve[,1:2] to curve[1,2] and removing cbind

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


  # Transform to WGS84 latitude and longitude
  transform = rgdal::rawTransform(roi@proj4string@projargs, wgs84, n=nrow(waypoints), x=waypoints[,1], y=waypoints[,2])
  lats = transform[[2]]
  lngs = transform[[1]]
  photos = waypoints[,4]
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
  dfLitchi$altitudemode = 1 # Above ground (Fixed AGL)
  dfLitchi$speed.m.s. = flightSpeedMs
  dfLitchi$heading.deg. = c(finalHeading, 90)
  dfLitchi$curvesize.m. = 0
  dfLitchi$curvesize.m.[waypoints$isCurve==1] = flightLineDistance*0.5
#  dfLitchi$photo_distinterval = flight.params@photo.interval * flightSpeedMs
  dfLitchi$photo_timeinterval = flight.params@photo.interval * photos
  dfLitchi$gimbalmode = 2 # Interpolated between waypoints
  dfLitchi$gimbalpitchangle = gimbal.pitch.angle
#  dfLitchi$actiontype1 = 5
#  dfLitchi$actionparam1 = gimbal.pitch.angle

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
  message('Distance between flyin lines: ', appendLF = FALSE)
  message(round(flightLineDistance, 3))
  message('', appendLF = TRUE)

  return (waypoints)
}



#' Function to calculate flight parameters
#'
#' This function will calculate the flight parameters by providing the camera settings
#' target flight height or gsd, front and side overlap.
#'
#' @rdname flight.parameters
#'
#' @param gsd target ground resolution in centimeters, must provide either `gsd` or `height`
#' @param height target flight height, default NA
#' @param focal.length35 numeric. Camera focal length 35mm equivalent, default 20
#' @param image.width.px numeric. Image width in pixels, default 4000
#' @param image.height.px numeric. Image height in pixels, default 3000
#' @param side.overlap desired width overlap between photos, default 0.8
#' @param front.overlap desired height overlap between photos, default 0.8
#' @param flight.speed.kmh flight speed in km/h, default 54.
#' @param max.gsd maximum ground resolution
#'
#' @examples
#' params = flight.parameters(
#'   gsd = 4,
#'   side.overlap = 0.8,
#'   front.overlap = 0.8,
#'   flight.speed.kmh = 54,
#'   max.gsd = 0
#' )
#'
#' @export
flight.parameters = function(
  height = NA,
  gsd = NA,
  focal.length35 = 20,
  image.width.px = 4000,
  image.height.px = 3000,
  side.overlap = 0.8,
  front.overlap = 0.8,
  flight.speed.kmh = 54,
  max.gsd = 0) {

  if (is.na(gsd) & is.na(height)) {
    stop("You must specify either gsd or height!")
  }

  image.diag.px = sqrt(image.width.px^2 + image.height.px^2)
  if (is.na(gsd)) {
    mult.factor = (height / focal.length35)
    diag.ground = DIAG_35MM * mult.factor
    gsd = diag.ground / image.diag.px * 100
    if ((max.gsd != 0) && (gsd > max.gsd)) {
      height = height * max.gsd / gsd
      warning(paste0("GSD of ", gsd, " is above target of ", max.gsd, " so adjusting height down to ", height))
      # Repeat as a Warning message because warnings are not always getting through
      message("WARNING: GSD of ", gsd, " is above target of ", max.gsd, " so adjusting height down to ", height)
      mult.factor = (height / focal.length35)
      diag.ground = DIAG_35MM * mult.factor
      gsd = diag.ground / image.diag.px * 100
      message("Final GSD is ", gsd)
    }
    groundWidth = image.width.px * gsd / 100
  } else {
    groundWidth = image.width.px * gsd / 100
    diag.ground = image.diag.px * gsd / 100
    mult.factor = diag.ground / DIAG_35MM
    height = mult.factor * focal.length35
  }

  flightLineDistance = groundWidth - side.overlap * groundWidth

  flightSpeedMs = flight.speed.kmh / 3.6
  speedPxPerSecond = flightSpeedMs / (gsd*0.01)

  # FIGUEIREDO, E. O. et al.
  # Planos de Voo Semiautônomos para Fotogrametria
  # com Aeronaves Remotamente Pilotadas de Classe 3
  maxPixelRoll = 1.2
  minimumShutterSpeed = paste("1/",round(speedPxPerSecond/maxPixelRoll), sep="")

  groundHeight = image.height.px * gsd / 100
  groundHeightOverlap = groundHeight * front.overlap
  groundAllowedOffset = groundHeight - groundHeightOverlap
  photoInterval = groundAllowedOffset / flightSpeedMs
  if (photoInterval < MIN_PHOTO_INTERVAL) {
    photoInterval = MIN_PHOTO_INTERVAL
    flightSpeedMs = groundAllowedOffset / photoInterval
    flight.speed.kmh = flightSpeedMs * 3.6
    warning(paste0("Speed had to be lowered because frequency of photos would be too high
        New speed: ", flight.speed.kmh, " km/h"))
    # Repeat as a Warning message because warnings are not always getting through
    message("WARNING: Speed had to be lowered because frequency of photos would be too high
        New speed: ", flight.speed.kmh, " km/h")
  } else if ((photoInterval %% .1) > 1e-4) {
    # Allow 0.1s resolution because integer seconds blocks useful drone speeds
    photoInterval = ceiling(photoInterval * 10) / 10
    flightSpeedMs = groundAllowedOffset / photoInterval
    flight.speed.kmh = flightSpeedMs*3.6
    warning(paste0("Speed lowered to ", flight.speed.kmh, " km/h to round up photo interval time
                   to ", photoInterval, " seconds"))
    # Repeat as a Warning message because warnings are not always getting through
    message("WARNING: Speed lowered to ", flight.speed.kmh, " km/h to round up photo interval
            time to ", photoInterval, " seconds")
  }
  params = methods::new("Flight Parameters")
  params@height = height
  params@gsd = gsd
  params@flight.line.distance = flightLineDistance
  params@minimum.shutter.speed = minimumShutterSpeed
  params@photo.interval = photoInterval
  params@ground.height = groundHeight
  params@front.overlap = front.overlap
  params@flight.speed.kmh = flight.speed.kmh

  return (params)
}

# # Example
# #
# (params = flight.parameters(
#   gsd=5,
#   flight.speed.kmh=52,
#   side.overlap = 0.8, #mudar para overlapFront
#   front.overlap = 0.8 #mudar para overpSide
# ))
