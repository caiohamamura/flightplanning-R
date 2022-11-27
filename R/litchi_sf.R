MAX_WAYPOINTS = 99


#'  Function to generate Litchi csv flight plan
#'
#' @rdname litchi_sf
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
#' @param launch list(0,0) launch point coordinates (x, y), has to be provided in the same metric CRS as roi
#' @param grid logical (default FALSE). Change direction of the fly lines over polygon from parallel to perpendicular
#' @param distancemethod logical (default FALSE). Change shutter interval from time to distance
#'
#' @examples
#' library(flightplanning)
#'
#' exampleBoundary =
#'   sf::st_read(
#'     system.file("extdata", "exampleBoundary.shp", package="flightplanning")) |>
#'     sf::st_as_sf()
#'
#' outPath = tempfile(fileext=".csv")
#'
#' flight.params = flightplanning::flight.parameters(
#'   gsd = 4,
#'   side.overlap = 0.8,
#'   front.overlap = 0.8,
#'   flight.speed.kmh = 54
#' )
#'
#' litchi_sf(exampleBoundary,
#'             outPath,
#'             flight.params,
#'             flight.lines.angle = -1,
#'             max.waypoints.distance = 2000,
#'             max.flight.time = 15)
#'
#'
#' @import sf
#' @importFrom graphics text
#' @importFrom shotGroups getMinBBox
#' @importFrom methods slot
#' @importFrom utils data read.csv write.csv
#'
#' @export
litchi_sf = function(roi,
                     output,
                     flight.params,
                     gimbal.pitch.angle = -90,
                     flight.lines.angle = -1,
                     max.waypoints.distance = 2000,
                     max.flight.time = 15,
                     starting.point = 1,
                     launch = list(0, 0),
                     grid = FALSE,
                     distancemethod = FALSE) {

  # Check parameters
  if (class(roi)[1] != "sf") {
    roi <- sf::st_as_sf(roi)
  }

  if(nrow(roi) > 1) {
    roi <- roi[1,] |>
      sf::st_as_sf()
  }

  if(!sf::st_geometry_type(roi)[[1]] %in% c("POLYGON", "MULTIPOLYGON")) {
    stop("ROI is neither POLYGON nor MULTIPOLYGON")
  }
  if (!grepl("LENGTHUNIT[\"metre\",1]", sf::st_crs(roi)[2], fixed = TRUE))
    stop("ROI is not in a metric projection")
  if (methods::is(flight.params)[1] != "Flight Parameters")
    stop("Flight parameters is not an instance returned from flight.parameters()")
  if (!is.logical(grid)) {
    stop("grid has to be TRUE or FALSE")
  }
  if (!is.logical(distancemethod)) {
    stop("distancemethod has to be TRUE or FALSE")
  }

  # Parameters calculated
  flight.speed.kmh = flight.params@flight.speed.kmh
  flightSpeedMs = flight.speed.kmh / 3.6
  height = flight.params@height
  groundHeight = flight.params@ground.height
  groundHeightOverlap = groundHeight * flight.params@front.overlap
  flightLineDistance = flight.params@flight.line.distance
  vertices <- sf::st_coordinates(roi)[,1:2]
  roiCRS <- sf::st_crs(roi)

  # Get bounding box parameters
  if (flight.lines.angle != -1) {
    minBbox = getBBoxAngle(vertices, flight.lines.angle)
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
  # Then need to update flightLineDistance to avoid offset/quantization errors
  flightLineDistance = width / (nLines - 1)
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
  if (starting.point == 0) {
    # In this case we will automatically pick the best starting point
    # TODO check if launch is valid and not (0,0)
    # TODO figure out closest corner in shape to launch point and then set starting.point to 1-4
    starting.point == 1
  }
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

  # ---------------------------------------------------------------------------------------------
  inter <-suppressWarnings(sf::st_intersection(
    sf::st_buffer(sf::st_as_sf(roi), flightLineDistance),
    lines) |>
      sf::st_cast(to = "LINESTRING"))

  # ---------------------------------------------------------------------------------------------
  waypoints <- sf::st_coordinates(inter)[,1:2]

  # ---------------------------------------------------------------------------------------------

  # Calculate curves points to allow smooth curves
  curvedPoints = outerCurvePoints(waypoints = waypoints,
                                  angle = alpha,
                                  flightLineDistance = flightLineDistance)
  row.names(curvedPoints) <- curvedPoints$index

  # Adjust curve points position to avoid acute angles
  adjustedCurves = adjustAcuteAngles(xy = curvedPoints,
                                     angle = alpha,
                                     minAngle = 80)
  row.names(adjustedCurves) <- adjustedCurves$index

  # Concatenate regular waypoints with curve waypoints
  wptsMatrix = as.data.frame(matrix(nrow=nrow(waypoints)+nrow(adjustedCurves),ncol=4))
  colnames(wptsMatrix) = colnames=c("x", "y", "isCurve", "takePhoto")
  mat_pos = 1
  for (i in seq_len(nrow(waypoints))) {
    #    i <- 6
    #    mat_pos <- 11
    curve = as.vector(adjustedCurves[as.character(i), ])
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

  # from https://github.com/caiohamamura/flightplanning-R/pull/4/commits/c36501d8ea0cab6cdcb3e5123452bb4c18599aac
  #
  # Break if distance greater than the maxWaypointDistance
  # A single pass only adds one intermediate waypoint even if a leg is longer than max,
  # but more than one intermediate point may be needed.
  # We can iterate this process as a temp fix but we may be adding more intermediate
  # waypoints than strictly necessary--e.g. when 2 intermediate points will suffice we will get 3.
  retest = TRUE
  while (retest) {
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
    else {
      retest = FALSE
    }
  }
  # from: https://github.com/caiohamamura/flightplanning-R/pull/4/commits/fed8f6928aef89eb5d7884b6001b7e16f0ec4bfd
  #
  # Check if launch point has been specified before inserting it as way-point 1
  hasCustomLaunch = (launch[1] != 0) || (launch[2] != 0)
  if (hasCustomLaunch) {
    message("Launch point specified: ", launch[1], ',', launch[2])
    MAX_WAYPOINTS = MAX_WAYPOINTS - 1
  } else {
    message("No launch point specified")
  }

  # ---------------------------------------------------------------------------------------------
  t <- waypoints |>
    sf::st_as_sf(coords = c("x", "y"), crs = roiCRS) |>
    sf::st_transform(crs = "EPSG:4326")
  lngs <- as.numeric(sf::st_coordinates(t)[,1])
  lats <- as.numeric(sf::st_coordinates(t)[,2])
  photos = t$takePhoto
  # ---------------------------------------------------------------------------------------------

  graphics::plot(waypoints[,1:2])
  graphics::polygon(sf::st_coordinates(roi))

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
  dfLitchi$altitudemode = 1 # AGL handling
  dfLitchi$speed.m.s. = flightSpeedMs
  dfLitchi$heading.deg. = c(finalHeading, 90)
  dfLitchi$curvesize.m. = 0
  dfLitchi$curvesize.m.[waypoints$isCurve==1] = flightLineDistance*0.5
  if (!isTRUE(distancemethod)) {
    dfLitchi$photo_timeinterval = flight.params@photo.interval * photos
  } else {
    dfLitchi$photo_distinterval = flight.params@photo.interval * flightSpeedMs * photos
  }
  dfLitchi$gimbalmode = 2 # Interpolate gimbal position between waypoints
  dfLitchi$gimbalpitchangle = gimbal.pitch.angle

  # Split the flight if is too long
  dists = sqrt(diff(waypoints[,1])**2+diff(waypoints[,2])**2)
  distAcum = c(0,cumsum(dists))
  flightTime = distAcum / (flightSpeedMs*0.75) / 60
  finalSize = nrow(dfLitchi)
  totalFlightTime = flightTime[finalSize]
  dfLitchi$split = 1
  if ((totalFlightTime > max.flight.time) || (nrow(waypoints) > MAX_WAYPOINTS)) {
    indexes = seq_len(finalSize)
    nBreaks = max(ceiling(totalFlightTime/max.flight.time), ceiling(nrow(waypoints)/MAX_WAYPOINTS))
    breaks = seq(0, flightTime[finalSize], length.out = nBreaks+1)[c(-1, -nBreaks-1)]
    endWaypointsIndex = indexes[waypoints$isCurve & (seq_len(finalSize) %% 2 == 0)]
    endWaypoints = flightTime[waypoints$isCurve & (seq_len(finalSize) %% 2 == 0)]
    selected = sapply(breaks, function(x) which.min(abs(endWaypoints-x)))
    waypointsBreak = endWaypointsIndex[indexes[selected]]

    dfLitchi$split = rep(1:nBreaks, diff(c(0, waypointsBreak, finalSize)))
    splits = split.data.frame(dfLitchi, f = dfLitchi$split)


    if (hasCustomLaunch) {
      p0x = launch[[1]][1]
      p0y = launch[[2]][1]

      message("adding custom launch point to submissions")

      launch84 <- sf::st_point(x = c(launch[[1]], launch[[2]])) |>
        sf::st_sfc(crs = roiCRS) |>
        sf::st_transform(crs = "EPSG:4326") |>
        sf::st_coordinates()

      overage = NULL

      for (i in 1:length(splits)) {
        message("starting ", i)
        if (!is.null(overage)) {
          message("setting ", i, " to ", "overage (", nrow(overage), ") and ", nrow(splits[[i]]))
          splits[[i]] = rbind(overage, splits[[i]])
        }

        mercator = splits[[i]] |>
          sf::st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326") |>
          sf::st_transform(crs = roiCRS) |>
          subset(select = "geometry") |>
          sf::st_coordinates()

        p1x = as.numeric(mercator[1, 1])
        p1y = as.numeric(mercator[1, 2])

        dx = p1x - p0x
        dy = p1y - p0y
        distance = earth.dist(launch84[[1]][1], launch84[[2]][1], splits[[i]]$longitude[1], splits[[i]]$latitude[1])

        interpPtsToAdd = floor(distance / max.waypoints.distance)
        nPtsToAdd = 1 + interpPtsToAdd

        message("adding ", nPtsToAdd, " points")

        ptsToAdd = rbind(splits[[1]][1:nPtsToAdd,])
        ptsToAdd$split <- i
        ptsToAdd$curvesize.m. <- 0
        ptsToAdd$photo_distinterval <- 0
        ptsToAdd$photo_timeinterval <- 0

        toConvert = data.frame(
          lon = numeric(nPtsToAdd),
          lat = numeric(nPtsToAdd)
        )

        toConvert[1,] = c(p0x, p0y)
        if (nPtsToAdd > 1) {
          for (j in 2:nPtsToAdd) {
            toConvert[j,] <- c(p0x + ((j - 1) / nPtsToAdd) * dx, p0y + ((j - 1) / nPtsToAdd) * dy)
          }
        }

        wgs84D = toConvert |>
          sf::st_as_sf(coords = c("lon", "lat"), crs = roiCRS) |>
          sf::st_transform(crs = "EPSG:4326") |>
          subset(select = "geometry") |>
          sf::st_coordinates()

        ptsToAdd$latitude = wgs84D[ ,2]
        ptsToAdd$longitude = wgs84D[ ,1]

        splitSize = nrow(splits[[i]])
        totalSize = splitSize + nPtsToAdd
        rem = 0
        if (totalSize > MAX_WAYPOINTS + 1) {
          rem = totalSize - (MAX_WAYPOINTS + 1)
        }

        if (rem > 0) {
          message("setting overage to ", splitSize + 1 - rem, " : ", splitSize)
          message(class(splits[[i]]))
          message(splits[[i]][splitSize + 1 - rem:splitSize,])
          message(colnames(splits[[i]]))
          message(splits[[i]])
          message(colnames(splits[[i]][splitSize + 1 - rem:splitSize,]))
          message(rownames(splits[[i]][splitSize + 1 - rem:splitSize,]))
          overage = rbind(splits[[i]][splitSize + 1 - rem:splitSize,])
          message("overage has ", nrow(overage))
          message(overage)
          message(overage[splitSize + 1 - rem: splitSize,])
        } else {
          message("setting overage to NULL")
          overage = NULL
        }

        message("setting ", i, " to ptsToAdd (", nrow(ptsToAdd), ") + splits 1 : ", splitSize - rem)
        splits[[i]] = rbind(ptsToAdd, splits[[i]][1:splitSize - rem,])
      }

      if (!is.null(overage)) {
        newIdx = length(splits) + 1
        splits[[newIdx]] = rbind(overage)
        splits[[newIdx]]$split = newIdx
      }
    }

    if (nrow(waypoints) > MAX_WAYPOINTS) {
      message("Your flight was split into ", length(splits), " sub-flights, because the number of waypoints ", nrow(waypoints), " exceeds the maximum of ", MAX_WAYPOINTS, ".")
    }
    else {
      message("Your flight was split into ", length(splits), " sub-flights, because the total flight time of ", round(totalFlightTime, 2), " minutes exceeds the max of ", max.flight.time, " minutes.")
    }
    message("The flights were saved as:")
    first = paste0(substr(output, 1, nchar(output)-4), "_")
    second = substr(output, nchar(output)-3, nchar(output))
    for (dataSplit in splits) {
      i = dataSplit[1, ]$split
      output2 = paste0(first, i, second)
      write.csv(dataSplit[,-ncol(dataSplit)], output2, row.names = FALSE)
      message(output2)
    }
    output2 = paste0(first, "entire", second)
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

  #update flight.params
  flight.params@number.of.waypoints <- nWaypoints
  flight.params@alpha <- alpha
  flight.params@total.flight.time <- totalFlightTime

  # messages
  flight.summary(flight.params)

  return (waypoints)
}
