MIN_PHOTO_INTERVAL = 2

#'  Function to generate Litchi csv flight plan
#'
#'
#' @rdname litchi.plan
#'
#' @param roi range of interest loaded as an OGR layer, must be in
#' a metric units projection for working properly
#' @param output output path for the csv file
#' @param sensor.width numeric. sensor width in mm, default 6.17
#' @param focal.length35 numeric. focal length equivalent to 35mm, default 20
#' @param aspect.ratio character. Aspect ratio of the camera, default 4:3
#' @param image.width.px integer. The number of pixels captured in the width direction, default 4000
#' @param gsd numeric. Target ground resolution in centimeters, must specify either gsd or height, default NA
#' @param height numeric. target flight height, default NA
#' @param flight.speed.kmh flight speed in km/h, default 54
#' @param side.overlap desired width overlap between photos, default 0.8
#' @param front.overlap desired height overlap between photos, defualt 0.8
#' @param gimbal.pitch.angle gimbal angle for taking photos, default -90 (can be overriden at flight time)
#' @param flight.lines.angle angle for the flight lines, default -1 (auto set based on larger dimensions)
#' @param max.waypoints.distance maximum distance between waypoints in meters,
#' default 2000 (some issues have been reported with distances > 2km)
#' @param max.flight.time maximum flight time. If mission is greater than the estimated
#' time, it will be splitted into smaller missions.
#' @param starting.point numeric (1, 2, 3 or 4). Change position from which to start the flight, default 1
#'
#' @note this function will feed the csv flight plan with the `gimbal.pitch.angle`
#' and the photo time interval for each waypoint, but those are not supported
#' by Litchi yet, although they are present in the exported csv inside the
#' Litchi hub platform, though it may be supported in the future; when it does
#' the function will already work with this feature.
#'
#' @examples
#' library(flightplanning)
#' data(exampleBoundary)
#' outPath = tempfile(fileext=".csv")
#' litchi.plan(roi = exampleBoundary,
#'             output = outPath,
#'             gsd = 5,
#'             flight.speed.kmh = 54,
#'             side.overlap = 0.8,
#'             front.overlap = 0.8,
#'             flight.lines.angle = -1,
#'             max.waypoints.distance = 2000,
#'             max.flight.time = 15)
#'
#'
#' @export
#' @import sp rgeos rgdal
#' @importFrom graphics text
#' @importFrom utils data read.csv write.csv
litchi.plan = function(roi, output, sensor.width = 6.17,
                              focal.length35=20, aspect.ratio = "4:3",
                              image.width.px = 4000, gsd = NA, height = NA,
                              flight.speed.kmh = 54, side.overlap = 0.8,
                              front.overlap = 0.8, gimbal.pitch.angle = -90,
                              flight.lines.angle = -1, max.waypoints.distance = 2000,
                              max.flight.time = 15, starting.point = 1) {
  # Check parameters
  if (class(roi)[1] != "SpatialPolygonsDataFrame")
    stop("roi is not a valid polygon layer")
  if (!grep("units=m", as.character(roi@proj4string@projargs)))
    stop("roi is not in a metric projection")

  # Parameters calculated from UAV
  params = flight.parameters(
    sensor.width = sensor.width,
    focal.length35 = focal.length35,
    aspect.ratio = aspect.ratio,
    image.width.px = image.width.px,
    gsd = gsd,
    side.overlap = side.overlap,
    front.overlap = front.overlap,
    flight.speed.kmh = flight.speed.kmh,
    height = height
  )

  flight.speed.kmh = params$flight.speed.kmh
  flightSpeedMs = flight.speed.kmh / 3.6
  height = params$height
  groundHeight = params$ground.height
  groundHeightOverlap = groundHeight*front.overlap
  flightLineDistance = params$flight.line.distance
  vertices = roi@polygons[[1]]@Polygons[[1]]@coords

  # Get bounding box parameters
  if (flight.lines.angle != -1) {
    minBbox = getBBoxAngle(vertices, flight.lines.angle)
  } else {
    # if angle not specified use minimum possible bounding box
    minBbox = getMinBBox(vertices)
  }
  width = minBbox$width
  height = minBbox$height
  alpha = minBbox$angle
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
  # Interleave one upper, two bottom
  # two upper, two bottom... until end
  yHeights = c(rep(c(yHeights, rev(yHeights)), nLines/2+1))
  yHeights = yHeights[1:length(xWidths)]

  # Switch position of the first point
  if (starting.point == 2) {
    yHeights = rev(yHeights)
  } else if (starting.point == 3) {
    xWidths = rev(xWidths)
    yHeights = rev(yHeights)
  } else if (starting.point == 4) {
    xWidths = rev(xWidths)
  }


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
  wktLines = paste(apply(waypoints, 1, paste, collapse=" "), collapse=", ")
  wktLines = paste("LINESTRING(", wktLines,")")
  gLines = rgeos::readWKT(wktLines, p4s = roi@proj4string)
  inter = rgeos::gIntersection(rgeos::gBuffer(roi, width = flightLineDistance), gLines)
  nLines = length(inter@lines[[1]]@Lines)

  flightLines = t(sapply(inter@lines[[1]]@Lines, function(x) x@coords))
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
    hasCurve = !anyNA(curve)
    if (hasCurve) {
      if (curve$before) {
        wptsMatrix[mat_pos,] = c(curve[,1:2], TRUE, FALSE)
        mat_pos = mat_pos + 1
        wptsMatrix[mat_pos,] = c(waypoints[i, 1:2], FALSE, i %% 2 == 1)
        mat_pos = mat_pos + 1
      } else {
        wptsMatrix[mat_pos,] = c(waypoints[i, 1:2], FALSE, i %% 2 == 1)
        mat_pos = mat_pos + 1
        wptsMatrix[mat_pos,] = cbind(curve[,1:2], TRUE, FALSE)
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
  graphics::plot(waypoints[,1:2])
  graphics::polygon(roi@polygons[[1]]@Polygons[[1]]@coords)
  graphics::lines(waypoints[,1:2], lty=2)
  graphics::text(waypoints[,1], waypoints[,2], seq_along(waypoints[,1]), pos=3)


  # Calculate heading
  nWaypoints = nrow(waypoints)
  latDiff = lats[-1]-lats[-nWaypoints]
  lngDiff = lngs[-1]-lngs[-nWaypoints]
  headDegree = atan(latDiff/lngDiff)/pi*180
  finalHeading = 270-headDegree
  finalHeading[lngDiff > 0] = 90-headDegree[lngDiff > 0]

  # Set parameters of the flight in the CSV
  dfLitchi = flightplanning::litchi
  dfLitchi = dfLitchi[rep(1, length(lats)),]
  dfLitchi$latitude = lats
  dfLitchi$longitude = lngs
  dfLitchi$altitude.m. = params$height
  dfLitchi$speed.m.s. = flightSpeedMs
  dfLitchi$heading.deg. = c(finalHeading, 90)
  dfLitchi$curvesize.m. = 0
  dfLitchi$curvesize.m.[waypoints$isCurve==1] = flightLineDistance*0.5
  dfLitchi$photo_timeinterval[waypoints$takePhoto==1] = params$photo.interval
  dfLitchi$actiontype1[waypoints$takePhoto==1] = 1
  dfLitchi$actiontype1[waypoints$takePhoto==0 & waypoints$isCurve==0] = 1
  dfLitchi$gimbal.pitch.angle = gimbal.pitch.angle


  # Split the flight if is too long
  dists = sqrt(diff(waypoints[,1])**2+diff(waypoints[,2])**2)
  distAcum = c(0,cumsum(dists))
  flightTime = distAcum / (flightSpeedMs*0.75) / 60
  finalSize = nrow(dfLitchi)
  totalFlightTime = flightTime[finalSize]
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
    message("Your flight was splitted in ", length(splits), "splits,
because the total time would be ", totalFlightTime, " minutes.")
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

  cat("#####################\n")
  cat("## Flight settings ## \n")
  cat("#####################\n")
  cat("Min shutter speed: ")
  cat(params$minimum.shutter.speed)
  cat("\nPhoto interval:    ")
  cat(params$photo.interval)
  cat(" s")
  cat("\nFlight speed:      ")
  cat(params$flight.speed.kmh)
  cat(" km/h")
  cat("\nFlight lines angle: ")
  cat(alpha)
  cat('\nTotal flight time: ')
  cat(totalFlightTime)
  cat("\n")
}


p3 = list(sensorWD = 6.17,
          focal.length35 = 20,
          aspect.ratio = "4:3",
          image.width.px = 4000)

p4adv = list(sensorWD = 13.2,
          focal.length35 = 24,
          aspect.ratio = "3:2",
          image.width.px = 5472)


#' Function to calculate flight parameters
#'
#' This function will calculate the flight parameters by providing the camera settings
#' target flight height or gsd, front and side overlap.
#'
#' @rdname flight.parameters
#'
#' @param sensor.width numeric. Camera sensor width in milimeters, default 6.17
#' @param focal.length35 numeric. Camera focal length 35mm equivalent, default 20
#' @param aspect.ratio character. Aspect ratio of the picture, default "4:3"
#' @param image.width.px numeric. Width of the image in number of pixels, default 4000
#' @param gsd target ground resolution in centimeters, must provide either `gsd` or `height`
#' @param height target flight height, default NA
#' @param side.overlap desired width overlap between photos
#' @param front.overlap desired height overlap between photos
#' @param flight.speed.kmh flight speed in km/h
#'
#' @examples
#' params = flight.parameters(
#'   gsd=4.325,
#'   flight.speed.kmh=30,
#'   side.overlap = 0.8,
#'   front.overlap = 0.8
#'  )
#'
#' @export
flight.parameters = function(
  sensor.width = 6.17,
  focal.length35 = 20,
  aspect.ratio = "4:3",
  image.width.px = 4000,
  gsd = NA,
  side.overlap = 0.8,
  front.overlap = 0.8,
  flight.speed.kmh = NA,
  height = NA) {

  if (is.na(gsd) == is.na(height)) {
    stop("You must specify either gsd or height!")
  }

  # Size factor to divide
  ratio = 3/4
  sizeFactor = 34.6
  if (aspect.ratio == "3:2") {
    sizeFactor = 36.0
    ratio = 2/3
  }

  realFocalLength = (sensor.width * focal.length35) / sizeFactor
  if (is.na(gsd)) {
    gsd = height * sensor.width*100 / realFocalLength / image.width.px
    groundWidth = image.width.px * gsd / 100
  } else {
    groundWidth = image.width.px * gsd / 100
    height = (groundWidth / sensor.width) * realFocalLength
  }


  flightLineDistance = groundWidth - side.overlap * groundWidth

  flightSpeedMs = flight.speed.kmh / 3.6
  speedPxPerSecond = flightSpeedMs / (gsd*0.01)

  # FIGUEIREDO, E. O. et al.
  # Planos de Voo Semiautônomos para Fotogrametria
  # com Aeronaves Remotamente Pilotadas de Classe 3
  maxPixelRoll = 1.2
  minimumShutterSpeed = paste("1/",round(speedPxPerSecond/maxPixelRoll), sep="")

  groundHeight = groundWidth * ratio
  groundHeightOverlap = groundHeight * front.overlap
  groundAllowedOffset = groundHeight - groundHeightOverlap
  photoInterval = groundAllowedOffset / flightSpeedMs
  if (photoInterval < MIN_PHOTO_INTERVAL) {
    photoInterval = 2
    flightSpeedMs = groundAllowedOffset / photoInterval
    flight.speed.kmh = flightSpeedMs*3.6
    warning(paste0("Speed had to be lowered because frequency of photos would be too high
  New speed: ", flight.speed.kmh, "km/h"))
  } else if ((photoInterval %% 1) > 1e-4) {
    photoInterval = ceiling(photoInterval)
    flightSpeedMs = groundAllowedOffset / photoInterval
    flight.speed.kmh = flightSpeedMs*3.6
    cat(paste0("Speed lowered to ", flight.speed.kmh, "km/h to round up photo interval time\n"))
  }

  return (list(
    height = height,
    gsd = gsd,
    flight.line.distance=flightLineDistance,
    minimum.shutter.speed=minimumShutterSpeed,
    photo.interval=photoInterval,
    ground.height = groundHeight,
    flight.speed.kmh = flight.speed.kmh))
}

# # Example
# #
# (params = flight.parameters(
#   gsd=5,
#   flight.speed.kmh=52,
#   side.overlap = 0.8, #mudar para overlapFront
#   front.overlap = 0.8 #mudar para overpSide
# ))


#TODO
#Make DOUBLE GRID (perpendicular)
