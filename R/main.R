MIN_PHOTO_INTERVAL = 2

#'  Function to generate Litchi csv flight plan
#'
#'
#' @rdname generateLitchiPlan
#'
#' @param ogrROI range of interest loaded as an OGR layer, must be in
#' a metric units projection for working properly
#' @param outputPath output path for the csv file
#' @param uav either "p3" or "p4adv" for loading Phantom 3-4std or Phanton4-adv/pro
#' camera profiles, default "p3"
#' @param GSD target ground resolution in centimeters, must specify either GSD or flightHeight, default NA
#' @param flightHeight flight height in meters, default NA
#' @param flightSpeedKmH flight speed in km/h, default 54
#' @param sideOverlap desired width overlap between photos, default 0.8
#' @param frontOverlap desired height overlap between photos, defualt 0.8
#' @param gimbalPitchAngle gimbal angle for taking photos, default -90 (can be overriden at flight time)
#' @param flightLinesAngle angle for the flight lines, default -1 (auto set based on larger dimensions)
#' @param maxWaypointsDistance maximum distance between waypoints in meters,
#' default 2000 (some issues have been reported with distances > 2km)
#' @param maxFlightTime maximum flight time. If mission is greater than the estimated
#' time, it will be splitted into smaller missions.
#' @param startingPoint numeric (1, 2, 3 or 4). Change position from which to start the flight, default 1
#'
#' @examples
#' data(exampleBoundary)
#' outPath = tempfile(fileext=".csv")
#' generateLitchiPlan(ogrROI = exampleBoundary,
#'                    outputPath = outPath,
#'                    uav = "p3",
#'                    GSD = 4.325,
#'                    flightSpeedKmH = 30,
#'                    sideOverlap = 0.8,
#'                    frontOverlap = 0.8,
#'                    gimbalPitchAngle = -90,
#'                    flightLinesAngle = -1,
#'                    maxWaypointsDistance = 2000,
#'                    maxFlightTime = 15)
#'
#'
#' @export
#' @import sp rgeos rgdal
#' @importFrom graphics text
#' @importFrom utils data read.csv write.csv
generateLitchiPlan = function(ogrROI, outputPath, sensorWidth = 6.17,
                              focalLength35=20, aspectRatio = "4:3",
                              imageWidthPx = 4000, GSD = NA, flightHeight = NA,
                              flightSpeedKmH = 54, sideOverlap = 0.8,
                              frontOverlap = 0.8, gimbalPitchAngle = -90,
                              flightLinesAngle = -1, maxWaypointsDistance = 2000,
                              maxFlightTime = 15, startingPoint = 1) {
  # Check parameters
  if (class(ogrROI)[1] != "SpatialPolygonsDataFrame")
    stop("ogrROI is not a valid polygon layer")
  if (!grep("units=m", as.character(ogrROI@proj4string@projargs)))
    stop("ogrROI is not in a metric projection")

  # Parameters calculated from UAV
  params = flightParameters(
    sensorWidth = sensorWidth,
    focalLength35 = focalLength35,
    aspectRatio = aspectRatio,
    imageWidthPx = imageWidthPx,
    GSD = GSD,
    sideOverlap = sideOverlap,
    frontOverlap = frontOverlap,
    flightSpeedKmH = flightSpeedKmH,
    flightHeight = flightHeight
  )

  flightSpeedKmH = params$flightSpeedKmH
  flightSpeedMs = flightSpeedKmH / 3.6
  flightHeight = params$flightHeight
  groundHeight = params$groundHeight
  groundHeightOverlap = groundHeight*frontOverlap
  flightLineDistance = params$flightLineDistance
  vertices = ogrROI@polygons[[1]]@Polygons[[1]]@coords

  # Get bounding box parameters
  if (flightLinesAngle != -1) {
    minBbox = getBBoxAngle(vertices, flightLinesAngle)
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
  if (startingPoint == 2) {
    yHeights = rev(yHeights)
  } else if (startingPoint == 3) {
    xWidths = rev(xWidths)
    yHeights = rev(yHeights)
  } else if (startingPoint == 4) {
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
  gLines = rgeos::readWKT(wktLines, p4s = ogrROI@proj4string)
  inter = rgeos::gIntersection(rgeos::gBuffer(ogrROI, width = flightLineDistance), gLines)
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

  # Transform to WGS84 latitude and longitude
  transform = rgdal::rawTransform(ogrROI@proj4string@projargs, wgs84, n=nrow(waypoints), x=waypoints[,1], y=waypoints[,2])
  lats = transform[[2]]
  lngs = transform[[1]]
  graphics::plot(waypoints[,1:2])
  graphics::polygon(ogrROI@polygons[[1]]@Polygons[[1]]@coords)
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
  dfLitchi$altitude.m. = params$flightHeight
  dfLitchi$speed.m.s. = flightSpeedMs
  dfLitchi$heading.deg. = c(finalHeading, 90)
  dfLitchi$curvesize.m. = 0
  dfLitchi$curvesize.m.[waypoints$isCurve==1] = flightLineDistance*0.5
  dfLitchi$photo_timeinterval[waypoints$takePhoto==1] = params$photoInterval
  dfLitchi$actiontype1[waypoints$takePhoto==1] = 1
  dfLitchi$actiontype1[waypoints$takePhoto==0 & waypoints$isCurve==0] = 1
  dfLitchi$gimbalpitchangle = gimbalPitchAngle


  # Split the flight if is too long
  dists = sqrt(diff(waypoints[,1])**2+diff(waypoints[,2])**2)
  distAcum = c(0,cumsum(dists))
  flightTime = distAcum / (flightSpeedMs*0.75) / 60
  finalSize = nrow(dfLitchi)
  totalFlightTime = flightTime[finalSize]
  if (totalFlightTime > maxFlightTime) {
    indexes = seq_len(finalSize)
    nBreaks = ceiling(totalFlightTime/maxFlightTime)
    breaks = seq(0, flightTime[finalSize], length.out = nBreaks+1)[c(-1, -nBreaks-1)]
    endWaypointsIndex = indexes[waypoints$isCurve & (seq_len(finalSize) %% 2 == 0)]
    endWaypoints = flightTime[waypoints$isCurve & (seq_len(finalSize) %% 2 == 0)]
    selected = sapply(breaks, function(x) which.min(abs(endWaypoints-x)))
    waypointsBreak = endWaypointsIndex[indexes[selected]]


    dfLitchi$split = rep(1:nBreaks, diff(c(0, waypointsBreak, finalSize)))
    splits=split.data.frame(waypoints, f = dfLitchi$split)
    for (dataSplit in split.data.frame(dfLitchi, f = dfLitchi$split)) {
      first = substr(outputPath, 1, nchar(outputPath)-4)
      second = substr(outputPath, nchar(outputPath)-3, nchar(outputPath))
      outputPath2 = paste(first, dataSplit[1, ]$split, second, sep = "")
      write.csv(dataSplit[,-ncol(dataSplit)], outputPath2, row.names = FALSE)
    }
  } else {
    write.csv(dfLitchi, outputPath, row.names = FALSE)
  }

  cat("#####################\n")
  cat("## Flight settings ## \n")
  cat("#####################\n")
  cat("Min shutter speed: ")
  cat(params$minimumShutterSpeed)
  cat("\nPhoto interval:    ")
  cat(params$photoInterval)
  cat(" s")
  cat("\nFlight speed:      ")
  cat(params$flightSpeedKmH)
  cat(" km/h")
  cat("\nFlight line angle: ")
  cat(alpha)
  cat('\n')
}


p3 = list(sensorWD = 6.17,
          focalLength35 = 20,
          aspectRatio = "4:3",
          imageWidthPx = 4000)

p4adv = list(sensorWD = 13.2,
          focalLength35 = 24,
          aspectRatio = "3:2",
          imageWidthPx = 5472)


#'  Function to generate Litchi csv flight plan
#'
#'
#' @rdname flightParameters
#'
#' @param sensorWidth numeric. Camera sensor width in milimeters, default 6.17
#' @param focalLength35 numeric. Camera focal length 35mm equivalent, default 20
#' @param aspectRatio character. Aspect ratio of the picture, default "4:3"
#' @param imageWidth numeric. Width of the image in number of pixels, default 4000
#' @param GSD target ground resolution in centimeters
#' @param flightSpeedKmH flight speed in km/h
#' @param sideOverlap desired width overlap between photos
#' @param frontOverlap desired height overlap between photos
#'
#' @examples
#' params = flightParameters(
#'   GSD=4.325,
#'   flightSpeedKmH=30,
#'   sideOverlap = 0.8,
#'   frontOverlap = 0.8
#'  )
#'
#' @export
flightParameters = function(
  sensorWidth = 6.17,
  focalLength35 = 20,
  aspectRatio = "4:3",
  imageWidthPx = 4000,
  GSD = 4,
  sideOverlap = 0.8,
  frontOverlap = 0.8,
  flightSpeedKmH = NA,
  flightHeight = NA) {

  if (is.na(GSD) && is.na(flightHeight)) {
    stop("You must specify either GSD or flightHeight!")
  }

  # Size factor to divide
  ratio = 3/4
  sizeFactor = 34.6
  if (aspectRatio == "3:2") {
    sizeFactor = 36.0
    ratio = 2/3
  }

  realFocalLength = (sensorWidth * focalLength35) / sizeFactor
  if (is.na(GSD)) {
    GSD = flightHeight * sensorWidth*100 / realFocalLength / imageWidthPx
    groundWidth = imageWidthPx * GSD / 100
  } else {
    groundWidth = imageWidthPx * GSD / 100
    flightHeight = (groundWidth / sensorWidth) * realFocalLength
  }


  flightLineDistance = groundWidth - sideOverlap * groundWidth

  flightSpeedMs = flightSpeedKmH / 3.6
  speedPxPerSecond = flightSpeedMs / (GSD*0.01)

  # FIGUEIREDO, E. O. et al.
  # Planos de Voo Semiautônomos para Fotogrametria
  # com Aeronaves Remotamente Pilotadas de Classe 3
  maxPixelRoll = 1.2
  minimumShutterSpeed = paste("1/",round(speedPxPerSecond/maxPixelRoll), sep="")

  groundHeight = groundWidth * ratio
  groundHeightOverlap = groundHeight * frontOverlap
  groundAllowedOffset = groundHeight - groundHeightOverlap
  photoInterval = groundAllowedOffset / flightSpeedMs
  if (photoInterval < MIN_PHOTO_INTERVAL) {
    photoInterval = 2
    flightSpeedMs = groundAllowedOffset / photoInterval
    flightSpeedKmH = flightSpeedMs*3.6
    warning(paste0("Speed had to be lowered because frequency of photos would be too high
  New speed: ", flightSpeedKmH, "km/h"))
  } else if ((photoInterval %% 1) > 1e-4) {
    photoInterval = ceiling(photoInterval)
    flightSpeedMs = groundAllowedOffset / photoInterval
    flightSpeedKmH = flightSpeedMs*3.6
    cat(paste0("Speed lowered to ", flightSpeedKmH, "km/h to round up photo interval time\n"))
  }

  return (list(
    flightHeight = flightHeight,
    flightLineDistance=flightLineDistance,
    minimumShutterSpeed=minimumShutterSpeed,
    photoInterval=photoInterval,
    groundHeight = groundHeight,
    flightSpeedKmH = flightSpeedKmH))
}

# # Example
# #
# (params = flightParameters(
#   GSD=4.325,
#   flightSpeedKmH=30,
#   sideOverlap = 0.8, #mudar para overlapFront
#   frontOverlap = 0.8 #mudar para overpSide
# ))


#TODO
#Possibilidade de DOUBLE GRID (perpendicular ao mais comprido)
