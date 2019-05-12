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
#' @param GSD target ground resolution in centimeters, default 4
#' @param flightSpeedKmH flight speed in km/h, default 30
#' @param overlapWidth desired width overlap between photos, default 0.8
#' @param overlapHeight desired height overlap between photos, defualt 0.8
#' @param gimbalPitchAngle gimbal angle for taking photos, default -90 (can be overriden at flight time)
#' @param flightLinesAngle angle for the flight lines, default -1 (auto set based on larger dimensions)
#' @param maxWaypointsDistance maximum distance between waypoints in meters,
#' default 2000 (some issues have been reported with distances > 2km)
#'
#' @examples
#' data(exampleBoundary)
#' outPath = tempfile(fileext=".csv")
#' generateLitchiPlan(ogrROI = exampleBoundary,
#'                    outputPath = outPath,
#'                    uav = "p3",
#'                    GSD = 4.325,
#'                    flightSpeedKmH = 30,
#'                    overlapWidth = 0.8,
#'                    overlapHeight = 0.8,
#'                    gimbalPitchAngle = -90,
#'                    flightLinesAngle = -1,
#'                    maxWaypointsDistance = 2000)
#'
#'
#' @export
#' @import sp rgeos rgdal
#' @importFrom graphics text
#' @importFrom utils data read.csv write.csv
generateLitchiPlan = function(ogrROI, outputPath, uav = "p3",
                              GSD = 4, flightSpeedKmH = 30,
                              overlapWidth = 0.8, overlapHeight = 0.8,
                              gimbalPitchAngle = -90, flightLinesAngle = -1,
                              maxWaypointsDistance = 2000) {
  if (summary(ogrROI)[2] != "SpatialPolygonsDataFrame")
    stop("ogrROI is not a valid polygon layer")
  if (!grep("units=m", as.character(ogrROI@proj4string@projargs)))
    stop("ogrROI is not in a metric projection")
  params = flightParameters(
    uav = uav,
    GSD = GSD,
    overlapWidth = overlapWidth,
    overlapHeight = overlapHeight,
    flightSpeedKmH = flightSpeedKmH
  )
  flightSpeedMs = flightSpeedKmH / 3.6
  groundHeight = params$groundHeight
  groundHeightOverlap = groundHeight*overlapHeight
  flightLineDistance = params$flightLineDistance

  vertices = ogrROI@polygons[[1]]@Polygons[[1]]@coords
  minBbox = getMinBBox(vertices)
  width = minBbox$width
  height = minBbox$height
  alpha = minBbox$angle
  if (flightLinesAngle != -1) {
    alpha = flightLinesAngle
  }
  rads=alpha*pi/180

  # install.packages("rgeos")
  # library(rgeos)
  pol = sp::Polygon(rbind(minBbox$pts, minBbox$pts[1,]))@coords
  polDef=paste(apply(pol, 1, paste, collapse=" "), collapse=", ")
  minBboxGeos = rgeos::readWKT(paste("POLYGON((", polDef, "))", sep=""))
  centroid = rgeos::gCentroid(minBboxGeos)@coords

  if (flightLinesAngle != -1) {
    height = width
  }
  nLines = ceiling(height / flightLineDistance)+1
  xWidths = (-nLines/2):(nLines/2) * flightLineDistance
  xWidths = rep(xWidths, each=2)

  heightDistance = groundHeight-groundHeightOverlap
  heightAdjusted = width + 2*heightDistance
  heightMHalf = -heightAdjusted/2
  heightPHalf = heightAdjusted/2
  yHeights = c(heightMHalf, heightPHalf)
  yHeights = c(rep(c(yHeights, rev(yHeights)), nLines/2+1))
  yHeights = yHeights[1:length(xWidths)]

  xys = data.frame(x=-xWidths*sin(rads)+yHeights*cos(rads), y=xWidths*cos(rads)+yHeights*sin(rads))

  # Initial waypoints to intersect waypoints
  waypoints = xys + rep(centroid, each=nrow(xys))

  # For some reason gIntersection with MULTILINESTRING
  # returns in inconsistent order though needs to be
  # done in a for loop
  lineList = matrix(ncol=2, nrow=0)
  for (i in 1:(nLines+1)) {
    pt = i*2-1
    lineCoords = waypoints[pt:(pt+1),]
    Li = rgeos::readWKT(paste("LINESTRING(",paste(apply(lineCoords, 1, paste, collapse=" "), collapse=", "), ")", sep=""), p4s = ogrROI@proj4string)
    inter = rgeos::gIntersection(ogrROI, Li)
    if (!is.null(inter)) {
      for (l in inter@lines[[1]]@Lines) {
        l2 = sp::SpatialLines(list(sp::Lines(list(l), ID=1)))
        LiLength = rgeos::gLength(l2)
        if (LiLength > maxWaypointsDistance) {
          splitNum = ceiling(LiLength/maxWaypointsDistance)+1
          splitLocations = seq(0, LiLength, length.out = splitNum)
          interpolates = rgeos::gInterpolate(l2, splitLocations)
          inter = rgeos::readWKT(paste("LINESTRING(", paste(apply(interpolates@coords, 1, paste, collapse=" "), collapse=", "), ")", sep=""))
          for (l3 in inter@lines[[1]]@Lines) {
            lineList = append(lineList, t(l3@coords))
          }
        } else {
          lineList = append(lineList, t(l@coords))
        }
      }
    }
  }
  waypoints = matrix(lineList, ncol=2, byrow=T)


  # Waypoints
  wgs84 = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  transform = rgdal::rawTransform(ogrROI@proj4string@projargs, wgs84, n=nrow(waypoints), x=waypoints[,1], y=waypoints[,2])
  lats = transform[[2]]
  lngs = transform[[1]]
  sp::plot(ogrROI)
  graphics::lines(waypoints, lty=2)
  graphics::points(waypoints)
  graphics::text(waypoints[,1], waypoints[,2], seq_len(length(lats)), pos=3)


  # Calculate heading
  nWaypoints = nrow(waypoints)
  latDiff = lats[-1]-lats[-nWaypoints]
  lngDiff = lngs[-1]-lngs[-nWaypoints]
  headDegree = atan(latDiff/lngDiff)/pi*180
  finalHeading = 270-headDegree
  finalHeading[lngDiff > 0] = 90-headDegree[lngDiff > 0]

  #
  # dfLitchi=read.csv("litchi.csv", header=TRUE)
  # dfLitchi$gimbalpitchangle = -90
  # write.csv(dfLitchi, "litchi.csv", row.names=FALSE)
  dfLitchi = flightplanning::litchi
  dfLitchi = dfLitchi[rep(1, length(lats)),]
  dfLitchi$latitude = lats
  dfLitchi$longitude = lngs
  dfLitchi$altitude.m. = params$flightHeight
  dfLitchi$speed.m.s. = flightSpeedMs
  dfLitchi$heading.deg. = c(finalHeading, 90)
  dfLitchi$curvesize.m. = 0
  dfLitchi$photo_timeinterval[seq(1,nWaypoints, 2)] = params$photoInterval
  dfLitchi$gimbalpitchangle = gimbalPitchAngle
  write.csv(dfLitchi, outputPath, row.names = FALSE)
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
#' @param uav either "p3" or "p4adv" for loading Phantom 3-4std or Phanton4-adv/pro
#' camera profiles
#' @param GSD target ground resolution in centimeters
#' @param flightSpeedKmH flight speed in km/h
#' @param overlapWidth desired width overlap between photos
#' @param overlapHeight desired height overlap between photos
#'
#' @examples
#' params = flightParameters(
#'   GSD=4.325,
#'   flightSpeedKmH=30,
#'   overlapWidth = 0.8,
#'   overlapHeight = 0.8
#'  )
#'
#' @export
flightParameters = function(
  uav="p3",
  GSD = 4,
  overlapWidth = 0.8,
  overlapHeight = 0.8,
  flightSpeedKmH = 30) {

  uavModel = p3

  if (uav == "p4adv") {
    uavModel = p4adv
  }

  # Size factor to divide
  ratio = 3/4
  sizeFactor = 34.6
  if (uavModel$aspectRatio == "3:2") {
    sizeFactor = 36.0
    ratio = 2/3
  }

  realFocalLength = (uavModel$sensorWD * uavModel$focalLength35) / sizeFactor
  groundWidth = uavModel$imageWidthPx * GSD / 100

  flightHeight = (groundWidth / uavModel$sensorWD) * realFocalLength
  flightLineDistance = groundWidth - overlapWidth * groundWidth

  flightSpeedMs = flightSpeedKmH / 3.6
  speedPxPerSecond = flightSpeedMs / (GSD*0.01)

  # FIGUEIREDO, E. O. et al.
  # Planos de Voo Semiautônomos para Fotogrametria
  # com Aeronaves Remotamente Pilotadas de Classe 3
  maxPixelRoll = 1.2
  minimumShutterSpeed = paste("1:",round(speedPxPerSecond/maxPixelRoll), sep="")

  groundHeight = groundWidth * ratio
  groundHeightOverlap = groundHeight * overlapHeight
  groundAllowedOffset = groundHeight - groundHeightOverlap
  photoInterval = groundAllowedOffset / flightSpeedMs

  return (list(
    flightHeight=flightHeight,
    flightLineDistance=flightLineDistance,
    minimumShutterSpeed=minimumShutterSpeed,
    photoInterval=photoInterval,
    groundHeight = groundHeight))
}

# Example
#
# (params = flightParameters(
#   GSD=4.325,
#   flightSpeedKmH=30,
#   overlapWidth = 0.8,
#   overlapHeight = 0.8
# ))
