MIN_PHOTO_INTERVAL = 2
DIAG_35MM = sqrt(36^2 + 24^2) # Classical 35mm film diagonal
MAX_WAYPOINTS = 99



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

  if ((is.na(gsd) & is.na(height)) || (!is.na(gsd) & !is.na(height)))  {
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
