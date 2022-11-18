params = flight.parameters(height = 100)
TOLERANCE = 1e-4

test_that("Flight parameters must have either gsd or height set up", {
  expect_error( flight.parameters(gsd=NA, height=NA) )
  expect_error( flight.parameters(gsd=5, height=100) )
})


test_that("Flight parameters return an S4 Flight Parameters", {
  expect_equal( typeof(params), "S4" )
  expect_equal( is(params), "Flight Parameters" )
})


test_that("GSD calculation from height is correct", {
  params = flight.parameters(height = 100,
                             focal.length35 = 20,
                             image.width.px = 5472,
                             image.height.px = 3648,
                             flight.speed.kmh = 43.2)
  expect_equal( params@gsd, 3.289473684210527, tolerance = TOLERANCE )
})


test_that("Height calculation from GSD is correct", {
  params = flight.parameters(gsd = 5,
                             focal.length35 = 20,
                             image.width.px = 5472,
                             image.height.px = 3648,
                             flight.speed.kmh = 54)
  expect_equal( params@height, 152, tolerance = TOLERANCE )
})

test_that("Flight parameters front overlap is the same as input", {
  front.overlap = 0.7
  params = flight.parameters(height = 100,
                             focal.length35 = 20,
                             image.width.px = 5472,
                             image.height.px = 3648,
                             flight.speed.kmh = 43,
                             front.overlap = front.overlap)
  expect_equal( params@front.overlap, front.overlap, tolerance = TOLERANCE )
})


test_that("Side overlap is correct", {
  side.overlap = 0.5
  gsd = 5
  image.width.px = 5472
  ground.width = image.width.px * (gsd/100)
  overlap.meters = side.overlap * ground.width
  params = flight.parameters(gsd = gsd,
                             focal.length35 = 20,
                             image.width.px = image.width.px,
                             image.height.px = 3648,
                             flight.speed.kmh = 43.2,
                             side.overlap = side.overlap)

  expect_equal( params@flight.line.distance, overlap.meters, tolerance = TOLERANCE )
})


test_that("Front overlap is correct", {
  front.overlap = 0.5
  gsd = 5
  image.height.px = 3648
  ground.height = image.height.px * (gsd/100)
  overlap.meters = front.overlap * ground.height
  speed.kmh = 46.9028571428571
  speed.ms = speed.kmh / 3.6
  interval = overlap.meters / speed.ms

  params = flight.parameters(gsd = gsd,
                             focal.length35 = 20,
                             image.width.px = 5472,
                             image.height.px = image.height.px,
                             flight.speed.kmh = speed.kmh,
                             front.overlap = front.overlap)

  expect_equal( params@photo.interval, interval, tolerance = TOLERANCE )
})


test_that("Photo time interval is rounded up and speed adjusted", {
  front.overlap = 0.5
  gsd = 5
  image.height.px = 3648
  ground.height = image.height.px * (gsd/100)
  overlap.meters = front.overlap * ground.height
  speed.kmh = 54
  speed.ms = speed.kmh / 3.6
  interval = overlap.meters / speed.ms
  rounded.interval = ceiling(interval)
  new.speed.ms = overlap.meters / rounded.interval
  new.speed.kmh = new.speed.ms * 3.6

  params = flight.parameters(gsd = gsd,
                             focal.length35 = 20,
                             image.width.px = 5472,
                             image.height.px = image.height.px,
                             flight.speed.kmh = speed.kmh,
                             front.overlap = front.overlap)

  expect_equal( params@photo.interval, rounded.interval, tolerance = TOLERANCE )
  expect_equal( params@flight.speed.kmh, new.speed.kmh, tolerance = TOLERANCE )
})


test_that("Ground height is properly calculated", {
  gsd = 5
  image.height.px = 3648
  ground.height = image.height.px * (gsd/100)

  params = flight.parameters(gsd = gsd,
                             focal.length35 = 20,
                             image.width.px = 5472,
                             image.height.px = image.height.px,
                             flight.speed.kmh = 54,
                             front.overlap = 0.5)

  expect_equal( params@ground.height, ground.height, tolerance = TOLERANCE )
})


test_that("Shutter speed calculation is correct", {
  gsd = 5
  speed.kmh = 46.9028571428571
  speed.ms = speed.kmh / 3.6
  speed.pxs = speed.ms / (gsd/100)

  MAX_PX_ROLL = 1.2

  time.roll = MAX_PX_ROLL / speed.pxs
  shutter.speed = paste0("1/", round(1/time.roll))

  params = flight.parameters(gsd = gsd,
                             focal.length35 = 20,
                             image.width.px = 5472,
                             image.height.px = 3648,
                             flight.speed.kmh = speed.kmh,
                             front.overlap = 0.5)

  expect_equal( params@minimum.shutter.speed, shutter.speed )
})


test_that("Litchi plan outputs the csv file", {
  exampleBoundary = sf::st_read(
    system.file("extdata", "exampleBoundary.shp", package="flightplanning")) |>
    sf::as_Spatial()

  outPath = tempfile(fileext=".csv")

  params = flight.parameters(
    gsd = 4,
    side.overlap = 0,
    front.overlap = 0,
    flight.speed.kmh = 54
  )

  litchi.plan(exampleBoundary,
              outPath,
              params)
  title("Defaults")

  expect_true(file.exists(outPath))
})


test_that("Different starting points are working", {
  exampleBoundary = sf::st_read(
    system.file("extdata", "exampleBoundary.shp", package="flightplanning")) |>
    sf::as_Spatial()
  outPath = tempfile(fileext=".csv")

  params = flight.parameters(
    gsd = 4,
    side.overlap = 0,
    front.overlap = 0,
    flight.speed.kmh = 54
  )

  litchi.plan(exampleBoundary,
              outPath,
              params,
              starting.point = 2)
  title("Starting point 2")
  litchi.plan(exampleBoundary,
              outPath,
              params,
              starting.point = 3)
  title("Starting point 3")
  litchi.plan(exampleBoundary,
              outPath,
              params,
              starting.point = 4)
  title("Starting point 4")
  succeed()
})


test_that("Different flight line angles are working", {
  exampleBoundary = sf::st_read(
    system.file("extdata", "exampleBoundary.shp", package="flightplanning")) |>
    sf::as_Spatial()
  outPath = tempfile(fileext=".csv")

  params = flight.parameters(
    gsd = 4,
    side.overlap = 0,
    front.overlap = 0,
    flight.speed.kmh = 54
  )

  litchi.plan(exampleBoundary,
              outPath,
              params,
              flight.lines.angle = 45)
  title("45 degrees")
  litchi.plan(exampleBoundary,
              outPath,
              params,
              flight.lines.angle = 90)
  title("90 degrees")
  litchi.plan(exampleBoundary,
              outPath,
              params,
              flight.lines.angle = 135)
  title("135 degrees")
  succeed()
})


test_that("Did not provide legal ROI", {
  outPath = tempfile(fileext=".csv")
  expect_error( litchi.plan(NA, outPath, NA) )
})


test_that("ROI is not in a metric projection", {
  outPath = tempfile(fileext=".csv")
  exampleBoundary = sf::st_read(
    system.file("extdata", "exampleBoundary.shp", package="flightplanning")) |>
    sf::as_Spatial()
  roi = exampleBoundary |>
    sf::st_as_sf() |>
    sf::st_transform(crs = "EPSG:4326") |>
    sf::as_Spatial()
  expect_error( litchi.plan(roi, outPath, NA) )
})


test_that("Did not provide Flight Parameters", {
  exampleBoundary = sf::st_read(
    system.file("extdata", "exampleBoundary.shp", package="flightplanning")) |>
    sf::as_Spatial()
  outPath = tempfile(fileext=".csv")
  expect_error( litchi.plan(exampleBoundary, outPath, NA) )
})


test_that("Break waypoints too far", {
  outPath = tempfile(fileext=".csv")
  exampleBoundary = sf::st_read(
    system.file("extdata", "exampleBoundary.shp", package="flightplanning")) |>
    sf::as_Spatial()
  params = flight.parameters(
    gsd = 4,
    side.overlap = 0,
    front.overlap = 0,
    flight.speed.kmh = 54
  )

  litchi.plan(exampleBoundary, outPath, params,
                            max.waypoints.distance = 1000)
  title("Break waypoints farther than 1000 meters")

  succeed()
})


test_that("Break flight if exceeds max flight time", {
  outPath = tempfile(fileext=".csv")
  exampleBoundary = sf::st_read(
    system.file("extdata", "exampleBoundary.shp", package="flightplanning")) |>
    sf::as_Spatial()
  params = flight.parameters(
    gsd = 4,
    side.overlap = 0,
    front.overlap = 0,
    flight.speed.kmh = 54
  )

  litchi.plan(exampleBoundary, outPath, params,
              max.flight.time = 10)
  title("Break into multiple flights")
  expect_equal(length(Sys.glob(paste0(tools::file_path_sans_ext(outPath), "*.csv"))), 3)
})
