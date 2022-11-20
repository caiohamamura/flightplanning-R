library(flightplanning)
f <- "mytest/lasek.gpkg"
f <- "/home/sapi/projekty/flightplanning-R/mytest/lasek.gpkg"
roi <- sf::st_read(f)

output <- "mytest/fly.csv"
output <- "/home/sapi/projekty/flightplanning-R/mytest/fly.csv"

params <- flight.parameters(
  height = 120,
  focal.length35 = 24,
  flight.speed.kmh = 24,
  side.overlap = 0.8,
  front.overlap = 0.8
)

litchi_sf(roi,
  output,
  params,
  gimbal.pitch.angle = -90,
  flight.lines.angle = -1,
  max.waypoints.distance = 400,
  max.flight.time = 18,
  grid = FALSE
)

if(nrow(roi) > 1) {
  roi <- sf::st_union(roi)
}
litchi_sf(roi,
          output,
          params,
          gimbal.pitch.angle = -90,
          flight.lines.angle = -1,
          max.waypoints.distance = 400,
          max.flight.time = 18,
          grid = FALSE
)



# Create the csv plan
flightplanning::litchi.plan(roi,
                            output,
                            params,
                            gimbal.pitch.angle = -90,
                            flight.lines.angle = -1,
                            max.waypoints.distance = 400,
                            max.flight.time = 18,
                            grid = FALSE)
