library(flightplanning)
f <- "mytest/uav_bug.gpkg"
roi <- sf::st_read(f, layer = "Zalew")
str(roi)
if(nrow(roi) > 1) {
  roi <- sf::st_union(roi)
}

output <- "mytest/fly.csv"

params <- flight.parameters(
  height = 120,
  focal.length35 = 24,
  flight.speed.kmh = 24,
  side.overlap = 0.7,
  front.overlap = 0.7
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

litchi_sf(roi,
          output,
          params,
          gimbal.pitch.angle = -90,
          flight.lines.angle = -1,
          max.waypoints.distance = 400,
          max.flight.time = 18,
          grid = TRUE
)



litchi.plan(roi,
          output,
          params,
          gimbal.pitch.angle = -90,
          flight.lines.angle = -1,
          max.waypoints.distance = 400,
          max.flight.time = 18,
          grid = FALSE
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



# Create the csv plan
flightplanning::litchi.plan(roi,
                            output,
                            params,
                            gimbal.pitch.angle = -90,
                            flight.lines.angle = -1,
                            max.waypoints.distance = 400,
                            max.flight.time = 18,
                            grid = FALSE)
