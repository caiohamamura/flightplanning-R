flightplanning-R
================================
[![Build Status](https://travis-ci.com/caiohamamura/flightplanning-R.svg)](https://travis-ci.com/caiohamamura/flightplanning-R)
[![codecov](https://codecov.io/gh/caiohamamura/flightplanning-R/branch/master/graph/badge.svg)](https://codecov.io/gh/caiohamamura/flightplanning-R)
![license](https://img.shields.io/badge/license-MIT-green.svg) 

An R package for generating UAV flight plans, specially for Litchi.


## Installation

This package should be installed using the devtools.

```r
# install.packages("devtools")
devtools::install_github("caiohamamura/flightplanning-R")
```

## Usage
There are two main functions available:
 * `flight.parameters()`: this will calculate the flight parameters given desired settings for GSD/height, target overlap, flight speed and camera specifications.
 * `litchi.plan()`: it depends on the `flight.parameters()` return object to generate the CSV flight plan ready to import into the Litchi Hub.
 
### flight.parameters
 - `height`: target flight height, default NA
 - `gsd`: target ground resolution in centimeters, must provide either `gsd` or `height`
 - `focal.length35`: numeric. Camera focal length 35mm equivalent, default 20
 - `image.width.px`: numeric. Image width in pixels, default 4000
 - `image.height.px`: numeric. Image height in pixels, default 3000
 - `side.overlap`: desired width overlap between photos, default 0.8
 - `front.overlap`: desired height overlap between photos, default 0.8
 - `flight.speed.kmh`: flight speed in km/h, default 54.
 
 ### litchi.plan
  - `roi`: range of interest loaded as an OGR layer, must be in
a metric units projection for working properly
 - `output`: output path for the csv file
 - `flight.params`: Flight Parameters. parameters calculated from flight.parameters()
 - `gimbal.pitch.angle`: gimbal angle for taking photos, default -90 (overriden at flight time)
 - `flight.lines.angle`: angle for the flight lines, default -1 (auto set based on larger direction)
 - `max.waypoints.distance`: maximum distance between waypoints in meters,
default 2000 (some issues have been reported with distances > 2 Km)
 - `max.flight.time`: maximum flight time. If mission is greater than the estimated time, 
 it will be splitted into smaller missions.
 - `starting.point`: numeric (1, 2, 3 or 4). Change position from which to start the flight, default 1
 
## Authors
 - Caio Hamamura
 - Danilo Roberti Alves de Almeida
 - Daniel de Almeida Papa
 - Hudson Franklin Pessoa Veras
 - Evandro Orfanó Figueiredo

This package was developed by author and its contributors which helped providing the calculations and testing.

## Example
``` R
# Install and load the package
install.packages("devtools")
require(devtools)
install_github("caiohamamura/flightplanning-R", "v0.7.1")
library(flightplanning)

params = flight.parameters(height=100,
                          flight.speed.kmh=54,
                          side.overlap = 0.8,
                          front.overlap = 0.8)
                          
params

## Slot "flight.line.distance": 
## [1] 34.61329 
##  
## Slot "flight.speed.kmh": 
## [1] 46.72794 
## 
## Slot "front.overlap":
## [1] 0.8
## 
## Slot "gsd":
## [1] 4.326662
## 
## Slot "height":
## [1] 100
## 
## Slot "ground.height":
## [1] 129.7998
## 
## Slot "minimum.shutter.speed":
## [1] "1/289"
## 
## Slot "photo.interval":
## [1] 2

# Load example SpatialDataFrame polygon
data(exampleBoundary)

# Set the output
output = "output.csv"

# Create the csv plan 
litchi.plan(exampleBoundary,
            output,
            params,
            flight.lines.angle = -1,
            max.waypoints.distance = 2000,
            max.flight.time = 15)
            
# Your flight was splitted in 2splits,
# because the total time would be 27.6237815915755 minutes.
# They were saved as:
# output1.csv
# output2.csv
# The entire flight plan was saved as:
# output_entire.csv
# #####################
# ## Flight settings ## 
# #####################
# Min shutter speed: 1/289
# Photo interval:    2 s
# Flight speed:      46.7279 km/h
# Flight lines angle: 39.0626
# Total flight time: 27.6239
```
<img src="https://github.com/caiohamamura/flightplanning-R/blob/master/man/images/plot_flightplan.png" alt="Flight plan plot" align="center"/>

## References
FIGUEIREDO, E. O. et al. Planos de Voo Semiautônomos para Fotogrametria com Aeronaves Remotamente Pilotadas de Classe 3. Acre: EMBRAPA. 2018
