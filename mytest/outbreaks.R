#gLines
tmp <- sf::st_sfc()
class(tmp)[1] <- "sf::sfc_LINESTRING"
df <- sf::st_sf(ID = integer(0), geometry = tmp)
roiCRS <- sf::st_crs(roi)
sf::st_crs(df) <- roiCRS
nrows <- nrow(waypoints)-1
for (i in 1:nrows) {
  df[i,1] <- i
  df[i,2] <-sf::st_as_sf(
    sf::st_as_sfc(
      paste("LINESTRING(", waypoints[1,1], waypoints[1,2], ",", waypoints[1+1,1], waypoints[1+1,2], ")\"")
    ),
    crs = roiCRS)
}


# funkcja do zmiany nazwy kolumny z geometrią -------------------------------------------------

rename_geometry <- function(df, name){
  current = attr(df, "sf_column")
  names(df)[names(df) == current] = name
  sf::st_geometry(df) = name
  df
}
