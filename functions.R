haversine <- function(coord1, coord2) {
  R <- 6371.01  # Earth's radius in kilometers
  
  lat1 <- as.numeric(coord1[2]) * pi / 180
  long1 <- as.numeric(coord1[1]) * pi / 180
  lat2 <- as.numeric(coord2[2]) * pi / 180
  long2 <- as.numeric(coord2[1]) * pi / 180
  
  dlat <- lat2 - lat1
  dlong <- long2 - long1
  
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlong/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  
  d <- R * c
  return(d)
}
