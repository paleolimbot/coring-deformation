
`%do%` <- foreach::`%do%`

layers <- list.files("data-raw/layers", pattern=".csv", full.names=TRUE)
deformed_layers <- foreach::foreach(layer=layers, .combine=rbind) %do% {
  df <- read.csv(layer)
  names(df) <- c("x", "y")
  layer <- basename(layer)
  df$photo <- strsplit(layer, ".", fixed=TRUE)[[1]][1]
  df$layer <- strsplit(layer, ".", fixed=TRUE)[[1]][2]
  df$layercode <- substr(layer, 1, nchar(layer)-4)
  df
}

deformed_core_photos <- read.csv("data-raw/photos/photo_info.csv")

