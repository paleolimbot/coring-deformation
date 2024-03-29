
`%do%` <- foreach::`%do%`
library(dplyr)

layer_files <- list.files("data-raw/layers", pattern=".txt", full.names=TRUE)

deformations <- foreach::foreach(layer=layer_files, .combine=rbind) %do% {
  df <- read.delim(layer, header=FALSE)
  names(df) <- c("x", "y")
  layer <- basename(layer)
  df$photo <- strsplit(layer, ".", fixed=TRUE)[[1]][1]
  df$layer <- strsplit(layer, ".", fixed=TRUE)[[1]][2]
  df$layercode <- substr(layer, 1, nchar(layer)-4)
  df
}
rm(df, layer, layer_files)

deformed_core_photos <- read.csv("data-raw/photos/photo_info.csv")
deformations <- merge(deformations, deformed_core_photos, by="photo")
deformations$x <- deformations$x / deformations$scale
deformations$y <- deformations$y / deformations$scale
layers <- deformations %>% 
  group_by(layercode) %>% 
  summarise(x0=mean(range(x)), y0=min(y))
deformations <- merge(deformations, layers, by="layercode")
deformations$r <- deformations$x-deformations$x0
deformations$d <- deformations$y-deformations$y0
deformations <- deformations %>% select(layercode, photo, layer, x, y, r, d)

# define the regression function
create_logarithmic_model <- function(df) {
  # from Acton et al. 2002
  Z <- function(r, R, b) -b*(log(1-r/R) + r/R)
  # start with 0, increment by 0.01 until residuals decrease
  resids <- NA
  r <- abs(df$r)
  R <- max(r)
  r <- r[2:(length(r)-1)]
  d <- df$d[2:(length(df$d)-1)]
  for(b in seq(0, 10, 0.01)) {
    vals <- Z(r, b=b, R=R)
    mask <- is.nan(vals) | is.infinite(vals)
    residnew <- sum((vals[!mask] - d[!mask])^2)
    if(is.na(resids) || residnew < resids) {
      resids <- residnew
    } else {
      break
    }
  }
  return(data.frame(b=b, 
                    r2=1 - (resids / (var(d)*nrow(df)))
                    ))
}

# calculate quadratic models
modelparams <- deformations %>% 
  group_by(photo, layercode) %>%
  do(create_logarithmic_model(.))

pcounts <- modelparams %>% group_by(photo) %>% summarise(nlayers=length(layercode))
deformed_core_photos <- merge(deformed_core_photos, pcounts, by="photo")

deformed_layer_data <- deformations; rm(deformations)
deformed_layers <- modelparams; rm(layers, modelparams, pcounts)

rm(`%do%`, create_logarithmic_model)
