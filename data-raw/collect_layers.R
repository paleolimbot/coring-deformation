
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
create_quadratic_model <- function(df) {
  model <- lm(data=df, formula=d ~ poly(r, 2, raw=TRUE))
  return(data.frame(a=model$coefficients[3],
                    r2=summary(model)$r.squared,
                    df=model$df.residual))
}

# calculate quadratic models
modelparams <- deformations %>% 
  group_by(photo, layercode) %>%
  do(create_quadratic_model(.))

pcounts <- modelparams %>% group_by(photo) %>% summarise(nlayers=length(layercode))
deformed_core_photos <- merge(deformed_core_photos, pcounts, by="photo")

deformed_layer_data <- deformations; rm(deformations)
deformed_layers <- modelparams; rm(layers, modelparams, pcounts)

rm(`%do%`, create_quadratic_model)
