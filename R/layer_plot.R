


plot_layer <- function(layercode, zoom=0.6,
                       lwd=3, col="white", lty=2, title=NULL, ...) {
  layer <- deformed_layers[deformed_layers$layercode==layercode,]
  if(nrow(layer) == 0) stop("No such layercode: ", layercode)
  photo <- deformed_core_photos[deformed_core_photos$photo == layer$photo,]
  coords <- deformed_layer_data[deformed_layer_data$layercode == layercode,]
  files <- file.path("data-raw", "photos", paste0(layer$photo, c(".png", ".jpg")))
  if(file.exists(files[1])) {
    img <- png::readPNG(files[1])
  } else if(file.exists(files[2])) {
    img <- jpeg::readJPEG(files[2])
  } else {
    stop("No photo for layercode ", layercode)
  }
  # do the plotting
  dims <- dim(img)
  imgBoundsX <- c(0, dims[2]/photo$scale)
  imgBoundsY <- c(0, dims[1]/photo$scale)
  bbox <- rbind(range(coords$x), range(coords$y))
  bbox <- prettymapr::zoombbox(bbox, zoom)
  coordsBoundsX <- bbox[1,]
  coordsBoundsY <- bbox[2,]
  
  plot(coordsBoundsX, imgBoundsY[2]-coordsBoundsY, pch=".", asp=1,
       axes = FALSE, frame.plot=FALSE, ann=FALSE)
  rasterImage(img, imgBoundsX[1], imgBoundsY[1], imgBoundsX[2], imgBoundsY[2])
  lines(coords$x, imgBoundsY[2]-coords$y, col=col, lwd=lwd, lty=lty)
  prettymapr::addscalebar(plotunit="cm", ...)
  if(!is.null(title)) {
    graphics::title(title)
  }
}

