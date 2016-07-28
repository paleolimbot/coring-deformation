
#' Create a deformation model
#' 
#' Create a model based on a function \code{d0(d, r)}.
#'
#' @param slicesize The z range to be used.
#' @param d0 A vectorized function taking two parameters (d and r)
#' @param width The diameter of the core barrel to be modelled
#' @param cellsize The cellsize at which to model the results.
#'
#' @return A \code{data.frame} with columns x, y, z, and d0 
#' @export
#'
#' @examples
#' library(ggplot2)
#' library(dplyr)
#' model <- deformation_model(1)
#' # overhead view
#' ggplot(model %>% filter(z==0), aes(x, y)) + geom_raster(aes(fill=d0)) + 
#'   scale_fill_gradient2() + scale_y_reverse() + 
#'   coord_fixed()
#' 
#' # side on view
#' ggplot(model %>% filter(y==0), aes(x, z)) + geom_raster(aes(fill=d0)) + 
#'   scale_fill_gradient2() + scale_y_reverse() + 
#'   coord_fixed()
#' 
deformation_model <- function(slicesize, d0=function(d, r) d - 0.1 * r^2,
                              width=6, cellsize=0.1) {
  coords <- expand.grid(x=seq(-width/2, width/2, cellsize),
                        y=seq(-width/2, width/2, cellsize),
                        z=seq(-slicesize/2, slicesize/2, cellsize))
  coords$r <- sqrt(coords$x^2 + coords$y^2)
  coords <- coords[coords$r <= width/2,]
  coords$d0 <- floor(d0(coords$z, coords$r) / cellsize) * cellsize
  coords$r <- NULL
  return(coords)
}


#' Model extrusion of a deformation model
#' 
#' Simulate extruded data from "real" data.
#'
#' @param model A model as created by \link{deformation_model}.
#' @param original_data A \code{data.frame} with columns \code{d0}, \code{vals}, and \code{density}.
#'   The \code{original_data} must contain rows with a \code{d0} for every \code{cellsize/2} in the
#'   \code{model}.
#'
#' @return A \code{data.frame} with columns \code{d} and \code{vals}
#' @export
#'
#' @examples
#' library(ggplot2)
#' model <- deformation_model(1)
#' rdata <- create_random_data()
#' extruded <- extrude(model, rdata)
#' ggplot(extruded, aes(y=d, x=vals)) + geom_path() + scale_y_reverse()
#' 
extrude <- function(model, original_data) {
  odenv <- new.env()
  sapply(1:nrow(original_data), function(i) {
    odenv[[as.character(original_data$d0[i])]] <- i
  })
  paramsfunc <- function(d0) {
    return(original_data[unlist(lapply(as.character(d0[!is.na(d0)]), function(x) {
      odenv[[x]]
    }), use.names = F),])
  }
  slice_size <- max(model$z) - min(model$z)
  out <- data.frame(d=seq(min(original_data$d0)+slice_size/2, max(original_data$d0), slice_size))
  out$vals <- sapply(out$d, function(deltad0) {
    params <- paramsfunc(model$d0 + deltad0)
    params <- params[!is.na(params$d0),] # NAs are an issue for the first/last slices
    return(sum(params$vals * params$density) / sum(params$density))
  })
  out
}


#' Create random paleolimnological data
#' 
#' Creates random data for the purposes of evaluating the effect of extrusion on data.
#'
#' @param seed Seed for the random generator
#' @param smoothing Number of times to gaussion smooth the data
#' @param func The random function to use
#' @param transform The transformation to be applied to the random data
#' @param density A vector containing the two values density should range between
#' @param maxdepth The maximum depth to calculate to
#' @param cellsize The cellsize (data will be generated with resolution cellsize)
#'
#' @return A \code{data.frame} with columns \code{d0}, \code{vals}, and \code{density}.
#' @export
#'
#' @examples
#' rdata <- create_random_data()
#' ggplot(rdata, aes(y=d0, x=vals)) + geom_path() + scale_y_reverse()
#' 
create_random_data <- function(seed=2500, smoothing=10, func=rlnorm,
                               transform=function(x) 1.5 ^ x * 10,
                               density=c(0.1, 0.5),
                               maxdepth=20, cellsize=0.1) {
  set.seed(seed)
  original_data <- data.frame(d0=seq(0, maxdepth, cellsize), vals=func(maxdepth/cellsize+1))
  original_data$vals <- transform(original_data$vals)
  for(s in 1:smoothing) {
    original_data$vals <- ksmooth(original_data$d0, original_data$vals, bandwidth=cellsize*2)$y
  }
  # define density for each (linearly increasing with range)
  original_data$density <- approx(x=c(0, maxdepth), y=density, xout=original_data$d0)$y
  return(original_data)
}
