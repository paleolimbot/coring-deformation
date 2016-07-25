``` r
library(ggplot2)
library(reshape2)
library(dplyr)
```

Introduction
------------

That deformation and compression of lake sediment occurs during coring has long been known (Martin and Miller 1982; Wright 1993), and designs of new coring devices have strived to minimize the conditions that promote deformation during coring (Martin and Miller 1982; Lane and Tafts 2002). Compression of sediment occurs during coring is a widely accepted phenomenon (Glew et al. 2001), however convex upwards deformation, while widely observed (citations including Wright 1993; Rosenbaum et al. 2010; Figure 1), is infrequently discussed. The idea that horizontal sectioning (extrusion) of deformed sediment is undesirable has been previously noted (Rosenbaum et al. 2010), however the degree to which this deformation occurs and the effect that deformation has on paleolimnological data derived from horizontal sectioning has never been investigated quantitatively.

Rather than suggest that deformation does not occur (Glew et al. 2001 Figure 3), or that a particular coring method prevents this from happening (Smol 2009 p35), we suggest that acknowledging deformation and its effect on paleolimnological data is a more reasonable approach. We suspect, given the innumerable paleolimnological studies that use coring and extrusion to produce reasonable and reproducable results, that either deformation or its effect on the data is minimal. This paper is our attempt to quantify and constrain the degree to which convex upwards deformation adds bias to horizontally sectioned paleolimnological data.

If we take a slice through a deformed core (as we would while extruding), what is the distribution of original depths? What is the effect of typical deformation on data obtained from extruded samples? What is the effect of maximum likely deformation in extruded samples? Is there a minimum extrusion interval based on typical deformation in sediment cores?

![Figure 1: pictures of deformed cores](deformed%20cores/def_cores_summaryfig.png)

Methods
-------

``` r
# define parameters
barrel_width <- 6 # use 8 cm for barrel width
slice_sizes <- c(0.1, 0.5, 1)
cellsize <- 0.1
```

### Core photo analysis

To calculate parameters for the deformation model, we loaded XXX scale photos of deformed cores into ImageJ software and digitized deformed strata. Coordinates were transformed to `ri` and `di` values for individual strata by subtracting the minimum `d` value from the rest of the values, and subtracting the average `x` value from the rest of the values. Power regression (quadratic) was performed on the data to obtain the coefficients for minimum, maximum, and mean levels of deformation. Corer type and core barrel diameter were recorded with these data.

### Deformation model

We modeled horizontal sections with height `H` and diameter `D` as a 3-dimensional raster grid with a cell size of 0.25 mm (Figure 2). For each cell `i`, an original depth `oi` was calculated with the minimum, maximum, and mean parameters obtained from digitized strata. Density histograms were then obtained to estimate the percentage contribution of each original depth `o` to the slice. For each slice, `d=0` refers to the middle of the slice. We produced these models for `D`=8 cm, as this represents a frequently used barrel diameter for horizontally sectioned cores.

![model vars](figs/fig2_deformation_vars.png)

### Effect on paleolimnological data

To model the concentration we would obtain by sectioning and homogenizing a sample with variable concentration (`c`) and density (`sigma`), we need to calculate total mass of the target substance divided by the mass of the slice. In sigma notation with our deformation model (using `n` cells), this looks like the following:

![](figs/math1.jpg)

Because the volume for each cell is constant, we can remove the `Vi` from the summation in both the numerator and denominator. Thus we are left with:

![](figs/math2.jpg)

Our initial assumption is that `c` and `sigma` are a function of (and constatnt for each) `d0i`, our notation collapses to:

![](figs/math3.jpg)

Thus we can use our deformation model to simulate the effect of core barrel width, horizontal sectioning interval, and degree of deformation, on idealized high-resolution data. We used fictional paleolimnological data (random log-normal data) to visualize the effect of deformation on stratigraphic data. These data were inspired by ITRAX core scanner data that has a high sample resolution.

``` r
# generate data at resolution cellsize/2 (200 samples for 20 cm). smooth random normal data
set.seed(2500)
original_data <- data.frame(d0=seq(0, 20, cellsize/2), vals=rlnorm(20/cellsize*2+1))
original_data$vals <- 1.5 ^ original_data$vals * 10
smoothing = 10
for(s in 1:smoothing) {
  original_data$vals <- ksmooth(original_data$d0, original_data$vals, bandwidth=0.1)$y
}
# define density for each (produces curve from 0.1 to 0.5 ish)
original_data$density <- (original_data$d0 * (5-1) / 20 + 1) * 0.1

# ggplot(original_data, aes(y=d0, x=vals)) + geom_path() + scale_y_reverse()

original_data_params <- function(d0, tol=1e-13) {
    return(original_data[sapply(d0, function(x) {
      which(abs(original_data$d0-x) < tol)[1]
    }),])
}
```

Results
-------

### Core photo analysis

![pictures of deformed cores](deformed%20cores/def_cores_summaryfig.png)

``` r
# load
deformations <- read.csv("deformed cores/deformed_cores_summary.csv")
# calc the x0 and y0 to use as the base coordinates (to transform x and y to r and d)
layers <- deformations %>% 
  group_by(layercode) %>% 
  summarise(x0=mean(range(x)), y0=min(y))
deformations <- merge(deformations, layers, by="layercode")
deformations$r <- deformations$x-deformations$x0
deformations$d <- deformations$y-deformations$y0
deformations <- deformations %>% select(layercode, core, layer, r, d)

# plot cores with quadratic smoothing
ggplot(deformations, aes(x=r, y=d)) + 
  geom_point(aes(col=factor(layer))) + 
  stat_smooth(method=lm, formula=y ~ poly(x, 2, raw=TRUE), se=FALSE) + 
  scale_y_reverse() + facet_wrap(~core, scales="free")
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
# define the regression function
create_quadratic_model <- function(df) {
  model <- lm(data=df, formula=d ~ poly(r, 2, raw=TRUE))
  return(data.frame(a=model$coefficients[3],
                    r2=summary(model)$r.squared,
                    df=model$df.residual))
}

# calculate quadratic models
modelparams <- deformations %>% 
  group_by(layercode) %>%
  do(create_quadratic_model(.))

# these are used in the model later
mina <- round(min(modelparams$a), 3)
maxa <- round(max(modelparams$a), 2)
meana <- round(mean(modelparams$a), 2)
coeffs <- c(0, mina, meana, maxa)
```

We digitized 26 deformed layers from 7 scale photos of split cores. The quadratic regression performed produced an excellent fit of the data (r2 from 0.83 to 1). Coefficients for `x^2` range from 0.083 to 0.51, with a mean of 0.25.

``` r
knitr::kable(modelparams, digits=2)
```

| layercode          |     a|    r2|   df|
|:-------------------|-----:|-----:|----:|
| longlake\_pc1/1    |  0.51|  0.87|   17|
| menounos\_cheak1/1 |  0.28|  0.94|    6|
| menounos\_cheak1/2 |  0.20|  0.95|    9|
| menounos\_cheak1/3 |  0.16|  0.93|    9|
| menounos\_cheak1/4 |  0.13|  0.85|    8|
| menounos\_cheak1/5 |  0.15|  0.87|    8|
| menounos\_cheak1/6 |  0.17|  0.92|    7|
| menounos\_cheak1/7 |  0.18|  0.99|    5|
| menounos\_cheak1/8 |  0.14|  0.92|    9|
| menounos\_cheak2/1 |  0.30|  0.99|    6|
| menounos\_cheak2/2 |  0.40|  1.00|    7|
| menounos\_cheak2/3 |  0.37|  1.00|    8|
| menounos\_cheak2/4 |  0.29|  0.96|    6|
| menounos\_cheak2/5 |  0.38|  0.99|    7|
| menounos\_cheak2/6 |  0.38|  0.99|    7|
| menounos\_cheak2/7 |  0.24|  0.96|    6|
| menounos\_cheak2/8 |  0.26|  0.99|    6|
| sl2/1              |  0.13|  0.98|    8|
| sl2/2              |  0.12|  0.99|   11|
| sl2/3              |  0.08|  0.96|    9|
| suzie1/1           |  0.16|  0.93|    7|
| suzie1/2           |  0.27|  0.83|   11|
| suzie1/3           |  0.39|  1.00|   12|
| suzie1/4           |  0.45|  0.93|   15|
| whistler\_gc4/1    |  0.16|  0.99|    8|
| whistler\_gc8/1    |  0.10|  1.00|    9|

### Deformation model

``` r
# This function takes a `slicesize` (for example, 1 cm) and a function `d0` that
# returns a the original depth prior to deformation given a depth (`d`) and a
# radius from the centre of the core (`r`). The `cellsize` parameter specifies
# the resolution at which to create the model (low values for higher resolution)
# and the `width` specifies the interior diameter of the core barrel. The `d`
# values in the model are relative to the centre of the slice. For a single
# function `d0` and slice size, we can then plot a cross-section and histogram
# of the depths prior to deformation.

# output is a data.frame with columsn x, y, z, r, and d0 (which we should refactor to oi...)

deformation_model3d <- function(slicesize, d0, width=2, cellsize=0.025, snap=NULL) {
  if(is.null(snap)) {
    snap <- cellsize
  }
  coords <- expand.grid(x=seq(-width/2, width/2, cellsize),
                        y=seq(-width/2, width/2, cellsize),
                        z=seq(-slicesize/2, slicesize/2, cellsize))
  coords$r <- sqrt(coords$x^2 + coords$y^2)
  coords <- coords[coords$r <= width/2,]
  coords$d0 <- floor(d0(coords$z, coords$r) / cellsize) * cellsize
  return(coords)
}

# expand to grid
params <- expand.grid(slicesize=slice_sizes, coeff=coeffs)

# define a wrapper function that will create the model based on the
# parameters we would like to vary
# include stat of concentration of our value
deformation_model_wrapper <- function(slicesize, coeff) {
  deformation_model3d(slicesize, 
                      d0=function(d, r) d - coeff * r^2, 
                      width = barrel_width, 
                      cellsize = cellsize)
}

# create the models in long form so we can plot with facets
all <- params %>% 
  group_by(slicesize, coeff) %>% 
  do(deformation_model_wrapper(.$slicesize, .$coeff))

# overhead view (same for all sizes of slice)
ggplot(all %>% filter(abs(z)==min(abs(z)), slicesize==all$slicesize[1]), aes(x=x, y=y)) + 
  geom_raster(aes(fill=d0)) + 
  scale_fill_gradient2() + scale_y_reverse() + 
  coord_fixed() + facet_wrap(~coeff)
```

![](README_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
# side on view
ggplot(all %>% filter(abs(y)==min(abs(y))), aes(x=x, y=z)) + geom_raster(aes(fill=d0)) + 
  scale_fill_gradient2() + scale_y_reverse() + 
  coord_fixed() + facet_grid(slicesize ~ coeff)
```

![](README_files/figure-markdown_github/unnamed-chunk-6-2.png)

``` r
# create smoothing models (histograms with binwidth of cellsize, except with the $density param)
create_smoothing_model <- function(vals) {
  limit <- round(max(abs(vals)) / cellsize) * cellsize #snap to cellsize
  breaks <- seq(-limit-cellsize, limit+cellsize, cellsize)
  h <- hist(vals, breaks=breaks, plot=FALSE)
  # fill smoothing model with zeros to make symmetrical
  return(data.frame(d=breaks[-1], density=h$density))
}

models <- all %>% 
  group_by(slicesize, coeff) %>% 
  do(create_smoothing_model(.$d0))

# plot (same as density histogram but this is faster)
ggplot(models %>% filter(density != 0), aes(x=d, y=density)) + 
  geom_bar(width=cellsize, stat="identity") + 
  facet_grid(slicesize ~ coeff)
```

![](README_files/figure-markdown_github/unnamed-chunk-6-3.png)

Smoothing effects on stratigraphic data
---------------------------------------

The distribution (density) acts as a smoothing filter on geochem data.

``` r
# we need a function that simulates figs/math3.jpg
extrude <- function(model, paramsfunc, slice_size) {
    out <- data.frame(d=seq(0+slice_size/2, max(original_data$d0), slice_size))
    out$vals <- sapply(out$d, function(deltad0) {
      params <- paramsfunc(model$d0 + deltad0)
      params <- params[!is.na(params$d0),] # NAs are an issue for the first/last slices
      return(sum(params$vals * params$density) / sum(params$density))
    })
    out
}

sliced <- all %>%
  group_by(slicesize, coeff) %>% 
  do(extrude(., original_data_params, unique(.$slicesize)))

ggplot(sliced, aes(y=d, x=vals)) + geom_path() + scale_y_reverse() + 
  facet_grid(slicesize ~ coeff, scales="free_x")
```

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)

Conclusions
-----------

There is a limit to how small extrusion intervals can get based on deformation (flat topped distributions are bad because they don't have a 'mode' depth that they represent!).
