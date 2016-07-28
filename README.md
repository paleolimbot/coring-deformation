``` r
library(ggplot2)
library(dplyr)
library(multidplyr) # devtools::install_github("hadley/multidplyr")
source("R/deformation_library.R") # functions
source("data-raw/collect_layers.R") # layers data
source("R/layer_plot.R") # layer plotting
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
barrel_width <- 6 # use 6 cm for barrel width
slice_sizes <- c(0.1, 0.5, 1)
cellsize <- 0.05
```

### Core photo analysis

To calculate parameters for the deformation model, we loaded 13 scale photos of deformed cores from 3 sources into ImageJ software and digitized deformed strata (Table 1). Coordinates were transformed to `ri` and `di` values for individual strata by subtracting the minimum `d` value from the rest of the values, and subtracting the central `x` value from the rest of the values. Power regression (quadratic) was performed on the data to obtain the coefficients for minimum, maximum, and mean levels of deformation. Corer type and core barrel diameter were recorded with these data.

``` r
knitr::kable(deformed_core_photos, digits=2)
```

| photo            |   scale| reference                                 |
|:-----------------|-------:|:------------------------------------------|
| longlake\_pc1    |  102.60| White 2012                                |
| menounos\_cheak1 |   88.00| Menounos and Clague 2008                  |
| menounos\_cheak2 |   46.00| Menounos and Clague 2008                  |
| suzielake\_1     |   47.03| Spooner et al. 1997                       |
| suzielake\_2     |   26.70| Spooner et al. 1997                       |
| whistler\_gc4    |  290.25| Dunnington 2015                           |
| whistler\_gc8    |  310.23| Dunnington 2015                           |
| crevice\_lake    |   14.77| Rosenbaum et al. 2010                     |
| menounos\_cheak3 |   49.00| Menounos et al. 2005                      |
| ds\_unpubl1      |   87.74| Dunnington and Spooner (unpublished data) |
| ds\_unpubl2      |   68.36| Dunnington and Spooner (unpublished data) |
| ds\_unpubl3      |   82.00| Dunnington and Spooner (unpublished data) |
| ds\_unpubl4      |   61.87| Dunnington and Spooner (unpublished data) |

### Deformation model

We modeled horizontal sections with height `H` and diameter `D` as a 3-dimensional raster grid with a cell size of 0.25 mm (Figure 2). For each cell `i`, an original depth `d0i` was calculated with the minimum, maximum, and mean parameters obtained from digitized strata. Density histograms were then obtained to estimate the percentage contribution of each original depth `d0` to the slice. For each slice, `d=0` refers to the middle of the slice. We produced these models for `D`=6 cm, as this represents a frequently used barrel diameter for horizontally sectioned cores.

![model vars](figs/fig2_deformation_vars.png)

### Effect on paleolimnological data

To model the concentration we would obtain by sectioning and homogenizing a sample with variable concentration (`c`) and density (`density`), we need to calculate total mass of the target substance divided by the mass of the slice. In sigma notation with our deformation model (using `n` cells), this looks like the following:

![](figs/math1.jpg)

Because the volume for each cell is constant, we can remove the `Vi` from the summation in both the numerator and denominator. Thus we are left with:

![](figs/math2.jpg)

Our initial assumption is that `c` and `density` are a function of (and constatnt for each) `d0i`, our notation collapses to:

![](figs/math3.jpg)

Thus we can use our deformation model to simulate the effect of core barrel width, horizontal sectioning interval, and degree of deformation, on idealized high-resolution data. We used fictional paleolimnological data (random log-normal data) to visualize the effect of deformation on stratigraphic data. These data were inspired by ITRAX core scanner data that has a high sample resolution.

``` r
# generate data
original_data <- create_random_data(seed=2500, smoothing=10, func=rlnorm,
                                    transform=function(x) 1.5 ^ x * 10,
                                    density=c(0.1, 0.5),
                                    maxdepth=20, cellsize=cellsize)

# ggplot(original_data, aes(y=d0, x=vals)) + geom_path() + scale_y_reverse()
```

Results
-------

### Core photo analysis

``` r
ggplot(deformed_layers, aes(a)) + geom_histogram(bins=30)
```

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png)<!-- -->

``` r
# plot core photos
op <- par(mfrow=c(2, 3), oma=c(0,0,0,0), mar=c(0.1,0.1,1,0.1))
plot_layer("crevice_lake.5", title="crevice_lake.5 (a=0.10)") # 0.1
plot_layer("crevice_lake.12", title="crevice_lake.12 (a=0.20)") # 0.2
plot_layer("menounos_cheak2.2", title="menounos_cheak2.2 (a=0.40)") # 0.4

plot_layer("whistler_gc8.1", title="whistler_gc8.1 (a=0.10)") # 0.1
plot_layer("menounos_cheak1.7", title="menounos_cheak1.7 (a=0.19)") # 0.2
plot_layer("suzielake_2.9", title="suzielake_2.9 (a=0.41)") # 0.4
```

![](README_files/figure-markdown_github/unnamed-chunk-5-2.png)<!-- -->

``` r
par(op); rm(op)

# these are used in the model later
mina <- round(min(deformed_layers$a), 3)
maxa <- round(max(deformed_layers$a), 2)
meana <- round(mean(deformed_layers$a), 2)
coeffs <- c(0, 0.1, 0.2, 0.4)
```

We digitized 0 deformed layers from 0 scale photos of split cores. The quadratic regression performed produced an excellent fit of the data (r2 from 0.58 to 1). Coefficients for `x^2` range from 0.054 to 0.51, with a mean of 0.21.

``` r
knitr::kable(deformed_layers, digits=2)
```

| photo            | layercode          |     a|    r2|   df|
|:-----------------|:-------------------|-----:|-----:|----:|
| crevice\_lake    | crevice\_lake.1    |  0.06|  1.00|    6|
| crevice\_lake    | crevice\_lake.10   |  0.14|  0.98|    9|
| crevice\_lake    | crevice\_lake.11   |  0.17|  0.89|    8|
| crevice\_lake    | crevice\_lake.12   |  0.20|  0.91|    6|
| crevice\_lake    | crevice\_lake.2    |  0.05|  0.99|    7|
| crevice\_lake    | crevice\_lake.3    |  0.07|  0.98|    9|
| crevice\_lake    | crevice\_lake.4    |  0.08|  0.97|    7|
| crevice\_lake    | crevice\_lake.5    |  0.10|  0.98|    7|
| crevice\_lake    | crevice\_lake.6    |  0.21|  0.92|    9|
| crevice\_lake    | crevice\_lake.7    |  0.19|  1.00|    6|
| crevice\_lake    | crevice\_lake.8    |  0.25|  0.98|    8|
| crevice\_lake    | crevice\_lake.9    |  0.23|  0.98|    7|
| ds\_unpubl1      | ds\_unpubl1.1      |  0.08|  1.00|    4|
| ds\_unpubl2      | ds\_unpubl2.1      |  0.16|  0.74|    7|
| ds\_unpubl2      | ds\_unpubl2.2      |  0.14|  0.58|    7|
| ds\_unpubl3      | ds\_unpubl3.1      |  0.16|  0.71|    8|
| ds\_unpubl4      | ds\_unpubl4.1      |  0.21|  0.74|    9|
| longlake\_pc1    | longlake\_pc1.1    |  0.51|  0.87|   17|
| menounos\_cheak1 | menounos\_cheak1.1 |  0.27|  0.94|    6|
| menounos\_cheak1 | menounos\_cheak1.2 |  0.20|  0.94|    9|
| menounos\_cheak1 | menounos\_cheak1.3 |  0.16|  0.93|    9|
| menounos\_cheak1 | menounos\_cheak1.4 |  0.12|  0.85|    8|
| menounos\_cheak1 | menounos\_cheak1.5 |  0.15|  0.87|    8|
| menounos\_cheak1 | menounos\_cheak1.6 |  0.17|  0.93|    7|
| menounos\_cheak1 | menounos\_cheak1.7 |  0.19|  0.99|    5|
| menounos\_cheak1 | menounos\_cheak1.8 |  0.14|  0.92|    9|
| menounos\_cheak2 | menounos\_cheak2.1 |  0.30|  0.99|    6|
| menounos\_cheak2 | menounos\_cheak2.2 |  0.40|  1.00|    7|
| menounos\_cheak2 | menounos\_cheak2.3 |  0.37|  1.00|    8|
| menounos\_cheak2 | menounos\_cheak2.4 |  0.29|  0.96|    6|
| menounos\_cheak2 | menounos\_cheak2.5 |  0.38|  0.99|    7|
| menounos\_cheak2 | menounos\_cheak2.6 |  0.38|  0.99|    7|
| menounos\_cheak2 | menounos\_cheak2.7 |  0.24|  0.96|    6|
| menounos\_cheak2 | menounos\_cheak2.8 |  0.26|  0.99|    6|
| suzielake\_1     | suzielake\_1.1     |  0.16|  0.93|    7|
| suzielake\_1     | suzielake\_1.2     |  0.26|  0.83|   11|
| suzielake\_1     | suzielake\_1.3     |  0.39|  1.00|   12|
| suzielake\_1     | suzielake\_1.4     |  0.45|  0.94|   15|
| suzielake\_2     | suzielake\_2.1     |  0.16|  0.94|    7|
| suzielake\_2     | suzielake\_2.2     |  0.13|  0.98|    7|
| suzielake\_2     | suzielake\_2.3     |  0.07|  0.98|    6|
| suzielake\_2     | suzielake\_2.4     |  0.14|  0.95|    8|
| suzielake\_2     | suzielake\_2.5     |  0.29|  0.97|    8|
| suzielake\_2     | suzielake\_2.6     |  0.24|  0.85|    8|
| suzielake\_2     | suzielake\_2.7     |  0.26|  0.99|    7|
| suzielake\_2     | suzielake\_2.8     |  0.25|  0.91|    8|
| suzielake\_2     | suzielake\_2.9     |  0.41|  0.99|    9|
| whistler\_gc4    | whistler\_gc4.1    |  0.16|  0.99|    8|
| whistler\_gc8    | whistler\_gc8.1    |  0.10|  1.00|    9|

### Deformation model

``` r
# expand to grid
params <- expand.grid(slicesize=slice_sizes, coeff=coeffs)

# define a wrapper function that will create the model based on the
# parameters we would like to vary
deformation_model_wrapper <- function(slicesize, coeff) {
  deformation_model(slicesize, 
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

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)<!-- -->

``` r
# side on view
ggplot(all %>% filter(abs(y)==min(abs(y))), aes(x=x, y=z)) + geom_raster(aes(fill=d0)) + 
  scale_fill_gradient2() + scale_y_reverse() + 
  coord_fixed() + facet_grid(slicesize ~ coeff)
```

![](README_files/figure-markdown_github/unnamed-chunk-7-2.png)<!-- -->

``` r
# create histograms with binwidth of cellsize, except with the $density param)
create_histogram <- function(vals) {
  limit <- round(max(abs(vals)) / cellsize) * cellsize #snap to cellsize
  breaks <- seq(-limit-cellsize, limit+cellsize, cellsize)
  h <- hist(vals, breaks=breaks, plot=FALSE)
  # fill smoothing model with zeros to make symmetrical
  return(data.frame(d=breaks[-1], density=h$density))
}

histograms <- all %>% 
  group_by(slicesize, coeff) %>% 
  do(create_histogram(.$d0))

# plot (same as density histogram but this is faster)
ggplot(histograms %>% filter(density != 0), aes(x=d, y=density)) + 
  geom_bar(width=cellsize, stat="identity") + 
  facet_grid(slicesize ~ coeff)
```

![](README_files/figure-markdown_github/unnamed-chunk-7-3.png)<!-- -->

### Effect on stratigraphic data

Using the formulas, if we model extrusion, this is the effect on the data:

``` r
# use multicore to speed this up
cluster <- create_cluster(3)
cluster_assign_value(cluster, "extrude", extrude)
cluster_assign_value(cluster, "original_data", original_data)
sliced <- all %>%
  group_by(slicesize, coeff) %>% 
  partition(slicesize, coeff, cluster=cluster) %>%
  do(extrude(., original_data)) %>%
  collect()

ggplot(sliced, aes(y=d, x=vals)) + geom_path() + scale_y_reverse() + 
  facet_grid(slicesize ~ coeff, scales="free_x")
```

![](README_files/figure-markdown_github/unnamed-chunk-8-1.png)<!-- -->

Conclusions
-----------

There is a limit to how small extrusion intervals can get based on deformation. For minor deformation, even small extrusion intervals are ok.
