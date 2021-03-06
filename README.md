simrec
======

<!-- badges: start -->
[![R build status](https://github.com/federicomarini/simrec/workflows/R-CMD-check/badge.svg)](https://github.com/federicomarini/simrec/actions)
<!-- badges: end -->

## An R-Package for Simulation of Recurrent Events

**simrec** allows simulation of recurrent event data following the multiplicative intensity model described in Andersen and Gill (1982), with the baseline hazard being a function of the total/calendar time. 

To induce between-subject-heterogeneity, a random effect covariate (frailty term) can be incorporated.

Furthermore, via **simreccomp** simulation of data following a multistate model with recurrent event data of one type and a competing event is possible.  

**simrecint** gives the possibility to additionally simulate a recruitment time for each individual and cut the data to an interim data set. 

With **simrecPlot** and **simreccompPlot** the data can be plotted.

### Installation

To install the development version for the package **simrec**, please start a current version of R and type (using `remotes`):

```r
# currently this can be done via github
install.packages("remotes") # if needed
remotes::install_github("federicomarini/simrec") # or if needed:
remotes::install_github("federicomarini/simrec", build_vignettes=TRUE, dependencies = TRUE)
```

### Vignette

To inspect the vignette and the code used in it, type:

```r
vignette("simrec-vignette-Rhelp")
## and/or
browseVignettes("simrec")
```

### Contact

For additional details regarding the functions of **simrec**, please consult the documentation.
