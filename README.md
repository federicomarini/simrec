simrec
======

#### An R-Package for Simulation of Recurrent Event
**simrec** allows simulation of recurrent event data following the multiplicative intensity model described in
Andersen and Gill (1982) with the baseline hazard being a function of the total/calendar time. To induce 
between-subject-heterogeneity a random effect covariate (frailty term) can be incorporated.

### Installation
To install the development version for the package **simrec**, please start a current version of R and type (using `devtools`):

```r 
# currently this can be done via github
install.packages("devtools") # if needed
devtools::install_github("katharinaingel/simrec") # or if needed:
devtools::install_github("katharinaingel/simrec", build_vignettes=TRUE)
```

### Vignette
To inspect the vignette and the code used in it, type:

```r
vignette("simrec-vignette-Rhelp")
## and/or
browseVignettes("simrec")
```

### Contact
For additional details regarding the functions of **simrec**, please consult the documentation or write an email to ingel@uni-mainz.de. 
