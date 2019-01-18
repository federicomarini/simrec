#' simrec
#' 
#' Simulation of recurrent event data for non-constant baseline hazard
#' (total-time model)
#' 
#' Simulation of recurrent event data for non-constant baseline hazard 
#' in the total time model with risk-free intervalls and possibly a competing 
#' event. The simrec package enables to cut the data to an interim data set, 
#' and provides functionality to plot.
#' 
#' @importFrom graphics axis lines par plot
#' @importFrom stats na.omit pnorm qnorm rbinom rgamma rnorm runif
#' 
#' @author
#' Katharina Ingel, Stella Preussler, Antje Jahn-Eimermacher, Federico Marini
#' 
#' Maintainer: Antje Jahn-Eimermacher \email{jahna@uni-mainz.de}
#' 
#' @name simrec-pkg
#' @docType package
NULL