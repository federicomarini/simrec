#' simrecint
#'
#' With this function previously simulated data (for example simulated by the use of \code{simrec} or \code{simreccomp})
#' can be cut to an interim data set.
#' The simulated data must be in patient time (i.e. time since the patient entered the study),
#' and must be in the counting process format. Furthermore the dataset must have the variables \code{id}, \code{start}, \code{stop} and \code{status},
#' like data simulated by the use of \code{simrec} or \code{simreccomp}.
#' Then for every individual additionally a recruitment time is generated in study time (i.e. time since start of the study),
#' which is uniformly distributed on \code{[0, tR]}.
#' The timing of the interim analysis \code{tI} is set in study time and
#' data are being cut to all data, that are available at the interim analysis.
#' For further explanations on study time and patient time see the vignette.
#' If you only wish to simulate a recruitment time, \code{tI} can be set to \code{tR + fu.max} or something bigger.
#'
#' @param data  Previously generated data (in patient time), that shall be cut to interim data
#' @param N     Number of individuals, for which \code{data} was generated
#' @param tR    Length of the recruitment period (in study time)
#' @param tI    Timing of the interim analysis (in study time)
#'
#' @return The output is a data.frame consisting of the columns, that were put into, and additionally the following columns:
#'    \item{rectime}{The recruitment time for each individual (in study time).}
#'    \item{interimtime}{The time of the interim analysis \code{tI} (in study time).}
#'    \item{stop_study}{The stopping time for each event in study time.}
#'  Individuals that are not already recruited at the interim analysis are left out here.
#' @author Katharina Ingel, Stella Preussler, Antje Jahn-Eimermacher.
#' Institute of Medical Biostatistics, Epidemiology and Informatics (IMBEI),
#' University Medical Center of the Johannes Gutenberg-University Mainz, Germany
#' @seealso simrec, simreccomp
#' @export
#' @examples
#' ### Example - see example for simrec
#' library(simrec)
#' N         <- 10
#' dist.x    <- c("binomial", "normal")
#' par.x     <- list(0.5, c(0,1))
#' beta.x    <- c(0.3, 0.2)
#' dist.z    <- "gamma"
#' par.z     <- 0.25
#' dist.rec  <- "weibull"
#' par.rec   <- c(1,2)
#' fu.min    <- 2
#' fu.max    <- 2
#' cens.prob <- 0.2
#'
#' simdata <- simrec(N, fu.min, fu.max, cens.prob, dist.x, par.x, beta.x, dist.z,
#'                   par.z, dist.rec, par.rec)
#'
#' ### Now simulate for each patient a recruitment time in [0,tR=2]
#' ### and cut data to the time of the interim analysis at tI=1:
#'
#' simdataint <- simrecint(simdata, N=N, tR=2, tI=1)
#' # print(simdataint)  # only run for small N!
simrecint <- function(data, 
                      N, 
                      tR, 
                      tI){
  
  tab     <- data                # previously generated data with columns including id, start, stop, status
  nr.rows <- length(tab$id)
  
  ev <- NULL           # number of "events" for each individual (including censoring)
  for (j in 1:N) {
    ev[j] <- length(which(tab$id==j))
  }
  
  rectime         <- runif(N, min = 0, max = tR)                       # generate recruitment time for each individual
  tab$rectime     <- rep(rectime, ev)
  tab$interimtime <- rep(tI, nr.rows)
  tab$stop_study  <- tab$stop + tab$rectime                            # stopping time in study time
  
  m <- (tab$stop_study > tab$interimtime)                           # if m[i]=TRUE -> event after interim analysis
  for (i in 1:nr.rows) {                                               # Then time of interim analysis minus time of recruitment gives new stop time (in patient time)
    if (m[i]) { tab$stop[i]   <- tab$interimtime[i] - tab$rectime[i]   # this is negative, if patient is not recruited yet (see also below)
    tab$status[i] <- 0                                     # subsequent lines now equal, if more than one event after interim analysis
    }                                                        # Events after interim analysis are now censored at time of interim (status = 0)
  }
  
  k <- NULL
  for (i in 1:nr.rows-1) {
    k[i] <- (tab$stop[i]==tab$stop[i+1] & tab$id[i]==tab$id[i+1])      # if stop-time and ID in the following line are the same, set to TRUE
  }
  
  for (i in 1:length(k)) {
    if (k[i]) {is.na(tab[i+1,]) <- TRUE }     # cuts all events after interim time (accordinlgy to k as generated above), by setting the whole line to NA
  }
  
  tab <- na.omit(tab)
  tab <- tab[tab$rectime < tab$interimtime,]      # only patients, that were already recruited at interim
  
  return(tab)
}
