#' simrecPlot
#'
#' This function allows plotting of recurrent event data.
#'
#' @param data       A data set of recurrent event data to be plotted.
#'                      The input-data must include columns corresponding to:
#'                      \code{id} (patient-ID), \code{start} (= beginning of an interval where the patient is at risk for an event),
#'                      \code{stop} (= end of the interval due to an event or censoring),
#'                      \code{status} (= an indicator of the patient status at \code{stop} with = 0 censoring, 1 = event)
#' @param id            the name of the \code{id} column, default is \code{"id"}
#' @param start         the name of the \code{start} column, default is \code{"start"}
#' @param Stop          the name of the \code{stop} column, default is \code{"stop"}
#' @param status        the name of the \code{status} column, default is \code{"status"}
#'
#' @return The output  is a plot of the data with a bullet indicating a recurrent event and a circle indicating censoring.
#' @author Katharina Ingel, Stella Preussler, Antje Jahn-Eimermacher.
#' Institute of Medical Biostatistics, Epidemiology and Informatics (IMBEI),
#' University Medical Center of the Johannes Gutenberg-University Mainz, Germany
#' @seealso simrec, simreccomp, simreccompPlot
#' @export
#' @examples
#' ### Example:
#' ### First simulate a sample of 10 individuals (for more details see the help of \code{simrec})
#' N <- 10
#' dist.x <- c("binomial", "normal")
#' par.x  <- list(0.5, c(0,1))
#' beta.x   <- c(0.3, 0.2)
#' dist.z <- "gamma"
#' par.z  <- 0.25
#' dist.rec <- "weibull"
#' par.rec  <- c(1,2)
#' fu.min    <- 2
#' fu.max    <- 2
#' cens.prob <- 0.2
#' dfree <- 30/365
#' pfree <- 0.5
#' simdata <- simrec(N, fu.min, fu.max, cens.prob, dist.x, par.x, beta.x,
#'                   dist.z, par.z, dist.rec, par.rec, pfree, dfree)
#' simrecPlot(simdata)
#' # or possibly:
#' # pdf("Plot.pdf")
#' # simrecPlot(simdata)
#' # dev.off()
simrecPlot <- function(data, 
                       id="id", 
                       start="start", 
                       Stop="stop", 
                       status="status"){
  
  if(!(id %in% colnames(data))){stop("Please give the name of the id-column")}
  if(!(start %in% colnames(data))){stop("Please give the name of the start-column")}
  if(!(Stop %in% colnames(data))){stop("Please give the name of the stop-column")}
  if(!(status %in% colnames(data))){stop("Please give the name of the status-column")}
  
  colnames(data)[colnames(data)==id]     <- "id"
  colnames(data)[colnames(data)==start]  <- "start"
  colnames(data)[colnames(data)==Stop]   <- "stop"
  colnames(data)[colnames(data)==status] <- "status"
  
  data    <- data[order(data$id),]       # data ordered by id
  t       <- table(data$id)              # the table entries will also be ordered by id
  idvec   <- names(t)                    # all occuring IDs just one time
  idnum   <- seq(along=idvec)            # number the patients consecutively (corresponding to the ordering by id)
  idtable <- data.frame(idvec, idnum)
  data$idnum   <- NULL
  for(i in idvec){
    data$idnum[data$id==i] <- idtable$idnum[idtable$idvec==i]   # new column with a numerical id for each patient
  }
  
  events <- data$stop[data$status==1]       # all event times in one vector
  cens   <- data$stop[data$status==0]       # all censoring times in one vector
  
  idevents <- data$idnum[data$status==1]    # all numerical patient-IDs corresponding to the event times in one vector
  idcens   <- data$idnum[data$status==0]    # all numerical patient-IDs corresponding to the censoring times in one vector
  
  nevents <- length(events)                 # how many event times    ...  if no events: events = numeric(0) and length(events)=0
  ncens   <- length(cens)                   # how many censoring times
  
  par(las=1, cex.axis=0.5)
  plot( c(events, cens), c(idevents, idcens),                              # plot time points vs. corresponding patient-id
        pch=c(rep(20, nevents), rep(1,ncens)),                             # symbol for event: filled circle / symbol for censoring: circle
        xlab = "time", ylab="patient", yaxt="n", #axes=FALSE,
        main="event history of patients")
  for(i in idvec){
    datai <- subset(data, id==i)
    for(j in seq_along(datai$id)){
      lines(c(datai$start[j],datai$stop[j]), c(datai$idnum[j],datai$idnum[j]), lty="solid")          # add lines for each start-stop intervall
    }
  }
  # axis(1, labels=TRUE, at=0:max(data$fu), tick=TRUE)
  axis(2, labels=idvec, at=idnum, tick=TRUE)
}


#' simreccompPlot
#'
#' This function allows plotting of recurrent event data with a competing event.
#'
#' @param data          A data set of recurrent event data to be plotted.
#'                      The input-data must include columns corresponding to:
#'                      \code{id} (patient-ID), \code{start} (= beginning of an interval where the patient is at risk for an event),
#'                      \code{stop} (= end of the interval due to an event or censoring),
#'                      \code{status} (= an indicator of the patient status at \code{stop} with = 0 censoring, 1 = event, 2 = competing event)
#' @param id            the name of the \code{id} column, default is \code{"id"}
#' @param start         the name of the \code{start} column, default is \code{"start"}
#' @param Stop          the name of the \code{stop} column, default is \code{"stop"}
#' @param status        the name of the \code{status} column, default is \code{"status"}
#' 
#' @return The output  is a plot of the data with a bullet indicating a recurrent event, an x indicating the competing event and a circle indicating censoring.
#' @author Katharina Ingel, Stella Preussler, Antje Jahn-Eimermacher.
#' Institute of Medical Biostatistics, Epidemiology and Informatics (IMBEI),
#' University Medical Center of the Johannes Gutenberg-University Mainz, Germany
#' @seealso simrec, simreccomp, simrecPlot
#' @export
#' @examples
#' ### Example:
#' ### First simulate a sample of 10 individuals (for more details see the help of \code{simreccomp})
#' N <- 10
#' dist.x <- c("binomial", "normal")
#' par.x  <- list(0.5, c(0,1))
#' beta.xr <- c(0.3, 0.2)
#' beta.xc <- c(0.5,0.25)
#' dist.zr <- "gamma"
#' par.zr  <- 0.25
#' dist.zc <- "gamma"
#' par.zc  <- 0.3
#' dist.rec <- "weibull"
#' par.rec  <- c(1,2)
#' dist.comp <- "weibull"
#' par.comp  <-c(1,2)
#' fu.min    <- 2
#' fu.max    <- 2
#' cens.prob <- 0.2
#' dfree <- 30/365
#' pfree <- 0.5
#' simdata <- simreccomp(N=N, fu.min=fu.min, fu.max=fu.max, cens.prob=cens.prob,
#'                       dist.x=dist.x, par.x=par.x, beta.xr=beta.xr, beta.xc=beta.xc,
#'                       dist.zr=dist.zr, par.zr=par.zr,dist.zc=dist.zc, par.zc=par.zc,
#'                       dist.rec=dist.rec, par.rec=par.rec, dist.comp=dist.comp, par.comp=par.comp,
#'                       pfree= pfree, dfree=dfree)
#' simreccompPlot(simdata)
#' # or possibly:
#' # pdf("Plot.pdf")
#' # simreccompPlot(simdata)
#' # dev.off()
simreccompPlot <- function(data, 
                           id="id", 
                           start="start", 
                           Stop="stop", 
                           status="status"){
  
  if(!(id %in% colnames(data))){stop("Please give the name of the id-column")}
  if(!(start %in% colnames(data))){stop("Please give the name of the start-column")}
  if(!(Stop %in% colnames(data))){stop("Please give the name of the stop-column")}
  if(!(status %in% colnames(data))){stop("Please give the name of the status-column")}
  
  colnames(data)[colnames(data)==id]     <- "id"
  colnames(data)[colnames(data)==start]  <- "start"
  colnames(data)[colnames(data)==Stop]   <- "stop"
  colnames(data)[colnames(data)==status] <- "status"
  
  data    <- data[order(data$id),]       # data ordered by id
  t       <- table(data$id)              # the table entries will also be ordered by id
  idvec   <- names(t)                    # all occuring IDs just one time
  idnum   <- seq(along=idvec)            # number the patients consecutively (corresponding to the ordering by id)
  idtable <- data.frame(idvec, idnum)
  data$idnum   <- NULL
  for(i in idvec){
    data$idnum[data$id==i] <- idtable$idnum[idtable$idvec==i]   # new column with a numerical id for each patient
  }
  
  events <- data$stop[data$status==1]    # all event times in one vector
  cens   <- data$stop[data$status==0]    # all censoring times in one vector
  compev <- data$stop[data$status==2]    # all competing event times in one vector
  
  idevents <- data$idnum[data$status==1]    # all numerical patient-IDs corresponding to the event times in one vector
  idcens   <- data$idnum[data$status==0]    # all numerical patient-IDs corresponding to the censoring times in one vector
  idcompev <- data$idnum[data$status==2]    # all numerical patient-IDs corresponding to the competing event times in one vector
  
  nevents <- length(events)              # how many event times                if no events: events = numeric(0) and length(events)=0
  ncens   <- length(cens)                # how many censoring times
  ncompev <- length(compev)              # how many competing events
  
  par(las=1, cex.axis=0.5)
  plot( c(events, cens, compev), c(idevents, idcens, idcompev),         # plot time points vs. corresponding patient-id
        pch=c(rep(20, nevents), rep(1,ncens), rep(4, ncompev)),         # symbol for event: filled circle / symbol for censoring: circle / symbol for comp. event: x
        xlab = "time", ylab="patient", yaxt="n", #axes=FALSE,
        main="event history of patients")
  for(i in idvec){
    datai <- subset(data, id==i)
    for(j in seq_along(datai$id)){
      lines(c(datai$start[j],datai$stop[j]), c(datai$idnum[j],datai$idnum[j]), lty="solid")        # add lines for each start-stop intervall
    }
  }
  # axis(1, labels=TRUE, at=0:max(data$fu), tick=TRUE)
  axis(2, labels=idvec, at=idnum, tick=TRUE)
}
