#' simrec
#'
#' This function allows simulation of recurrent event data following the multiplicative
#' intensity model described in Andersen and Gill [1] with the baseline hazard being a
#' function of the total/calendar time. To induce between-subject-heterogeneity a random
#' effect covariate (frailty term) can be incorporated. Data for individual \eqn{i} are generated
#' according to the intensity process \deqn{Y_i(t) * \lambda_0(t)* Z_i *exp(\beta^t X_i),}
#' where \eqn{X_i} defines the covariate vector and \eqn{\beta} the regression coefficient vector.
#' \eqn{\lambda_0(t)} denotes the baseline hazard, being a function of the total/calendar
#' time \eqn{t}, and \eqn{Y_i(t)} the predictable process
#' that equals one as long as individual \eqn{i} is under observation and at risk for experiencing events.
#' \eqn{Z_i} denotes the frailty variable with \eqn{(Z_i)_i} iid with \eqn{E(Z_i)=1} and
#' \eqn{Var(Z_i)=\theta}. The parameter \eqn{\theta} describes the degree of
#' between-subject-heterogeneity. Data output is in the counting process format.
#'
#' @param N          Number of individuals
#' @param fu.min     Minimum length of follow-up.
#' @param fu.max     Maximum length of follow-up. Individuals length of follow-up is
#'  generated from a uniform distribution on
#'   \code{[fu.min, fu.max]}. If \code{fu.min=fu.max}, then all individuals have a common
#'   follow-up.
#' @param cens.prob  Gives the probability of being censored due to loss to follow-up before
#'   \code{fu.max}. For a random set of individuals defined by a B(N,\code{cens.prob})-distribution,
#'   the time to censoring is generated from a uniform
#'   distribution on \code{[0, fu.max]}. Default is \code{cens.prob=0}, i.e. no censoring
#'   due to loss to follow-up.
#' @param dist.x     Distribution of the covariate(s) \eqn{X}. If there is more than one covariate,
#'   \code{dist.x} must be a vector of distributions with one entry for each covariate. Possible
#'   values are \code{"binomial"} and \code{"normal"}, default is \code{dist.x="binomial"}.
#' @param par.x      Parameters of the covariate distribution(s). For \code{"binomial", par.x} is
#'   the probability for \eqn{x=1}. For \code{"normal"}, \code{par.x=c(}\eqn{\mu, \sigma}\code{)}
#'   where \eqn{\mu} is the mean and \eqn{\sigma} is the standard deviation of a normal distribution.
#'   If one of the covariates is defined to be normally distributed, \code{par.x} must be a list,
#'   e.g. \code{ dist.x <- c("binomial", "normal")} and \code{par.x  <- list(0.5, c(1,2))}.
#'   Default is \code{par.x=0}, i.e. \eqn{x=0} for all individuals.
#' @param beta.x  Regression coefficient(s) for the covariate(s) \eqn{x}. If there is more than one
#'   covariate, \code{beta.x} must be a vector of coefficients with one entry for each covariate.
#'   \code{simrec} generates as many covariates as there are entries in \code{beta.x}. Default is
#'   \code{beta.x=0}, corresponding to no effect of the covariate \eqn{x}.
#' @param dist.z     Distribution of the frailty variable \eqn{Z} with \eqn{E(Z)=1} and
#'   \eqn{Var(Z)=\theta}. Possible values are \code{"gamma"} for a Gamma distributed frailty
#'    and \code{"lognormal"} for a lognormal distributed frailty.
#'    Default is \code{dist.z="gamma"}.
#' @param par.z      Parameter \eqn{\theta} for the frailty distribution: this parameter gives
#'   the variance of the frailty variable \eqn{Z}.
#'   Default is \code{par.z=0}, which causes \eqn{Z=1}, i.e. no frailty effect.
#' @param dist.rec   Form of the baseline hazard function. Possible values are \code{"weibull"} or
#'   \code{"gompertz"} or \code{"lognormal"} or \code{"step"}.
#' @param par.rec  Parameters for the distribution of the event data.
#'   If \code{dist.rec="weibull"} the  hazard function is \deqn{\lambda_0(t)=\lambda*\nu* t^{\nu - 1},}
#'   where \eqn{\lambda>0} is the scale and \eqn{\nu>0} is the shape parameter. Then
#'   \code{par.rec=c(}\eqn{\lambda, \nu}\code{)}. A special case
#'   of this is the exponential distribution for \eqn{\nu=1}.\\
#'   If \code{dist.rec="gompertz"}, the hazard function is \deqn{\lambda_0(t)=\lambda*exp(\alpha t),}
#'   where \eqn{\lambda>0} is the scale and \eqn{\alpha\in(-\infty,+\infty)} is the shape parameter.
#'   Then \code{par.rec=c(}\eqn{\lambda, \alpha}\code{)}.\\
#'   If \code{dist.rec="lognormal"}, the hazard function is
#'   \deqn{\lambda_0(t)=[(1/(\sigma t))*\phi((ln(t)-\mu)/\sigma)]/[\Phi((-ln(t)-\mu)/\sigma)],}
#'   where \eqn{\phi} is the probability density function and \eqn{\Phi} is the cumulative
#'   distribution function of the standard normal distribution, \eqn{\mu\in(-\infty,+\infty)} is a
#'   location parameter and \eqn{\sigma>0} is a shape parameter. Then \code{par.rec=c(}\eqn{\mu,\sigma}\code{)}.
#'   Please note, that specifying \code{dist.rec="lognormal"} together with some covariates does not
#'   specify the usual lognormal model (with covariates specified as effects on the parameters of the
#'   lognormal distribution resulting in non-proportional hazards), but only defines the baseline
#'   hazard and incorporates covariate effects using the proportional hazard assumtion.\\
#'   If \code{dist.rec="step"} the hazard function is \deqn{\lambda_0(t)=a, t<=t_1, and \lambda_0(t)=b, t>t_1}.
#'   Then \code{par.rec=c(}\eqn{a,b,t_1}\code{)}.
#' @param pfree Probability that after experiencing an event the individual is not at risk
#'   for experiencing further events for a length of \code{dfree} time units.
#'   Default is \code{pfree=0}.
#' @param dfree Length of the risk-free interval. Must be in the same time unit as \code{fu.max}.
#'  Default is \code{dfree=0}, i.e. the individual is continously at risk for experiencing
#'  events until end of follow-up.
#' @return The output is a data.frame consisting of the columns:
#'    \item{id}{An integer number for identification of each individual}
#'    \item{x}{or \code{x.V1, x.V2, ...} - depending on the covariate matrix. Contains the
#'       randomly generated value of the covariate(s) \eqn{X} for each individual.}
#'    \item{z}{Contains the randomly generated value of the frailty variable \eqn{Z} for each individual.}
#'    \item{start}{The start of interval \code{[start, stop]}, when the individual
#'       starts to be at risk for a next event.}
#'    \item{stop}{The time of an event or censoring, i.e. the end of interval
#'    \code{[start, stop]}.}
#'    \item{status}{An indicator of whether an event occured at time \code{stop} (\code{status=1})
#'        or the individual is censored at time \code{stop} (\code{status=0}).}
#'    \item{fu}{Length of follow-up period \code{[0,fu]} for each individual.}
#'    For each individual there are as many lines as it experiences events,
#'    plus one line if being censored.
#'    The data format corresponds to the counting process format.
#' @details
#' Data are simulated by extending the methods proposed by Bender et al [2]
#' to the multiplicative intensity model.
#' @references
#' \enumerate{
#' \item Andersen P, Gill R (1982): Cox's regression model for counting processes:
#'    a large sample study. The Annals of Statistics 10:1100-1120
#' \item Bender R, Augustin T, Blettner M (2005): Generating survival times to simulate Cox
#'   proportional hazards models. Statistics in Medicine 24:1713-1723
#' \item Jahn-Eimermacher A, Ingel K, Ozga AK, Preussler S, Binder H (2015): Simulating recurrent event data
#'   with hazard functions defined on a total time scale. BMC Medical Research Methodology 15:16
#' }
#' @author Katharina Ingel, Stella Preussler, Antje Jahn-Eimermacher.
#' Institute of Medical Biostatistics, Epidemiology and Informatics (IMBEI),
#' University Medical Center of the Johannes Gutenberg-University Mainz, Germany
#' @seealso simreccomp
#' @export
#' @examples
#' ### Example:
#' ### A sample of 10 individuals
#'
#' N <- 10
#'
#' ### with a binomially distributed covariate with a regression coefficient
#' ### of beta=0.3, and a standard normally distributed covariate with a
#' ### regression coefficient of beta=0.2,
#'
#' dist.x <- c("binomial", "normal")
#' par.x  <- list(0.5, c(0,1))
#' beta.x   <- c(0.3, 0.2)
#'
#' ### a gamma distributed frailty variable with variance 0.25
#'
#' dist.z <- "gamma"
#' par.z  <- 0.25
#'
#' ### and a Weibull-shaped baseline hazard with shape parameter lambda=1
#' ### and scale parameter nu=2.
#'
#' dist.rec <- "weibull"
#' par.rec  <- c(1,2)
#'
#' ### Subjects are to be followed for two years with 20% of the subjects
#' ### being censored according to a uniformly distributed censoring time
#' ### within [0,2] (in years).
#'
#' fu.min    <- 2
#' fu.max    <- 2
#' cens.prob <- 0.2
#'
#' ### After each event a subject is not at risk for experiencing further events
#' ### for a period of 30 days with a probability of 50%.
#'
#' dfree <- 30/365
#' pfree <- 0.5
#'
#' simdata <- simrec(N, fu.min, fu.max, cens.prob, dist.x, par.x, beta.x, dist.z, par.z,
#'                   dist.rec, par.rec, pfree, dfree)
#' # print(simdata)  # only run for small N!
simrec<- function(N, 
                  fu.min, 
                  fu.max, 
                  cens.prob=0, 
                  dist.x="binomial", 
                  par.x=0, 
                  beta.x=0, 
                  dist.z="gamma", 
                  par.z=0, 
                  dist.rec, 
                  par.rec, 
                  pfree=0, 
                  dfree=0) {
  
  ID <- c(1:N)
  
  # generating the follow-up
  # follow-up uniformly distributed in [fu.min, fu.max] if not censored
  # or uniformly distributed in [0, fu.max] if censored
  if(cens.prob<0 || cens.prob>1){stop("cens.prob must be a probability between 0 and 1")}
  if(fu.min>fu.max || fu.min<0){stop("fu.min must be a non-negative value smaller or equal fu.max")}
  fu <- rbinom(N, 1, cens.prob)           # 1 = censored
  nr.cens <- sum(fu)
  if(nr.cens==0){                         # nobody censored
    fu <- runif(N, min=fu.min, max=fu.max)
  } else{
    index.cens <- which(fu==1)
    fu[-index.cens] <- runif((N-nr.cens), min=fu.min, max=fu.max)
    fu[index.cens]  <- runif(nr.cens, min=0, max=fu.max)
  }
  
  if(length(beta.x)!=length(dist.x)){stop("dimensions of beta.x and dist.x differ")}
  if(length(beta.x)!=length(par.x)){stop("dimensions of beta.x and par.x differ")}
  
  # generating the covariate-matrix x
  nr.cov <- length(beta.x)       # number of covariates
  x<-matrix(0,N,nr.cov)        # matrix with N lines and one column for each covariate
  for(i in 1:nr.cov){
    dist.x[i] <- match.arg(dist.x[i], choices=c("binomial", "normal"))
    if (dist.x[i]=="binomial") {
      if(length(par.x[[i]])!=1){stop("par.x has wrong dimension")}
      if(par.x[[i]]<0 || par.x[[i]]>1){stop("par.x must be a probability between 0 and 1 for the binomial distributed covariate")}
      x[,i] <- c(rbinom(N,1,par.x[[i]]))
    } else {                                         # normally distributed covariate
      if(length(par.x[[i]])!=2){stop("par.x has wrong dimension")}
      mu.x    <- par.x[[i]][1]
      sigma.x <- par.x[[i]][2]
      x[,i]   <- c(rnorm(N, mean=mu.x, sd=sigma.x))
    }
  }
  
  # generating the frailty variable z
  z <- rep(1,N)
  dist.z <- match.arg(dist.z, choices=c("gamma", "lognormal"))
  if(length(par.z)!=1){stop("par.z has wrong dimension")}
  if(par.z < 0) {stop("par.z must be non-negative")}
  if(par.z!=0){                                       # if par.z=0 then frailty=1 for all
    if(dist.z=="gamma"){                              # gamma-frailty
      aGamma <- 1/par.z
      z      <- rgamma(N, shape=aGamma, scale=1/aGamma)
    } else {                                          # lognormal frailty
      mu.z    <- log(1/sqrt(par.z + 1))
      sigma.z <- sqrt(log(par.z+1))
      z       <- exp(rnorm(N, mean = mu.z, sd=sigma.z))
    }
  }
  
  # derivation of the distributional parameters for the recurrent event data
  dist.rec <- match.arg(dist.rec, choices=c("weibull", "lognormal", "gompertz", "step"))
  if (dist.rec=="lognormal") {                                       # lognormal
    if(length(par.rec)!=2){stop("par.rec has wrong dimension")}
    mu    <- par.rec[1]
    sigma <- par.rec[2]
    if(any(beta.x!=0)){
      warning("lognormal together with covariates specified: this does not define the usual lognormal model! see help for details")
    }
  } else if(dist.rec=="weibull"){                                    # weibull
    if(length(par.rec)!=2){stop("par.rec has wrong dimension")}
    lambda <- par.rec[1]
    nu     <- par.rec[2]
  } else if(dist.rec=="gompertz"){
    if(length(par.rec)!=2){stop("par.rec has wrong dimension")}      # gompertz
    lambdag <- par.rec[1]
    alpha   <- par.rec[2]
  } else if(dist.rec=="step"){                                       # step
    if(length(par.rec)!=3){stop("par.rec has wrong dimensions")}
    fc   <- par.rec[1]
    sc   <- par.rec[2]
    jump       <- par.rec[3]
    jumpinv    <- jump*fc
  }
  
  if(length(pfree)!=1){stop("pfree has wrong dimension")}
  if(length(dfree)!=1){stop("dfree has wrong dimension")}
  if(pfree<0 || pfree>1){stop("pfree must be a probability between 0 and 1")}
  
  # initial step: simulation of N first event times
  U <- runif(N)
  Y <- (-1)*log(U)*exp((-1)*x%*%beta.x)*1/z
  if (dist.rec=="lognormal") {                   # lognormal
    t <- exp(qnorm(1-exp((-1)*Y))*sigma + mu)
  } else if (dist.rec=="weibull"){               # weibull
    t <- ((lambda)^(-1)*Y)^(1/nu)
  } else if(dist.rec=="gompertz"){               # gompertz
    t <- (1/alpha)*log((alpha/lambdag)*Y+1)
  } else if(dist.rec=="step"){                   # step
    t <- rep(NA,N)
    indexTr1    <- which(Y <= jumpinv)
    if(length(indexTr1 >0)){ t[indexTr1] <- Y[indexTr1]/fc }
    indexTr2    <- which(Y>jumpinv)
    if(length(indexTr2 >0)){ t[indexTr2] <- (Y[indexTr2]-(fc-sc)*jump)/sc }
  }
  T  <- matrix(t,N,1)
  dirty <- rep(TRUE,N)
  T1 <- NULL
  
  # recursive step: simulation of N subsequent event times
  while (any(dirty)) {
    pd <- rbinom(N,1,pfree)
    U  <- runif(N)
    Y  <- (-1)*log(U)*exp((-1)*x%*%beta.x)*1/z
    t1 <- t+pd*dfree
    if (dist.rec=="lognormal") {                                                     # lognormal
      t <- (t1 + exp(qnorm(1-exp(log(1-pnorm((log(t1)-mu/sigma)))-Y))*sigma+mu)-(t1))
    } else if (dist.rec=="weibull"){                                                 # weibull
      t <- (t1 + ((Y+lambda*(t1)^(nu))/lambda)^(1/nu)-(t1))
    } else if(dist.rec=="gompertz"){                                                 # gompertz
      t <- (t1 + ((1/alpha)*log((alpha/lambdag)*Y+exp(alpha*t1))) - (t1))
    } else if(dist.rec=="step"){                                                      # step
      indexTr3 <- which((t1 <=jump) & (Y <= (jump-t1)*fc))
      if(length(indexTr3 >0)){ t[indexTr3] <- t1[indexTr3] + Y[indexTr3]/fc }
      indexTr4 <- which((t1 <= jump) & (Y > (jump-t1)*fc))
      if(length(indexTr4 >0)){ t[indexTr4] <- t1[indexTr4] + (Y[indexTr4]+(fc-sc)*(t1[indexTr4]-jump))/sc }
      indexTr5 <- which(t1 > jump)
      if(length(indexTr5 >0)){ t[indexTr5] <- t1[indexTr5] + Y[indexTr5]/sc }
    }
    T1 <- cbind(T1,ifelse(dirty,t1,NA))
    dirty <- ifelse(dirty,(t(t) < fu) & (t(t1) < fu),dirty)
    if(!any(dirty)) break
    T <- cbind(T,ifelse(dirty,t,NA))
  }
  
  # start times
  start.t     <- cbind(0,T1)
  start.t     <- as.vector(t(start.t))
  tab.start.t <- start.t[!is.na(start.t)]
  
  # stop times
  stop.t <- cbind(T,NA)
  d <- apply(!is.na(T),1,sum)       # number of events per individual
  f <- d+1
  for (i in 1:N){
    stop.t[i,f[i]] <- fu[i]
  }
  stop.t     <- as.vector(t(stop.t))
  tab.stop.t <- stop.t[!is.na(stop.t)]
  
  # deriving the censoring indicator variable and truncating stop times that are larger than FU
  e <- NULL
  for (i in 1:N) {
    e <- cbind(e,t(rep(1,d[i])),0)
  }
  
  tab.ID <- rep(ID,f)
  tab.X  <- x[rep(1:nrow(x),f),]
  tab.Z  <- rep(z,f)
  tab.Fu <- rep(fu,f)
  
  w <- tab.start.t>tab.stop.t
  v <- rep(0,length(w))
  for (i in 1:length(w)){
    if (w[i]) {v[i-1] <- 1}
  }
  
  l <- tab.stop.t>tab.Fu
  for (i in 1:length(l)) {
    if (l[i]) {tab.stop.t[i] <- tab.Fu[i]; e[i]<-0}
  }
  
  tab <- cbind(tab.ID,tab.X,tab.Z,tab.start.t,tab.stop.t,t(e),tab.Fu)
  for (i in 1:length(w)) {
    if(w[i]) {tab[i,] <- rep(NA, nr.cov+6)}
  }
  tab<- data.frame(id=tab[,1], x=tab[,2:(nr.cov+1)], z=tab[,(nr.cov+2)], start=tab[,(nr.cov+3)], stop=tab[,(nr.cov+4)], status=tab[,(nr.cov+5)], fu=tab[,(nr.cov+6)])
  tab<-na.omit(tab)
  return(tab)
}
