#' Randomly create a von Bertalanffy data set.
#'
#' Randomly creates a data set according to a typical von Bertalanffy growth function (VBGF).
#'
#' This function can be used to generate random growth data, in one of four formats (see \sQuote{return} section and examples below) controlled by \code{dataType}.  The essential steps for constructing the simulated data are as follows:
#' \enumerate{
#'   \item Use \code{Linf}, \code{K}, and \code{t0} to set \sQuote{true} values for the three parameters in a typical VBGF (see \code{\link[FSA]{vbModels}} in \pkg{FSA} for this model).
#'   \item For each fish, generate a random age-at-capture between \code{minage} and \code{maxage} from the probability distribution for each age given in \code{probdist}.
#'   \item For each fish, generate a random value for Linf, K, and t0 from a normal distribution with a mean from \code{Linf}, \code{K}, and \code{t0} and a standard deviation derived from the value in \code{Linf}, \code{K}, and \code{t0} times the value in \code{paramCV}.  This models different parameter values for individual fish.
#'   \item For each fish, generate the mean length at each age in the fish's life (i.e., if the fish is five years old, then a random length for ages one through five) according to the typical VBGF and the random parameters determined above for that fish.  These means are deterministic within a fish, but not among fish (because of the previous step).
#'   \item For each fish, add a random error to the mean lengths-at-age derived in the previous step.  These random errors are from a normal distribution with a mean of 0 and a standard error given by \code{SE} or \code{SECV} times the mean length-at-age.  This models random error in the mean lengths-at-age within an individual fish.
#'   \item For each fish, add a random error to the mean lengths-at-age derived in the previous step.  These random errors are from a normal distribution with a mean of 0 and a standard deviation given by \code{SD} or \code{SDCV} times the mean length-at-age.  This models random error of individual lengths around the mean lengths-at-age within an individual fish.
#'   \item Repeat the previous steps (except the first) for all fish to be simulated.
#'   \item Store the data for all fish in a data.frame, modify that data.frame to match the format required by \code{dataType}, and return the data.frame.
#' }
#' 
#' Note that \code{seed} can be used to control the random seed used by the function.  Use the same \code{seed} with different \code{dataType}s to return the same data in different formats.  See examples.
#'
#' @param n Number of individuals to simulate data for.
#' @param Linf The \sQuote{true} value of the Linf parameter.
#' @param K The \sQuote{true} value of the K parameter.
#' @param t0 The \sQuote{true} value of the t0 parameter.
#' @param paramCV The coefficient of variation around the parameters for individual fish.  See details.
#' @param SE The constant standard error around the mean lengths-at-age within a fish.  See details.
#' @param SECV The coefficient of variation around the mean lengths-at-age within a fish.  See details.
#' @param SD The constant standard deviation of individuals around the mean lengths-at-age within a fish.  See details.
#' @param SDCV The coefficient of variation for individuals around the mean lengths-at-age within a fish.  See details.
#' @param minAge The minimum possible age-at-capture for an individual in the data set.
#' @param maxAge The maximum possible age-at-capture for an individual in the data set.
#' @param agedist A vector of values that are proportional to the percentage of fish expected for each age.  This vector must have a length equal to \code{maxAge}-\code{minAge}+1.  See examples.
#' @param dataType The type of data frame that should be returned.  See details.
#' @param digits A single numeric to control the number of digits for the length variable.
#' @param seed An integer that can be used to control the seed for the random number generator used in the function.  See details.
#' 
#' @return A data frame in one of four formats as determined by \code{dataType}:
#' \enumerate{
#'   \item The \sQuote{atCapture} format returns the fish's ID (\code{id}) and the the age (\code{age}) and length (\code{len}) at the time of capture.
#'   \item The \sQuote{wide} format returns the fish's ID (\code{id}), age-at-capture (\code{agecap}), and lengths at age (\code{anuX}; both previous ages and age-at-capture).  In \sQuote{wide} format each row corresponds to a fish and may include several lengths (at previous ages and at-capture).
#'   \item In \sQuote{long} format the same data as in \sQuote{wide} format is returned but reshaped to long format where each row contains the fish's ID (\code{id}), age-at-capture (\code{agecap}), age (\code{prvAge}) and length (\code{anu}) in the fish's history.  In \sQuote{long} format the data for an individual fish will be found in as many rows as the age of the fish.  Thus, the data for a five-year-old fish will appear in five rows corresponding to the ages and lengths at each of the fish's five years of life.  The length-at-capture is in the row where \code{agecap} equals \code{anu}.
#'   \item The \sQuote{groupedData} format is the same as \sQuote{long} format but as a \code{\link[nlme]{groupedData}} object (from the \pkg{nlme} package).
#' }
#' 
#' @author Derek H. Ogle, \email{dogle@@northland.edu}.  This function was motivated by the simulation methodology shown in Box 2 of Vigliola and Meekan (2009).
#'
#' @seealso See \code{\link[FSA]{gReshape}} for how growth data is reshaped from \sQuote{wide} to \sQuote{long} format.  See \code{\link[nlme]{groupedData}} from \pkg{nlme} for how to analyze data in \sQuote{groupedData} format.
#'
#' @references Vigliola, L. and M.G. Meekan.  2009.  The back-calculation of fish growth from otoliths.  Chapter 6 (pp. 174-211) in B.S. Green et al. (eds.), Tropical Fish Otoliths: Information for Assessment, Management and Ecology, Reviews: Methods and Technologies in Fish Biology and Fisheries 11.
#'
#' @keywords manip
#' 
#' @examples
#' ## At-capture (observation) format
#' # example using CV as errors for means and individuals
#' vbcap <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,digits=1)
#' head(vbcap)
#' plot(len~age,data=vbcap,pch=19,col=rgb(0,0,0,1/5))
#' 
#' # example using constant errors for means and individuals
#' vbcap <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,SE=2,SD=0.5,digits=1)
#' head(vbcap)
#' plot(len~age,data=vbcap,pch=19,col=rgb(0,0,0,1/5))
#' # note the uniform (within randomization) age distribution
#' xtabs(~age,data=vbcap)
#' 
#' # example with non-uniform age distribution
#' #  (i.e., relatively few young and old fish)
#' ad <- c(1,1,2,4,4,4,2,1,1)
#' vbcap <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,agedist=ad,SECV=0.02,SDCV=0.04,digits=1)
#' head(vbcap)
#' plot(len~age,data=vbcap,pch=19,col=rgb(0,0,0,1/5))
#' xtabs(~age,data=vbcap)
#'  
#' ###########################################################
#' ## wide format ... repeated measures
#' vbwide <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,dataType="wide",digits=1)
#' head(vbwide)
#' 
#' ###########################################################
#' ## long (ungrouped format) ... repeated measures
#' vblong <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,dataType="long",digits=1)
#' head(vblong,n=15)
#' 
#' ###########################################################
#' ## long (grouped format) ... repeated measures
#' vbgrp <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,dataType="groupedData",digits=1)
#' head(vbgrp)
#' if (require(lattice)) plot(vbgrp)
#' 
#' # just errors in parameters (i.e., models unique parameters
#' # for each fish but no variability around the VBGF for an
#' # individual fish).
#' vbgrp <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SE=0,SD=0,dataType="groupedData",digits=1)
#' head(vbgrp)
#' if (require(lattice)) plot(vbgrp)
#' 
#' ###########################################################
#' ## demonstrate use of seed ... can get long, wide, and
#' ## atCapture formats for same "fish"
#' # set seed
#' sd <- 1534783
#' # generate three types of data using same seed
#' vbcap <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,digits=1,seed=sd)
#' vbwide <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,dataType="wide",digits=1,seed=sd)
#' vblong <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,dataType="long",digits=1,seed=sd)
#' # look at a couple of fish (same IDs) from each type
#' ids <- c(2,8)
#' vbcap[vbcap$id %in% ids,]
#' vbwide[vbwide$id %in% ids,]
#' vblong[vblong$id %in% ids,]
#' 
#' @export
vbDataGen <- function(n,Linf,K,t0,paramCV=0.1,
                      SE=NULL,SECV=NULL,SD=NULL,SDCV=NULL,
                      minAge=1,maxAge=9,agedist=rep(1,maxAge-minAge+1),
                      dataType=c("atCapture","wide","long","groupedData"),
                      digits=getOption("digits"),seed=NULL) {
  ## Check dataType
  dataType <- match.arg(dataType)
  ## get typical von Bertalanffy model from vbFuns
  LVB <- vbFuns("typical")
  ## set the random seed if one was sent
  if (!is.null(seed)) set.seed(seed)
  ## assign random ages to each fish according to distribution
  ## in agedist and min and max ages
  #    first some checks
  if (minAge<1) stop("'minAge' argument must be greater than or equal to 1.",call.=FALSE)
  if (minAge>maxAge) stop("'minAge' can NOT be greater than 'maxAge'.",call.=FALSE)
  if (length(agedist)!=(maxAge-minAge+1)) stop("'agedist' must be same length as 'maxAge'-'minAge'+1.",call.=FALSE)
  #    convert agedist to a distribution of proportions
  agedist <- agedist/sum(agedist)
  #    generate the ages
  ages <- sample(minAge:maxAge,prob=agedist,size=n,replace=TRUE) 
  ## Generate random parameters for each fish
   iLinf <- rnorm(n,Linf,Linf*paramCV)
   iK <- rnorm(n,K,K*paramCV)
   it0 <- rnorm(n,t0,abs(t0)*paramCV)
  ## Generate size-at-age data (note that number of ages is maximum age here)
  numAges <- maxAge
  #    suppress warnings about NAs
  options(warn= -1)
  #    fill out a matrix of past measurements
  L <- matrix(NA,nrow=n,ncol=numAges)
  for (j in 1:n) {
    i <- 1:ages[j]
    # add NAs for ages that won't exist for this fish
    if (length(i)<numAges) i <- c(i,rep(NA,numAges-length(i)))
    # generage mean lengths based on typical VBGF
    L[j,] <- LVB(i,iLinf[j],iK[j],it0[j])
    # add some error to model for the mean ....
    if (!is.null(SE)) L[j,] <- L[j,] + rnorm(numAges,0,SE)
      else L[j,] <- L[j,] + rnorm(numAges,0,SECV*L[j,])
    # ... and for the individual fish
    if (!is.null(SD)) L[j,] <- L[j,] + rnorm(numAges,0,SD)
      else L[j,] <- L[j,] + rnorm(numAges,0,SDCV*L[j,])
  }
  #    allow warnings to return
  options(warn=0)
  #    put together as a data.frame
  d <- data.frame(1:n,ages,round(L,digits))
  names(d) <- c("id","agecap",paste0("anu",1:maxAge))
  ## Return result in format chosen by user (in dataType)
  if (dataType %in% c("long","groupedData")) {
    # reshape
    d <- gReshape(d,in.pre="anu",val.name="anu")
    # re-order
    d <- d[order(d$id,d$prvAge),]
    # change row-numbers
    rownames(d) <- 1:nrow(d)
  }
  if (dataType=="groupedData") d <- groupedData(anu~prvAge|id,data=d,labels=list(x="Age",y="Length"))
  if (dataType=="atCapture") {
    # finds max value in anuX to be length-at-capture and appends to ID and agecap
    d <- data.frame(id=d$id,age=d$age,len=apply(d[,3:ncol(d)],1,max,na.rm=TRUE))
  }
  d
}
