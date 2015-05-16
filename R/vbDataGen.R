#' Randomly create a von Bertalanffy data set.
#'
#' Randomly creates a data set according to a typical von Bertalanffy growth function (VBGF).
#'
#' This function can be used to generate random growth data, in one of four formats (see \sQuote{return} section and examples below) controlled by \code{dataType}.  The essential steps for constructing the simulated data are as follows:
#' \enumerate{
#'   \item Use \code{Linf}, \code{K}, and \code{t0} to set \sQuote{true} values for the three parameters in a typical VBGF (see \code{\link[FSA]{vbModels}} in \pkg{FSA} for this model).
#'   \item For each fish, generate a random age-at-capture (completed growing seasons) between \code{minage} and \code{maxage} from the probability distribution for each age given in \code{probdist}.
#'   \item Add a fractional age to the previous step and store in \code{ageFrac}.  By defaut, the \code{sample.per} is the first day of the \code{growth.per} so \code{ageCap} and \code{ageFrac} will be the same.  However, if \code{sample.per} is not a single day equal to the first day of the \code{growth.per}, than the \code{ageFrac} will be greater than \code{ageCap}.  See details for more details.
#'   \item For each fish, generate a random value for Linf, K, and t0 from a normal distribution with a mean from \code{Linf}, \code{K}, and \code{t0} and a standard deviation derived from the value in \code{Linf}, \code{K}, and \code{t0} times the value in \code{paramCV}.  This models different parameter values for individual fish.
#'   \item For each fish, generate the mean length at each \code{ageFrac} in the fish's life (i.e., if the fish is five years old, then a random length for ages one through five) according to the typical VBGF and the random parameters determined above for that fish.  These means are deterministic within a fish, but not among fish (because of the previous step).
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
#' @param growth.per A vector of length 2 that contains the Julian dates for the beginning and ending dates for the period when the fish grows.  See details and examples.
#' @param sample.per A vector of length 1 or 2 that contains the Julian dates for the beginning and ending dates for sampling the fish.  If only one value is given then it is assumed that sampling occurred only on that date.  See details and examples.
#' @param dataType The type of data frame that should be returned.  See details.
#' @param lendigs A single numeric to control the number of digits for the length variable.
#' @param agedigs A single numeric to control the number of digits for the fractional age variable.
#' @param backcalc A logical that indicates whether the lengths at previous ages should be considered to be from back-calculation (\code{TRUE}) or repeated observed (e.g., of tagged fish; \code{FALSE}, DEFAULT).
#' @param seed An integer that can be used to control the seed for the random number generator used in the function.  See details.
#' 
#' @return A data frame in one of four formats as determined by \code{dataType}.  The variables that may occur in each data.frame are:
#' \itemize{
#'   \item \code{id}: a unique identification for individual fish.
#'   \item \code{ageCap}: the number of completed growing seasons at the time of capture.
#'   \item \code{ageFrac}: the age that includes the fractional age of the fish -- i.e., number of completed growing seasons plus a decimal representation of the fraction of the current growing seasons completed.
#'   \item \code{agePrev}: an age in the fish's growth history.  Will be a whole number if \code{backcalc=TRUE} (i.e., length was back-calculated to an age that represents the beginning of a growing season) and a fractional number if \code{backcalc=FALSE} (i.e., repeated observations of length may not be at the beginning of a growing season).
#'   \item \code{lenCap}: the observed length at the time of capture (corresponds to \code{ageCap}).
#'   \item \code{len}: the observed length that corresponds to \code{agePrev}.
#' }
#' 
#' Which of these variables is returned depends on the \code{dataTYpe} as follows:
#' \itemize{
#'   \item \sQuote{atCapture} format: Returns \code{id}, \code{ageCap}, \code{ageFrac}, and \code{lenCap}.  Each row corresponds to the observations for an individual fish at the time of capture.
#'   \item \sQuote{long} format: Returns \code{id}, \code{ageCap}, \code{agePrev}, \code{ageFrac}, \code{len}, and \code{lenCap}.  Each row corresponds to one age observation and individual fish are spread across as many rows as \code{ageCap}.  If \code{backcalc=TRUE} then there will be one \code{ageFrac} greater than \code{ageCap} (represents \dQuote{plus} or current season's growth) and all other \code{ageFrac}s will be whole numbers.
  #'   \item \sQuote{groupedData} format: Returns same as \sQuote{long} format but as a \code{\link[nlme]{groupedData}} object (from the \pkg{nlme} package).
  #'   \item \sQuote{wide} format: returns \code{id}, \code{ageCap}, \code{lenCap}, and lengths-at-age (\code{ageX}; both previous ages and age-at-capture).  In \sQuote{wide} format each row corresponds to a fish and may include several lengths (at previous ages and at-capture).
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
#' vbcap <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,lendigs=1)
#' head(vbcap)
#' plot(lenCap~ageCap,data=vbcap,pch=19,col=rgb(0,0,0,1/5))
#' 
#' # example using constant errors for means and individuals
#' vbcap <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,SE=2,SD=0.5,lendigs=1)
#' head(vbcap)
#' plot(lenCap~ageCap,data=vbcap,pch=19,col=rgb(0,0,0,1/5))
#' 
#' # note the uniform (within randomization) age distribution above
#' xtabs(~ageCap,data=vbcap)
#' 
#' # example with non-uniform age distribution
#' #  (i.e., relatively few young and old fish)
#' ad <- c(1,1,2,4,4,4,2,1,1)
#' vbcap <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,agedist=ad,SECV=0.02,SDCV=0.04,lendigs=1)
#' head(vbcap)
#' plot(lenCap~ageCap,data=vbcap,pch=19,col=rgb(0,0,0,1/5))
#' xtabs(~ageCap,data=vbcap)
#'  
#' ###########################################################
#' ## long (ungrouped format) ... repeated measures
#' vblong <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,dataType="long",lendigs=1)
#' head(vblong,n=15)
#' 
#' ###########################################################
#' ## long (grouped format) ... repeated measures
#' vbgrp <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,dataType="groupedData",lendigs=1)
#' head(vbgrp)
#' if (require(lattice)) plot(vbgrp)
#' 
#' # just errors in parameters (i.e., models unique parameters
#' # for each fish but no variability around the VBGF for an
#' # individual fish).
#' vbgrp <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SE=0,SD=0,dataType="groupedData",lendigs=1)
#' head(vbgrp)
#' if (require(lattice)) plot(vbgrp)
#'
#' ###########################################################
#' ## wide format ... repeated measures
#' vbwide <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,dataType="wide",lendigs=1)
#' head(vbwide)
#' 
#' ###########################################################
#' ## demonstrate use of seed ... can get long, wide, and
#' ## atCapture formats for same "fish"
#' # set seed
#' sd <- 1534783
#' # generate three types of data using same seed
#' vbcap <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,lendigs=1,seed=sd)
#' vbwide <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,dataType="wide",lendigs=1,seed=sd)
#' vblong <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,dataType="long",lendigs=1,seed=sd)
#' # look at a couple of fish (same IDs) from each type
#' ids <- c(2,8)
#' vbcap[vbcap$id %in% ids,]
#' vbwide[vbwide$id %in% ids,]
#' vblong[vblong$id %in% ids,]
#'
#' ###########################################################
#' ## Fractional ages
#' # age-at-capture format
#' growth.per.date <- c("1-Apr-2010","21-Oct-2010")              # Assumed fish growth period
#' ( growth.per <- strptime(growth.per.date,"%d-%b-%Y")$yday+1 ) # as Julian date
#' sample.per.date <- c("15-Apr-2010","1-Jun-2010")              # Assumed sampling period
#' ( sample.per <- strptime(sample.per.date,"%d-%b-%Y")$yday+1 ) # as Julian date
#' vbcap <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,lendigs=1,
#'                    growth.per=growth.per,sample.per=sample.per)
#' head(vbcap)
#' 
#' # long format, assuming repeated observations
#' vblong <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,lendigs=1,
#'                     growth.per=growth.per,sample.per=sample.per,dataType="long")
#' head(vblong,n=15)
#' 
#' # long format, assuming back-calculated ages
#' vblong <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,lendigs=1,
#'                     growth.per=growth.per,sample.per=sample.per,dataType="long",backcalc=TRUE)
#' head(vblong,n=15)
#' 
#' ###########################################################
#' # Different minimum age
#' vbcap <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,minAge=3,SECV=0.02,SDCV=0.04,lendigs=1)
#' head(vbcap)
#' plot(lenCap~ageCap,data=vbcap,pch=19,col=rgb(0,0,0,1/5))
#' 
#' @export
vbDataGen <- function(n,Linf,K,t0,paramCV=0.1,
                      SE=NULL,SECV=NULL,SD=NULL,SDCV=NULL,
                      minAge=1,maxAge=9,agedist=rep(1,maxAge-minAge+1),
                      growth.per=c(1,365),sample.per=1,
                      dataType=c("atCapture","long","groupedData","wide"),
                      lendigs=getOption("digits"),agedigs=2,
                      backcalc=FALSE,seed=NULL) {
  len <- NULL # to avoide "global bindings" warning in rcmd check
  ## Check dataType
  dataType <- match.arg(dataType)
  ## set the random seed if one was sent
  if (!is.null(seed)) set.seed(seed)
  ## Assign random WHOLE ages to each fish
  ages <- iCreateAges(n,minAge,maxAge,agedist) 
  ## repeat fish ID and age-at-capture for each ages
  id <- rep(1:length(ages),ages)
  ageCap <- rep(ages,ages)
  ## Create previous ages (including age-at-capture)
  agePrev <- unlist(apply(matrix(ages),MARGIN=1,function(x) seq(from=1,to=x)))
  ## Create fractional ages
  ageFrac <- round(iCreateFracAges(agePrev,growth.per,sample.per),agedigs)
  ## Deal with considering lengths as back-calculated
  # put ids and age variables in a data.frame
  d <- data.frame(id,ageCap,agePrev,ageFrac)
  if (backcalc) {
    d <- iHndlBackcalc(d)
    ages <- ages+1
  }
  ## Generate growth histories for each fish
  # get typical von Bertalanffy model from vbFuns
  LVB <- vbFuns("typical")
  # Generate random parameters for each fish
  iLinf <- rep(rnorm(n,Linf,Linf*paramCV),ages)
  iK <- rep(rnorm(n,K,K*paramCV),ages)
  it0 <- rep(rnorm(n,t0,abs(t0)*paramCV),ages)
  # Generage mean lengths-at-ages based on typical VBGF
  tmp <- LVB(d$ageFrac,iLinf,iK,it0)
  # add some error to model for the mean ....
  if (!is.null(SE)) tmp <- tmp + rnorm(length(tmp),0,SE)
    else tmp <- tmp + rnorm(length(tmp),0,SECV*tmp)
  # ... and for the individual fish
  if (!is.null(SD)) tmp <- tmp + rnorm(length(tmp),0,SD)
    else tmp <- tmp + rnorm(length(tmp),0,SDCV*tmp)
  ## Find length-at-capture, repeat for each agePrev
  if (!backcalc) lenCap <- round(rep(tmp[d$ageCap==floor(d$ageFrac)],ages),lendigs)
    else lenCap <- round(rep(tmp[d$ageCap==d$ageFrac],ages),lendigs)
  ## put together as a data.frame
  d <- data.frame(d,lenCap,len=round(tmp,lendigs))
  ## Return result in format chosen by user (in dataType)
  if (dataType=="atCapture") d <- subset(d,ageCap==agePrev,c("id","ageCap","ageFrac","lenCap"))
  if (dataType=="groupedData") d <- nlme::groupedData(len~agePrev|id,data=d,labels=list(x="Age",y="Length"))
  if (dataType=="wide") {
    d$agePrev <- paste0("age",d$agePrev)
    d <- tidyr::spread(d[,-which(names(d)=="ageFrac")],agePrev,len)
  }
  d
}


# ############################################################
# INTERNAL FUNCTIONS
# ############################################################
# ============================================================
# Create random WHOLE ages
# ============================================================
iCreateAges <- function(n,minAge,maxAge,agedist) {
  ## Assign random ages to each fish according to distribution
  ## in agedist and min and max ages
  # First some checks
  if (minAge<1) stop("'minAge' argument must be greater than or equal to 1.",call.=FALSE)
  if (minAge>maxAge) stop("'minAge' can NOT be greater than 'maxAge'.",call.=FALSE)
  if (length(agedist)!=(maxAge-minAge+1)) stop("'agedist' must be same length as 'maxAge'-'minAge'+1.",call.=FALSE)
  # Convert agedist to a distribution of proportions
  agedist <- agedist/sum(agedist)
  # Generate and return the ages
  sample(minAge:maxAge,prob=agedist,size=n,replace=TRUE)
}

# ============================================================
# Convert WHOLE ages to FRACTIONAL ages
# ============================================================
iCreateFracAges <- function(age,growth.per,sample.per) {
  # check growth.per
  if (length(growth.per)!=2) stop("'growth.per' must have beginning and ending dates.",call.=FALSE)
  if (any(growth.per<0 | growth.per>365)) stop("'growth.per' values must be between 0 and 365.",call.=FALSE)
  if (growth.per[1]>growth.per[2]) growth.per <- growth.per[c(2,1)]
  # check sample.per
  if (length(sample.per)==1) sample.per <- rep(sample.per,2)
  if (length(sample.per)!=2) stop("'sample.per' must be one date or have beginning and ending dates.",call.=FALSE)
  if (any(sample.per<0 | sample.per>365)) stop("'sample.per' values must be between 0 and 365.",call.=FALSE)
  if (sample.per[1]>sample.per[2]) sample.per <- sample.per[c(2,1)]
  # Choose a sample data from within the sample period for each age
  smpl.date <- sample(sample.per[1]:sample.per[2],length(age),replace=TRUE)
  # Convert the sample date to a fraction of the growth season.
  # If sampled before growth starts then fraction is 0, if sampled
  # after growth ceases then fraction is 1, if sampled within the
  # growth season then convert to a fraction completed of the the
  # growth seasons
  frac <- (smpl.date-growth.per[1])/(growth.per[2]-growth.per[1])
  frac[frac < 0] <- 0
  frac[frac > 1] <- 1
  # add fractions to ages and return
  age+frac
}

iHndlBackcalc <- function(d,ages) {
  # isolate all at-capture info
  tmp <- d[d$ageCap==d$agePrev,]
  tmp$agePrev <- NA
  # reduce ageFrac to agePrev
  d$ageFrac <- d$agePrev
  # row bind at-capture info
  d <- rbind(d,tmp)
  d <- d[order(d$id,d$ageFrac),]
  rownames(d) <- 1:nrow(d)
  d
}
