#' @title Randomly create length-at-age data according to a von Bertalanffy growth function.
#'
#' @description Randomly creates length-at-age data according to a typical von Bertalanffy growth function (VBGF).
#'
#' @param n Number of individuals to simulate data for.
#' @param Linf The \sQuote{true} value of the Linf parameter.
#' @param K The \sQuote{true} value of the K parameter.
#' @param t0 The \sQuote{true} value of the t0 parameter.
#' @param paramCV The coefficient of variation for the parameters for individual fish.  See details.
#' @param SE The constant standard error for the mean lengths-at-age within a fish.  See details.
#' @param SECV The coefficient of variation for the mean lengths-at-age within a fish.  See details.
#' @param SD The constant standard deviation of individuals for the mean lengths-at-age within a fish.  See details.
#' @param SDCV The coefficient of variation for individuals for the mean lengths-at-age within a fish.  See details.
#' @param minAge The minimum possible age-at-capture for an individual.
#' @param maxAge The maximum possible age-at-capture for an individual.
#' @param agedist A vector of values that are proportional to the percentage of fish expected for each age.  This vector must have a length equal to \code{maxAge}-\code{minAge}+1.  See examples.
#' @param growth.per A numeric vector of length 2 that contains the Julian dates for the beginning and ending dates for the period when the fish grows.  See details and examples.
#' @param sample.per A numeric vector of length 1 or 2 that contains the Julian beginning and ending date(s) for sampling the fish.  If only one value is given, then it is assumed that sampling occurred only on that date.  See details and examples.
#' @param dataType The type of data frame that should be returned.  See details.
#' @param lendigs A single numeric that controls the number of digits for the length variable.
#' @param agedigs A single numeric that controls the number of digits for the fractional age variable.
#' @param backcalc A logical that indicates whether the lengths at previous ages should be considered to be from back-calculation (\code{TRUE}) or repeated observations (e.g., of tagged fish; \code{FALSE}, DEFAULT).
#' @param seed An integer that can be used to control the seed for the random number generator used in the function.  See details.
#'
#' @details This function can be used to generate random growth data, in one of four formats (see \sQuote{return} section and examples below) controlled by \code{dataType}.
#' 
#' The essential steps for constructing the simulated data are as follows:
#' \enumerate{
#'   \item Use \code{Linf}, \code{K}, and \code{t0} to set \sQuote{true} values for the three parameters in a typical VBGF (see \code{\link[FSA]{vbModels}} in \pkg{FSA} for the equation of this parameterization).
#'   \item For each fish, generate a random age-at-capture (i.e., number of completed growing seasons) between \code{minage} and \code{maxage} from the probability distribution for each age in \code{probdist}.
#'   \item Add a fractional age to the ages from the previous step and store in \code{ageFracG}.  By defaut, the \code{sample.per} is the first day of the \code{growth.per}, so \code{ageCap} and \code{ageFracG} will be the same.  However, if \code{sample.per} is not a single day equal to the first day of the \code{growth.per}, then \code{ageFracG} will be greater than \code{ageCap}.
#'   \item For each fish, generate random values for Linf, K, and t0 from a normal distribution with a mean of \code{Linf}, \code{K}, and \code{t0} and a standard deviation derived from each of \code{Linf}, \code{K}, and \code{t0}, respectively, times \code{paramCV}.  This step is used to model different parameter values for each fish.
#'   \item For each fish, generate the mean length at each \code{ageFracG} in the fish's life (i.e., if the fish is five years old, then a random length for ages one through five) according to the typical VBGF and the random parameters determined in the previous step for that fish.  These means are deterministic within a fish, but not among fish (because of the previous step).
#'   \item For each fish, add a random error to the mean lengths-at-age from the previous step.  These random errors are from a normal distribution with a mean of 0 and a standard error given by \code{SE} or \code{SECV} times the mean length-at-age.  This step is used to model random error in the mean lengths-at-age within an individual fish.
#'   \item For each fish, add a random error to the mean lengths-at-age from the previous step.  These random errors are from a normal distribution with a mean of 0 and a standard deviation given by \code{SD} or \code{SDCV} times the mean length-at-age.  This step is used to model random error of individual lengths around the mean lengths-at-age within an individual fish.
#'   \item Repeat the previous all but the first step for all fish to be simulated.
#'   \item Store the data for all fish in a data.frame, modify that data.frame to match the format required by \code{dataType}, and return the data.frame.
#' }
#' 
#' Note that \code{seed} can be used to control the random seed used by the function.  Use the same \code{seed} with different \code{dataType}s to return the same data in different formats.  See examples.
#' 
#' @return A data frame in one of four formats as determined by \code{dataType}.  The variables that may occur in each data.frame are:
#' \itemize{
#'   \item \code{id}: a unique identification for individual fish.
#'   \item \code{lenCap}: the observed length at the time of capture (corresponds to \code{ageCap}).
#'   \item \code{ageCap}: the number of completed growing seasons at the time of capture.
#'   \item \code{ageFracG}: the fractional age of the fish -- i.e., number of completed growing seasons plus a decimal representation of the fraction of the current growing seasons completed at the time of capture.
#'   \item \code{ageFracY}: the fractional age of the fish as a proportion of the sampling year -- i.e., number of completed growing seasons plus the Julian sampling date divided by 365 (assumes all years are 365 days).
#'   \item \code{agePrev}: an age in the fish's growth history.  Will be a whole number if \code{backcalc=TRUE} (i.e., length was back-calculated to an age that represents the beginning of a growing season) or a decimal if \code{backcalc=FALSE} (i.e., repeated observations of length may not be at the beginning of a growing season).
#'   \item \code{len}: the observed length that corresponds to \code{agePrev}.
#' }
#' 
#' Which of these variables is returned depends on the \code{dataTYpe} as follows:
#' \itemize{
#'   \item \sQuote{atCapture} format: Returns \code{id}, \code{ageCap}, \code{ageFracG}, and \code{lenCap}.  Each row corresponds to the observations for an individual fish at the time of capture.
#'   \item \sQuote{long} format: Returns \code{id}, \code{ageCap}, \code{agePrev}, \code{ageFracG}, \code{len}, and \code{lenCap}.  Each row corresponds to one age observation and individual fish are spread across as many rows as \code{ageCap}.  If \code{backcalc=TRUE} then there will be one \code{ageFracG} greater than \code{ageCap} (represents \dQuote{plus} or current season's growth) and all other \code{ageFracG}s will be whole numbers.
  #'   \item \sQuote{groupedData} format: Returns same as \sQuote{long} format but as a \code{\link[nlme]{groupedData}} object (from the \pkg{nlme} package).
  #'   \item \sQuote{wide} format: returns \code{id}, \code{ageCap}, \code{lenCap}, and lengths-at-age (\code{ageX}; both previous ages and age-at-capture).  In \sQuote{wide} format each row corresponds to a fish and may include several lengths (at previous ages and at-capture).
#' }
#' 
#' @author Derek H. Ogle, \email{dogle@@northland.edu}.  This function was motivated by the simulation methodology shown in Box 2 of Vigliola and Meekan (2009).
#'
#' @seealso See \code{\link[nlme]{groupedData}} from \pkg{nlme} for how to analyze data in \sQuote{groupedData} format.
#'
#' @references Vigliola, L. and M.G. Meekan.  2009.  The back-calculation of fish growth from otoliths.  Chapter 6 (pp. 174-211) in B.S. Green et al. (eds.), Tropical Fish Otoliths: Information for Assessment, Management and Ecology, Reviews: Methods and Technologies in Fish Biology and Fisheries 11.
#'
#' @keywords manip
#' 
#' @examples
#' ###########################################################
#' ## At-capture (observation) format examples
#' # Use CV as errors for means and individuals, sample at
#' #   beginning of growth season, uniform age distribution
#' vbcap1 <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,lendigs=1)
#' head(vbcap1)
#' plot(lenCap~ageCap,data=vbcap1,pch=19,col=rgb(0,0,0,1/5))
#' 
#' # Use constant errors for means and individuals, sample at
#' #   beginning of growth season, uniform age distribution
#' vbcap2 <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,SE=2,SD=0.5,lendigs=1)
#' head(vbcap2)
#' plot(lenCap~ageCap,data=vbcap2,pch=19,col=rgb(0,0,0,1/5))
#' xtabs(~ageCap,data=vbcap2)
#' 
#' # Use constant errors for means and individuals, sample at
#' #   beginning of growth season, non-uniform age distribution
#' #   (i.e., relatively few young and old fish)
#' ad <- c(1,1,2,4,4,4,2,1,1)
#' vbcap3 <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,agedist=ad,SECV=0.02,SDCV=0.04,lendigs=1)
#' head(vbcap3)
#' plot(lenCap~ageCap,data=vbcap3,pch=19,col=rgb(0,0,0,1/5))
#' xtabs(~ageCap,data=vbcap3)
#' 
#' # Use CV as errors for means and individuals, uniform age
#' #   distribution, sample throughout the growth season,
#' #   leads to fractional ages
#' # Assumed fish growth period (make a Julian date)
#' growth.per.date <- c("1-Apr-2010","21-Oct-2010")
#' ( growth.per <- strptime(growth.per.date,"%d-%b-%Y")$yday+1 )
#' # Assumed sampling period (make a Julian date)
#' sample.per.date <- c("15-Apr-2010","1-Jun-2010")
#' ( sample.per <- strptime(sample.per.date,"%d-%b-%Y")$yday+1 )
#' vbcap4 <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,lendigs=1,
#'                     growth.per=growth.per,sample.per=sample.per)
#' head(vbcap4)
#' plot(lenCap~ageCap,data=vbcap4,pch=19,col=rgb(0,0,0,1/5))
#' points(lenCap~ageFracG,data=vbcap4,pch=19,col=rgb(1,0,0,1/5))
#' points(lenCap~ageFracY,data=vbcap4,pch=19,col=rgb(0,1,0,1/5))
#' 
#' # Use CV as errors for means and individuals, sample at
#' #   beginning of growth season, uniform age distribution
#' #   but set a non-default minimum age
#' vbcap5 <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,minAge=3,SECV=0.02,SDCV=0.04,lendigs=1)
#' head(vbcap5)
#' plot(lenCap~ageCap,data=vbcap5,pch=19,col=rgb(0,0,0,1/5))
#' 
#' ###########################################################
#' ## Long (ungrouped format) ... repeated measures
#' # Use CV as errors for means and individuals, sample at
#' #   beginning of growth season, uniform age distribution
#' vblong1 <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,
#'                      dataType="long",lendigs=1)
#' head(vblong1,n=15)
#' 
#' # Use CV as errors for means and individuals, uniform age
#' #   distribution, sample throughout the growth season
#' #   (seasons set above), leads to fractional ages
#' vblong <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,lendigs=1,
#'                     growth.per=growth.per,sample.per=sample.per,dataType="long")
#' head(vblong,n=15)
#' 
#' ###########################################################
#' ## long (grouped format) ... repeated measures
#' # Use CV as errors for means and individuals, sample at
#' #   beginning of growth season, uniform age distribution
#' vbgrp1 <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,
#'                     dataType="groupedData",lendigs=1)
#' head(vbgrp1)
#' if (require(lattice)) plot(vbgrp1)
#' 
#' # Just errors in parameters (i.e., models unique parameters
#' # for each fish but no variability around the VBGF for an
#' # individual fish).
#' vbgrp2 <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SE=0,SD=0,dataType="groupedData",lendigs=1)
#' head(vbgrp2)
#' if (require(lattice)) plot(vbgrp2)
#'
#' ###########################################################
#' ## wide format ... repeated measures
#' # Use CV as errors for means and individuals, sample at
#' #   beginning of growth season, uniform age distribution
#' vbwide1 <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,
#'                      dataType="wide",lendigs=1)
#' head(vbwide1)
#' 
#' ###########################################################
#' ## Tag-recapture format data
#' # Use CV as errors for means and individuals, sample at
#' #   beginning of growth season, uniform age distribution
#' vbtag1 <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,minAge=3,SECV=0.02,SDCV=0.04,
#'                     lendigs=1,dataType="tagrecap")
#' head(vbtag1)
#' 
#' # Use CV as errors for means and individuals, uniform age
#' #   distribution, sample throughout the growth season,
#' #   leads to fractional ages ... more interesting
#' vbtag2 <- vbDataGen(100,Linf=30,K=0.2,t0=-0.2,minAge=3,SECV=0.02,SDCV=0.04,lendigs=1,
#'                     growth.per=growth.per,sample.per=sample.per,dataType="tagrecap")
#' head(vbtag2)
#' 
#' ###########################################################
#' # Back-calculated-like format
#' # Use CV as errors for means and individuals, uniform age
#' #   distribution, sample throughout the growth season,
#' #   leads to fractional ages
#' vbBClong1 <- vbDataGen(10,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,lendigs=1,
#'                        growth.per=growth.per,sample.per=sample.per,
#'                        dataType="long",backcalc=TRUE)
#' head(vbBClong1,n=15)
#' 
#' ###########################################################
#' ## Use seed to get same "fish" in different formats
#' sd <- 1534756   # set seed
#' # generate three types of data using same seed
#' ( vbcapA <- vbDataGen(3,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,lendigs=1,
#'                       growth.per=growth.per,sample.per=sample.per,seed=sd) )
#' ( vblongA <- vbDataGen(3,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,dataType="long",
#'                        growth.per=growth.per,sample.per=sample.per,lendigs=1,seed=sd) )
#' ( vbwideA <- vbDataGen(3,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,dataType="wide",
#'                        growth.per=growth.per,sample.per=sample.per,lendigs=1,seed=sd) )
#' ( vbBCA <- vbDataGen(3,Linf=30,K=0.2,t0=-0.2,SECV=0.02,SDCV=0.04,dataType="long",
#'                      growth.per=growth.per,sample.per=sample.per,lendigs=1,
#'                      seed=sd,backcalc=TRUE) )
#' # Note that back-calculated results do not maintain same fish measurements
#' 
#' @export
vbDataGen <- function(n,Linf,K,t0,paramCV=0.1,
                      SE=NULL,SECV=NULL,SD=NULL,SDCV=NULL,
                      minAge=1,maxAge=9,agedist=rep(1,maxAge-minAge+1),
                      growth.per=c(1,365),sample.per=1,
                      dataType=c("atCapture","long","groupedData","wide","tagrecap"),
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
  tmp <- iCreateFracAges(agePrev,growth.per,sample.per)
  ## Put ids and age variables in a data.frame
  d <- data.frame(id,ageCap,agePrev,ageFracG=round(tmp$ageFracG,agedigs),
                  ageFracY=round(tmp$ageFracY,agedigs))
  ## Handle if back-calculated
  if (backcalc) {
    d <- iHndlBackcalc(d,ages)
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
  tmp <- LVB(d$ageFracG,iLinf,iK,it0)
  # add some error to model for the mean ....
  if (!is.null(SE)) tmp <- tmp + rnorm(length(tmp),0,SE)
    else tmp <- tmp + rnorm(length(tmp),0,SECV*tmp)
  # ... and for the individual fish
  if (!is.null(SD)) tmp <- tmp + rnorm(length(tmp),0,SD)
    else tmp <- tmp + rnorm(length(tmp),0,SDCV*tmp)
  ## Find length-at-capture, repeat for each agePrev
#  plusgrowth <- ifelse(max(sample.per)==growth.per[1],FALSE,TRUE)
#  if (plusgrowth) lenCap <- round(rep(tmp[d$ageFracY>=d$ageCap],ages),lendigs)
#    else lenCap <- round(rep(tmp[d$ageCap==d$ageFracY],ages),lendigs)
  if (!backcalc) lenCap <- round(rep(tmp[d$ageCap==d$agePrev],ages),lendigs)
    else lenCap <- round(rep(tmp[is.na(d$agePrev)],ages),lendigs)
    ## put together as a data.frame
  d <- data.frame(d,lenCap,len=round(tmp,lendigs))
  ## Return result in format chosen by user (in dataType)
  if (dataType=="atCapture") d <- subset(d,ageCap==agePrev,c("id","lenCap","ageCap","ageFracG","ageFracY"))
    else if (dataType=="groupedData") d <- nlme::groupedData(len~agePrev|id,data=d,labels=list(x="Age",y="Length"))
      else if (dataType=="tagrecap") d <- iCreateTaggedDF(d)
        else if (dataType=="wide") {
          d$agePrev <- paste0("age",d$agePrev)
          d <- tidyr::spread(d[,-which(names(d) %in% c("ageFracG","ageFracY"))],agePrev,len)
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
  # Choose a sample date from within the sample period for each age
  smpl.date <- sample(sample.per[1]:sample.per[2],length(age),replace=TRUE)
  # Convert the sample date to a fraction of the year
  fracY <- smpl.date/365
  # Convert the sample date to a fraction of the growth season.
  # If sampled before growth starts then fraction is 0, if sampled
  # after growth ceases then fraction is 1, if sampled within the
  # growth season then convert to a fraction completed of the the
  # growth seasons
  fracG <- (smpl.date-growth.per[1])/(growth.per[2]-growth.per[1])
  fracG[fracG < 0] <- 0
  fracG[fracG > 1] <- 1
  # add fractions to ages and return
  list(ageFracY=age+fracY,ageFracG=age+fracG)
}

# ============================================================
# Create the data.frame of tagged fish
# ============================================================
iCreateTaggedDF <- function(d) {
  ageCap <- NULL # to avoide "global bindings" warning in rcmd check
  id <- NULL # to avoide "global bindings" warning in rcmd check
  # reduce to those fish that are at least age-2
  d <- subset(d,ageCap>=2)
  # get list of fish
  ids <- unique(d$id)
  # create a new matrix to store the data
  newd <- matrix(NA,nrow=length(ids),ncol=4)
  # get two ages from each fish and put in new matrix
  for (i in 1:length(ids)) {
    tmpdf <- subset(d,id==ids[i])
    tmp <- sample(1:nrow(tmpdf),2)
    tmp <- tmp[order(tmp)]
    newd[i,] <- c(ids[i],tmpdf$len[tmp[1]],tmpdf$len[tmp[2]],tmpdf$ageFracG[tmp[2]]-tmpdf$ageFracG[tmp[1]])
  }
  # turn matrix into a data.frame
  d <- data.frame(newd)
  names(d) <- c("id","Lm","Lr","atLarge")
  d
}

iHndlBackcalc <- function(d,ages) {
  # isolate all at-capture info
  tmp <- d[d$ageFracG>d$ageCap,]
  tmp$agePrev <- NA
  # reduce ageFracG to agePrev
  d$ageFracG <- d$agePrev
  # row bind at-capture info
  d <- rbind(d,tmp)
  d <- d[order(d$id,d$ageFracG),]
  rownames(d) <- 1:nrow(d)
  d
}
