#'Sample fish for construction of an age-length key.
#'
#'Samples fish from a data frame of lengths and ages to simulate selecton of fish
#'to be used to construct an age-length key.  Fish can either be selected as a
#'fixed number per length category (e.g., 5 fish per 10-mm length category) or as 
#'a random sample from the population.
#'
#'@param formula A formula of the form \code{age~length} where \dQuote{age}
#'generically represents a variable in \code{data} that contain the observed
#'ages and \dQuote{length} generically represents a variable in \code{data} that
#'contains the observed length measurements.
#'@param data A data.frame that minimally contains the observed age and length
#'measurements in the \code{formula}.
#'@param type A string indicating the type of sampling to conduct.  See details.
#'@param n.fixed A single numeric that indicates the number of fish to sample per
#'length category when \code{type="fixed"}.
#'@param n.random A single numeric that indicates the total number of fish to sample
#'when \code{type="random"}.
#'@param w A single numeric that indicates the width of length categories to create
#' when \code{type="fixed"}.
#'@param breaks A numeric vector of lower values for the break points of the
#'length categories when \code{type="fixed"}.
#'@param startcat A single numeric that indicates the beginning of the first
#'length category when \code{type="fixed"}.
#'@return A data.frame similar to \code{data} but where some of the ages have been
#'replaced with \code{NA} to simulate fish that have NOT been sampled for ages.
#'@keywords misc
#'@export
#'@examples
#'## get FSA library for lencat
#'require(FSA)
#'
#'## set seed for repeatability
#'set.seed(5234734)
#'
#'## Simulated individual ages (random) -- see simAges functions
#'bg.ages <- simAges(N0=500,A=0.35)
#'
#'## Simulated lengths, given the above ages
#'bg.lens <- simLenFromAge(bg.ages,228,0.206,0,8)
#'
#'## Combine lengths and ages into one data.frame
#'df <- data.frame(age=bg.ages,len=bg.lens)
#'# length frequency (10-mm) of the "population"
#'df$LCat <- lencat(df$len,w=10)
#'( popLF <- xtabs(~LCat,df) )
#'round(prop.table(popLF),2)
#'
#'## Simulate some fish having ages -- random selection of 200 individuals
#'dfage1 <- sample4ALK(age~len,data=df,type="random",n.random=200)
#'# length frequency (10-mm) of the "aged" fish
#'( smpl1LF <- xtabs(~LCat,Subset(dfage1,!is.na(age))) )
#'round(prop.table(smpl1LF),2)
#'
#'## Simulate some fish having ages -- 10 fish per 10-mm length category
#'dfage2 <- sample4ALK(age~len,data=df,type="fixed",w=10,n.fixed=10)
#'# length frequency (10-mm) of the "aged" fish
#'( smpl2LF <- xtabs(~LCat,Subset(dfage2,!is.na(age))) )
#'round(prop.table(smpl2LF),2)
#'
#'## Example age-length key from the last sample
#'# get age-sample
#'df.age <- Subset(dfage2,!is.na(age))
#'al.raw <- xtabs(~LCat+age,data=df.age)
#'al.key <- prop.table(al.raw,margin=1)
#'round(al.key,2)
#'
sample4ALK <- function(formula,data,type=c("fixed","random"),
                       n.fixed=10,n.random=100,
                       w=1,startcat=NULL,breaks=NULL) {
  ## Internal function to generate a fixed number from a group
  fsample <- function(x,nwant) {
    ntmp <- length(x)
    ntrue <- min(ntmp,nwant)
    nfalse <- max(0,ntmp-nwant)
    tmp <- sample(c(rep(TRUE,ntrue),rep(FALSE,nfalse)))
    as.logical(tmp)
  } ## end internal function
  ## Start main function
  type <- match.arg(type)
  cl <- getVarFromFormula(formula,data)
  if (length(cl)!=2) stop("Formulas must be of 'age~len'.",call.=FALSE)
  else {
    ca <- cl[1]
    cl <- cl[2]
  }
  if (type=="random") {
    # randomly sample for age by putting NA in for age a random selection of
    #   those NOT to age
    n <- nrow(data)
    notAged <- sample(1:n,n-n.random)
  } else {
    LCat <- lencat(df[,cl],w=w,breaks=breaks,startcat=startcat)
    notAged <- !as.logical(ave(LCat,LCat,FUN=function(x) fsample(x,nwant=n.fixed)))
  }
  data[notAged,ca] <- NA
  data
}
