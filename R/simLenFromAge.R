#'Simulate fish length from given ages.
#'
#'Constructs simulated fish lengths from a set of given ages and parameters for
#'a von Bertalanffy growth model supplied by the user.
#'
#'This simulation simulates fish lengths by first predicing the mean
#'length-at-age given a fish's age using the von Bertalanffy growth model and
#'then randomly selecting a length from a normal distribution given this mean
#'length-at-age and the supplied value for sigma.
#'
#'@param ages A vector containg ages of individual fish.
#'@param Linf A numeric representing the asyptotic mean length (L_infinity) in
#'the von Bertalanffy growth model.
#'@param K A numeric representing the Brody growth coefficient (K) in the von
#'Bertalanffy growth model.
#'@param t0 A numeric representing the time when the mean length is zero (t_0)
#'in the von Bertalanffy growth model.
#'@param sigma A numeric representing the standard deviation (i.e., individual
#'error) around the von Bertalanffy growth model.
#'@param additive A logical indicating whether the standard deviation is for
#'the additive- (\code{TRUE}; default) or multiplicative-error version of the
#'von Bertalanffy growth model.
#'@param digits A numeric controlling the number of digits to which the length
#'data should be rounded before returning.
#'@return A vector containing the simulated lengths for individual fish.
#'@seealso \code{\link{simAges}}, \code{\link{simLenSelectP}},
#'\code{\link{simLenSelectM}}
#'@export
#'@keywords misc
#'@examples
#'## Load FSA package for Summarize()
#'library(FSA)
#'
#'## set seed for repeatability
#'set.seed(5234734)
#'
#'## Simulated individual ages (random)
#'#    see simAges functions
#'bg.ages <- simAges(N0=500,A=0.35)
#'
#'## Simulated lengths, given the above ages
#'bg.lens <- simLenFromAge(bg.ages,228,0.206,0,8)
#'
#'## Some summaries
#'df <- data.frame(age=bg.ages,len=bg.lens)
#'Summarize(len~age,data=df,digits=1)
#'plot(len~age,data=df,pch=16,col=rgb(0,0,0,0.1),xlab="Age",ylab="Length")
#'
simLenFromAge <- function(ages,Linf,K,t0,sigma,additive=TRUE,digits=0) {
 # some error catching
  if (length(Linf)>1 | missing(Linf)) stop("'Linf' must contain one and only one value",call.=FALSE)
  if (Linf<0) stop("'Linf' must be a postive number",call.=FALSE)
  if (length(K)>1 | missing(K)) stop("'K' must contain one and only one value",call.=FALSE)
  if (K<0) stop("'K' must be a postive number",call.=FALSE)
  if (length(t0)>1 | missing(t0)) stop("'t0' must contain one and only one value",call.=FALSE)
  if (length(sigma)>1 | missing(sigma)) stop("'sigma' must contain one and only one value",call.=FALSE)
  if (sigma<0) stop("'sigma' must be a postive number",call.=FALSE)

 # simulate mean length-at-age using a von Bert model
  mns <- Linf*(1-exp(-K*(ages-t0)))
  if (additive) {
   # simulate lengths assuming additive normal deviates
    lens <- rnorm(length(mns),mean=mns,sd=sigma) 
  } else {
   # simulate lengths assuming multiplicative normal deviates
    lens <- exp(rnorm(length(mns),mean=log(mns),sd=sigma))
  }
  round(lens,digits)
}
