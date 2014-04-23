#'Simulate an age structure.
#'
#'Constructs either a simulated age frequency table or individuals that follow
#'a simulated age frequency table.  The number of individuals at each age is
#'determined from an initial cohort size and an instantaneous or annual
#'mortality rate supplied by the user.
#'
#'The heart of this simulation is the assumption of an exponential decay model
#'that starts with \code{N0} individuals at age-0 and estimates Nt individuals
#'at each age t, from 0 or 1 (depending on the \code{incl.zero} value) to
#'\code{max.age}, with Nt=N0*exp(-Zt).  If \code{indivs=FALSE} then the ages
#'(t) and frequencies (Nt) are output as a data frame.  If \code{indivs=TRUE}
#'(default) and \code{rand=FALSE} then the age frequency values are used to
#'deterministically produce ages for individual fish -- i.e., essentially
#'\code{rep(t,Nt)}.  Finally, if \code{indivs=TRUE} (default) and
#'\code{rand=TRUE} (default) then the age frequency values are used to
#'stochastically produce ages for individual fish -- i.e., essentially
#'\code{sample(t,sum(Nt),replace=TRUE,prob=Nt/sum(Nt))}.
#'
#'@param N0 A numeric representing the initial cohort size to be used in the simulations.
#'@param A A numeric representing the annual mortality rate to be used in the simulations.
#'@param Z A numeric representing the instantaneous mortality rate to be used
#'in the simulations.  Will default to the value that corresponds to the user
#'supplied value for \code{A}.
#'@param incl.zero A logical indicating whether age-0 should be included in the
#'ouptut (\code{TRUE}) or not (\code{FALSE}; default).
#'@param max.age A numeric indicating the maximum age to use in the
#'simulations.  If \code{NULL} then the maximum age will be the age at which
#'the projected cohort size is approximately equal to 1.
#'@param indivs A logical indicating whether the output result should be ages
#'assigned to individual fish (\code{TRUE}; default) or an age-frequency table.
#'@param rand A logical indicating whether the individual ages should be chosen
#'randomly (\code{TRUE}; default) or deterministically.
#'@return Either a data frame containing the ages and the simulated number of
#'fish at each age (i.e., when \code{indivs=FALSE}) or a vector containing the
#'simulated ages for individual fish (i.e., when \code{indivs=TRUE}; default).
#'@seealso \code{\link{simLenFromAge}}, \code{\link{simLenSelectP}},
#'\code{\link{simLenSelectM}}
#'@keywords misc
#'@export
#'@examples
#'## set seed for repeatability
#'set.seed(5234734)
#'
#'## Simulated individual ages (random)
#'bg.ages <- simAges(N0=500,A=0.35)
#'str(bg.ages)
#'table(bg.ages)
#'
#'## Simulated age frequency (random)
#'bg.ages.sum <- simAges(N0=500,A=0.35,indivs=FALSE) 
#'str(bg.ages.sum)
#'t(bg.ages.sum)
#'
#'## Simulated individual ages (non-random)
#'bg.ages2 <- simAges(N0=500,A=0.35,rand=FALSE) 
#'table(bg.ages2)
#'
simAges <- function(N0=500,A=0.3,Z=-log(1-A),incl.zero=FALSE,max.age=NULL,indivs=TRUE,rand=TRUE) {
 # some error catching
  if (length(N0)>1 | missing(N0)) stop("'N0' must contain one and only one value",call.=FALSE)
  if (N0<0) stop("'N0' must be a postive number",call.=FALSE)
  if (length(A)>1 | missing(A)) stop("'A' must contain one and only one value",call.=FALSE)
  if (A<0 | A>1) stop("'A' must be a postive number less than 1",call.=FALSE)
  if (length(Z)>1) stop("'Z' must contain one and only one value",call.=FALSE)
  if (Z<0) stop("'Z' must be a postive number",call.=FALSE)

 # maximum age is the age when Nt falls to (approx.) 1, if none is given
  if (is.null(max.age)) max.age <- floor(-log(1/N0)/Z) 
  if (max.age<0) stop("'max.age' must be a postive number",call.=FALSE)
  if (!incl.zero & max.age<1) stop("'max.age' must be >1 if 'incl.zero=FALSE'",call.=FALSE)

 # vector of ages to model
  ifelse(incl.zero,ages <- 0:max.age,ages <- 1:max.age)
 # construct population sizes from exponential growth model
  Nt <- round(N0*exp(-Z*ages),0)
 # decide what to return 
  if (!indivs) data.frame(age=ages,num=Nt)
    else {
      if (!rand) age.dat <- rep(ages,Nt)
        else age.dat <- sample(ages,sum(Nt),replace=TRUE,prob=Nt/sum(Nt))
      age.dat
    }
}
