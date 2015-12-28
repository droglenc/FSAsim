#' @title Dynamic plot to explore numbers in a cohort over time.
#'
#' @description Constructs a plot of the hypothetical number of individuals in a cohort over time.  The initial size of the cohort and the instantaneous mortality rate are controlled by slider bars.
#'
#' @param age.max A single numeric that indicates the maximum age to use in the simulations.
#'
#' @return None.  An interactive graphic connected to slider controls is produced.
#' 
#' @author Derek H. Ogle, \email{dogle@@northland.edu}
#' 
#' @keywords iplot
#' 
#' @examples
#' if (interactive()) {
#'
#' cohortSim()
#'
#' } # end if interactive
#' 
#' @export
#'
cohortSim <- function(age.max=15) {
  # Trying to deal with no visible bindings problem
  No <- Z <- NULL
  if (!iCheckRStudio()) stop("'cohortSim' only works in RStudio.",call.=FALSE)
  if (iChk4Namespace("manipulate")) {
    manipulate::manipulate(
      {
        Nt <- No*exp(-Z*1:age.max)
        op <- graphics::par(mar=c(3.5,3.5,1,1),mgp=c(2,0.4,0),tcl=-0.2)
        graphics::plot(0:age.max,c(No,Nt),type="l",lwd=2,col="blue",
                       xlab="Time / Age",ylab="Population Size",ylim=c(0,10000))
        graphics::par(op)
      },
      No=manipulate::slider(1000,10000,step=500,initial=10000,label="Initial Numbers (No)"),
      Z=manipulate::slider(0.05,1,step=0.05,initial=0.40,label="Instantaneous Mortality (Z)")
    )
  }
}
