#' @title Evaluates and plots the traditional von Bertalanffy growth model for several sets of parameter choices.
#'
#' @description Evaluates and plots the traditional von Bertalanffy growth model for several parameter choices.  Typical use will be for comparing the shape of the model among parameter estimates for multiple species.
#'
#' @details The data frame, \code{df}, must have a column that contains the species name and columns for the Linf, K, and to parameters (with the columns having those names).  Each row contains values for a seperate model evaluation.  Typical exercises will select rows from the \code{data(VBGMlit)} data frame.
#'
#' All evaluations will be over ages that begin at zero and continue to the value or values provided in \code{ages}.  If \code{ages} has a length of 1 then all evaluations will extend from zero to this age.  However, if \code{ages} is a vector of the same length as the number of rows in \code{df} then the extent of the evaluations for each row in \code{df} will be from 0 to the corresponding value in \code{ages}.  Thus, \code{ages} can be used to provide model evaluations over different ages.
#'
#' @param df A data.frame with rows of different parameter choices.  See details.
#' @param ages A numeric or vector of the maximum age over which to evaluate the models.  See details.
#' @param \dots Additional arguments to pass to \code{plot}.
#'
#' @return None.  However a plot of the hypothetical lengths-at-age for each species shown in the rows of \code{df} is produced.
#'
#' @author Derek H. Ogle, \email{dogle@@northland.edu}
#'
#' @seealso \code{\link{VBGMlit}}.
#'
#' @keywords manip hplot
#'
#' @examples
#' # load VBGMlit data frame
#' data(VBGMlit)
#' # select three example species from VBGMlit
#' v <- VBGMlit[c(1,7,57),]
#' # default behavior for three example species
#' vbComp(v)
#' # same three species but with varying maximum ages
#' vbComp(v,ages=c(10,20,25))
#'
#' @export
#'
vbComp <- function(df,ages=20,...) {
  # put given data into a temporary df
  tmpdf <- df
  # Number of unique simulations
  num.sims <- dim(tmpdf)[1]
  # If only one age then create vector
  if (num.sims>1 & length(ages)==1) ages <- rep(ages,num.sims)
  # Find maximum age
  max.age <- max(ages)
  # Put a 0 in for to if it does not exist
  tmpdf$t0[which(is.na(tmpdf$t0))] <- 0
  # Set up matrix to hold sim results
  simdf <- matrix(NA,nrow=(max.age+1),ncol=num.sims)
  # Loop over species in tmpdf
  for (i in 1:num.sims) {
    #   Ages over which to model
    simages <- 0:ages[i]
    #   hypothetical length-at-age data
    simdf[simages+1,i] <- with(tmpdf,Linf[i]*(1-exp(-K[i]*(simages-t0[i]))))
  }
  # Modify plot parameters
  op <- graphics::par(mar=c(3.5,3.5,1,1),mgp=c(2,0.75,0))
  on.exit(graphics::par(op))
  # plot results
  graphics::matplot(0:max.age,simdf[,1:num.sims],type="b",lwd=3,
                    xlab="Age",ylab="Mean Length",...)
  # Create a "key" to the graph to be
  key <- cbind(c(1:9,0,letters)[1:num.sims],df,ages)
  #   printed in the console
  colnames(key) <- c("Sim","Species","Linf","K","t0","Max. Age")
  # Print the legend/key with graph
  print(key)
}
