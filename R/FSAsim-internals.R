#' @title Internal functions used in FSAsim.
#'
#' @description Internal functions used in FSAsim
#'
#' @rdname FSA-internals
#' @keywords internal
#' @aliases .onAttach getVarFromFormula iChk4Namespace iCheckRStudio
#'

##################################################################
## Sends a start-up message to the console when the package is loaded.
##################################################################
.onAttach <- function(lib,pkg,...) {
  ## Get version number -- basically code from readPkgVersion in SweaveListingUtils
  vers <- read.dcf(system.file("DESCRIPTION",package=pkg,lib.loc=lib),fields="Version")
  ## Send message
  msg <- paste("\n\n")
  msg <- paste(msg,"############################################\n")
  msg <- paste(msg,"##","FSAsim package, version",vers,"         ##\n")
  msg <- paste(msg,"##   by Derek H. Ogle, Northland College  ##\n")
  msg <- paste(msg,"##                                        ##\n")
  msg <- paste(msg,"## Type ?FSAsim for documentation.        ##\n")
  msg <- paste(msg,"## Type citation('FSAsim') for citation   ##\n")
  msg <- paste(msg,"##   please cite if used in publication.  ##\n")
  msg <- paste(msg,"############################################\n\n")
  packageStartupMessage(msg)
}


##################################################################
### Get variable names from formulas 
##################################################################
getVarFromFormula <- function(formula,data,expNumVars=NULL) {
  varNms <- names(stats::model.frame(formula,data=data))
  # don't "error" check the number of variables
  if (is.null(expNumVars)) varNms
  else if (length(varNms)!=expNumVars) stop("Function only works with formulas with ",expNumVars," variable",ifelse(expNumVars==1,".","s."))
  else varNms
}

##################################################################
# Check if a required namespace can be loaded and, if not, send
#   an error message.
##################################################################
iChk4Namespace <- function(pkg) {
  res <- requireNamespace(pkg,quietly=TRUE)
  if (!res) stop(paste0("The '",pkg," package must be installed."))
  res
}


##################################################################
# Try to determine if RStudio is being used.
##################################################################
iCheckRStudio <- function () "tools:rstudio" %in% search()
