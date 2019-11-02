#' @title Internal functions used in FSAsim.
#'
#' @description Internal functions used in FSAsim
#'
#' @rdname FSA-internals
#' @keywords internal
#' @aliases getVarFromFormula iChk4Namespace iCheckRStudio
#'

##################################################################
### Get variable names from formulas 
##################################################################
getVarFromFormula <- function(formula,data,expNumVars=NULL) {
  varNms <- names(stats::model.frame(formula,data=data))
  # don't "error" check the number of variables
  if (is.null(expNumVars)) varNms
  else if (length(varNms)!=expNumVars)
    FSA:::STOP("Function only works with formulas with ",expNumVars," variable",
               ifelse(expNumVars==1,".","s."))
  else varNms
}

##################################################################
# Check if a required namespace can be loaded and, if not, send
#   an error message.
##################################################################
iChk4Namespace <- function(pkg) {
  res <- requireNamespace(pkg,quietly=TRUE)
  if (!res) FSA:::STOP(paste0("The '",pkg," package must be installed."))
  res
}


##################################################################
# Try to determine if RStudio is being used.
##################################################################
iCheckRStudio <- function () "tools:rstudio" %in% search()
