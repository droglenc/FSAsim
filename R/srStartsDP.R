#' @title Uses dynamic plots to find reasonable starting values for parameters in specific paramaterizations of common stock-recruitment models.
#'
#' @description Uses dynamic plots to find reasonable starting values for parameters in specific parameterizations of the \dQuote{Beverton-Holt}, \dQuote{Ricker},  \dQuote{Shepherd}, or \dQuote{Saila-Lorda} stock-recruitment models. Use \code{srFunShow} described in \code{\link[FSA]{stockRecruitment}} for the equations of each model.
#'
#' @details Starting values can be obtained by plotting the data with the model superimposed but tied to slider bars for changing parameters. One can change the parameters until a reasonable fit is observed and then use those valeus as starting values. The initial parameters for the slider bars are the starting values constructed as described in \code{\link[FSA]{srStarts}}. The range for the sliders will have a minimum that is \code{min.prop} times the initial value and a maximum that is \code{max.mult} times the initial value. The step or interval of the slider bar is \code{delta.mult} times the initial value.
#'
#' @param formula A formula of the form \code{Recruits~Stock}.
#' @param data A data frame in which \code{Recruits} and \code{Stock} are found.
#' @param type A string that indicates the type of the stock-recruitment model. Must be one of \code{"BevertonHolt"}, \code{"Ricker"}, \code{"Shepherd"}, or \code{"SailaLorda"}.
#' @param param A numeric that indicates the parameterization of the stock-recruitment model type. This is ignored if \code{type="Shepherd"} or \code{type="SailaLorda"}
#' @param min.prop A single numeric that is used to set the minimum values for the slider bars in the dynamic plot. see details.
#' @param max.mult A single numeric that is used to set the maximum values for the slider bars in the dynamic plot. see details.
#' @param delta.prop A single numeric that is used to set the step value for the slider bars in the dynamic plots. See details.
#' @param \dots Further arguments passed to the methods.
#'
#' @return Nothing, but a dynamic plot is created.
#'
#' @author Derek H. Ogle, \email{derek@@derekogle.com}
#'
#' @section IFAR Chapter: 13-Recruitment.
#'
#' @seealso See \code{\link[FSA]{srStarts}} for related functionality.
#'
#' @references Ogle, D.H. 2016. \href{http://derekogle.com/IFAR}{Introductory Fisheries Analyses with R}. Chapman & Hall/CRC, Boca Raton, FL.
#'
#' @keywords manip
#'
#' @examples
#' ## Dynamic Plots Method -- ONLY RUN IN INTERACTIVE MODE
#' if (interactive()) {
#'   require(FSA)
#'   data(CodNorwegian)
#'   # Beverton-Holt
#'   srStartsDP(recruits~stock,data=CodNorwegian)
#'   srStartsDP(recruits~stock,data=CodNorwegian,param=2)
#'   srStartsDP(recruits~stock,data=CodNorwegian,param=3)
#'   srStartsDP(recruits~stock,data=CodNorwegian,param=4)
#'   # Ricker Models
#'   srStartsDP(recruits~stock,data=CodNorwegian,type="Ricker")
#'   srStartsDP(recruits~stock,data=CodNorwegian,type="Ricker",param=2)
#'   srStartsDP(recruits~stock,data=CodNorwegian,type="Ricker",param=3)
#'   # Shepherd, Saila-Lorda, and Independence Models
#'   srStartsDP(recruits~stock,data=CodNorwegian,type="Shepherd")
#'   srStartsDP(recruits~stock,data=CodNorwegian,type="SailaLorda")
#'   srStartsDP(recruits~stock,data=CodNorwegian,type="independence")
#' } ## END .. ONLY INTERACTIVE
#'
#' @rdname srStartsDP
#' @export
srStartsDP <- function(formula,data=NULL,
                       type=c("BevertonHolt","Ricker","Shepherd",
                              "SailaLorda","independence"),
                       param=1,min.prop=0.1,max.mult=3,delta.prop=0.005,...) {
  ## some checks
  type <- match.arg(type)
  ## get the R and S vectors
  tmp <- FSA:::iHndlFormula(formula,data,expNumR=1,expNumE=1)
  if (!tmp$metExpNumR)
    FSA:::STOP("'formula' must have only one LHS variable.")
  if (!tmp$Rclass %in% c("numeric","integer"))
    FSA:::STOP("LHS variable must be numeric.")
  if (!tmp$metExpNumE)
    FSA:::STOP("'formula' must have only one RHS variable.")
  if (!tmp$Eclass %in% c("numeric","integer")) 
    FSA:::STOP("RHS variable must be numeric.")
  R <- tmp$mf[,tmp$Rname[1]]
  S <- tmp$mf[,tmp$Enames[1]]
  ## find starting values for the type and param
  sv <- FSA::srStarts(formula,data,type,param)
  if (!iCheckRStudio()) FSA:::STOP("'srStartsDP' only works in RStudio.")
  if (iChk4Namespace("manipulate")) {
    iSRStartsDynPlot(S,R,type,param,sv,min.prop,max.mult,delta.prop)
  }
}

##############################################################
# INTERNAL FUNCTIONS
##############################################################
#=============================================================
# Dynamics plots for finding starting values -- main function
#=============================================================
iSRStartsDynPlot <- function(S,R,type,param,sv,min.prop,max.mult,delta.prop) {
  p1 <- p2 <- p3 <- NULL
  ## Set the minimum value as a proportion (from min.prop) of the starting
  ## values, the maximum value at a multiple (from max.mult) of the starting
  ## values, and the delta value at a proportion (from delta.prop) of the
  ## starting values. Unlist first to make as a vector.
  sl.defaults <- unlist(sv)
  sl.mins <- sl.defaults*min.prop
  sl.maxs <- sl.defaults*max.mult
  sl.deltas <- sl.defaults*delta.prop
  ## Grab names from the sv vector
  sl.names <- names(sl.defaults)
  ## Set up names that are specific to type and param
  if (type=="independence") {  ## just one slider
    manipulate::manipulate(
      { iSRDynPlot(S,R,type,param,p1,p2,p3) },
      p1=manipulate::slider(sl.mins[[1]],sl.maxs[[1]],step=sl.deltas[[1]],
                            initial=sl.defaults[[1]],label=sl.names[[1]])
    )
  } else if (type %in% c("BevertonHolt","Ricker")) {  ## two sliders
    manipulate::manipulate(
      { iSRDynPlot(S,R,type,param,p1,p2,p3) },
      p1=manipulate::slider(sl.mins[[1]],sl.maxs[[1]],step=sl.deltas[[1]],
                            initial=sl.defaults[[1]],label=sl.names[[1]]),
      p2=manipulate::slider(sl.mins[[2]],sl.maxs[[2]],step=sl.deltas[[2]],
                            initial=sl.defaults[[2]],label=sl.names[[2]])
    )
  } else {  ## three sliders
    manipulate::manipulate(
      { iSRDynPlot(S,R,type,param,p1,p2,p3) },
      p1=manipulate::slider(sl.mins[[1]],sl.maxs[[1]],step=sl.deltas[[1]],
                            initial=sl.defaults[[1]],label=sl.names[[1]]),
      p2=manipulate::slider(sl.mins[[2]],sl.maxs[[2]],step=sl.deltas[[2]],
                            initial=sl.defaults[[2]],label=sl.names[[2]]),
      p3=manipulate::slider(sl.mins[[3]],sl.maxs[[3]],step=sl.deltas[[3]],
                            initial=sl.defaults[[3]],label=sl.names[[3]])
    )
  }
}

#===============================================================================
# Constructs the actual plot in the dynamics plots for finding starting values
#===============================================================================
iSRDynPlot <- function(S,R,type,param,p1,p2=NULL,p3=NULL) {
  ## create a sequence of spawning stock values
  max.S <- max(S,na.rm=TRUE)
  x <- c(seq(0,0.2*max.S,length.out=500),seq(0.2*max.S,max.S,length.out=500)) 
  ## create a sequence of matching recruit values according to
  ## the model type and param
  y <- iSRDynPlot_makeR(x,type,param,p1,p2,p3)
  ## Construct the scatterplot with superimposed model
  withr::local_par(list(mfrow=c(1,2),mar=c(3.5,3.5,1.25,1.25),
                        mgp=c(2,0.4,0),tcl=-0.2,pch=19))
  graphics::plot(S,R,xlab="Parental (Spawner) Stock",ylab="Recruits",
       ylim=c(0,max(R,na.rm=TRUE)),xlim=c(0,max(S,na.rm=TRUE)))
  graphics::lines(x,y,lwd=2,lty=1,col="blue")
  graphics::plot(S,R/S,xlab="Parental (Spawner) Stock",ylab="Recruits/Spawner",
       ylim=c(0,max(R/S,na.rm=TRUE)),xlim=c(0,max(S,na.rm=TRUE)))
  graphics::lines(x,y/x,lwd=2,lty=1,col="blue")
}

#===============================================================================
# Construct values for R given values of S, a stock-recruit model and
# parameterization, and values of the parameters. Called by iSRDynPlot()
#===============================================================================
iSRDynPlot_makeR <- function(x,type,param,p1,p2,p3=NULL){
  if (type=="BevertonHolt") {
    if (param==1) {
      # p1=a,p2=b
      y <- p1*x/(1+p2*x)
    } else if (param==2) {
      # p1=a, p2=Rp
      y <- p1*x/(1+p1*(x/p2))
    } else if (param==3) {
      # p1=atilde, p2=btilde
      y <- x/(p1+p2*x)
    } else if (param==4) {
      # p1=atilde, p2=Rp
      y <- x/(p1+x/p2)
    } # end B-H switch
  } else if (type=="Ricker") {
    if (param==1) {
      # p1=a,p2=b
      y <- p1*x*exp(-p2*x)
    } else if (param==2) {
      # p1=atilde, p2=b
      y <- x*exp(p1-p2*x)
    } else if (param==3) {
      # p1=a, p2=Rp
      y <- p1*x*exp(-p1*x/(p2*exp(1)))
    } # end Ricker switch
  } else if (type=="Shepherd") {
      # p1=a, p2=b, p3=c
      y <- (p1*x)/(1+(p2*x)^p3)
  } else if (type=="SailaLorda") {
      # p1=a, p2=b, p3=c    
      y <- p1*(x^p3)*exp(-p2*x) 
  } else {  # Independence model
      # p1=a
      y <- p1*x
  }
  # return the value of y
  y
}