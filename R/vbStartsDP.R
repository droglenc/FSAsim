#' @title Uses dynamic plots to find reasonable starting values for a von Bertalanffy growth function.
#' 
#' @description Uses dynamic plots to find reasonable starting values for the parameters in a specific parameterization of the von Bertalanffy growth function.
#' 
#' @details Starting values can be obtained by plotting the data with the model superimposed but tied to slider bars that allow the parameters to be interactively changed.  One can change the parameters until a reasonable fit is observed and then use those valeus as starting values.  The initial parameters for the slider bars are the starting values constructed with \code{\link[FSA]{vbStarts}}.  It should be noted that the dynamic plot may show an error of \dQuote{[tcl] can't get device image}, but the plot will correctly update if the slider bar is adjusted.
#' 
#' @param formula A formula of the form \code{len~age}.
#' @param data A data frame that contains the variables in \code{formula}.
#' @param type A string that indicates the parameterization of the von Bertalanffy model.
#' @param ages2use A numerical vector of the two ages to be used in the Schnute or Francis paramaterizations.  See \code{\link[FSA]{vbStarts}}.
#' @param methEV A string that indicates how the lengths of the two ages in the Schnute paramaterization or the three ages in the Francis paramaterization should be derived.  See \code{\link[FSA]{vbStarts}}.
#' @param meth0 A string that indicates how the t0 and L0 paramaters should be derived.  See \code{\link[FSA]{vbStarts}}.
#' @param \dots Further arguments passed to the methods.
#' 
#' @return None, but a dynamic plot is constructed.
#' 
#' @note The \sQuote{original} and \sQuote{vonBertalanffy} and the \sQuote{typical} and \sQuote{BevertonHolt} parameterizations are synonymous.
#' 
#' @author Derek H. Ogle, \email{derek@@derekogle.com}
#' 
#' @section IFAR Chapter: 12-Individual Growth.
#' 
#' @seealso See \code{\link[FSA]{vbStarts}} for related functionality.
#' 
#' @references Ogle, D.H.  2016.  \href{http://derekogle.com/IFAR}{Introductory Fisheries Analyses with R}.  Chapman & Hall/CRC, Boca Raton, FL.
#' 
#' See references in \code{\link{vbFuns}}.
#' 
#' @keywords manip
#' @examples
#' ## Dynamic Plots Method -- ONLY RUN IN INTERACTIVE MODE
#' if (interactive()) {
#'   require(FSA)
#'   data(SpotVA1)
#'   vbStartsDP(tl~age,data=SpotVA1)
#'   vbStartsDP(tl~age,data=SpotVA1,type="original")
#'   vbStartsDP(tl~age,data=SpotVA1,type="GQ")
#'   vbStartsDP(tl~age,data=SpotVA1,type="Mooij")
#'   vbStartsDP(tl~age,data=SpotVA1,type="Weisberg")
#'   vbStartsDP(tl~age,data=SpotVA1,type="Francis",ages2use=c(0,5))
#'   vbStartsDP(tl~age,data=SpotVA1,type="Schnute",ages2use=c(0,5))
#' } 
#' 
#' @export vbStartsDP
vbStartsDP <- function(formula,data=NULL,
                       type=c("Typical","typical","BevertonHolt",
                              "Original","original","vonBertalanffy",
                              "GQ","GallucciQuinn","Mooij","Weisberg",
                              "Schnute","Francis","Somers","Somers2"),
                       ages2use=NULL,methEV=c("means","poly"),meth0=c("yngAge","poly"),
                       ...) {
  ## some checks of arguments
  type <- match.arg(type)
  type <- FSA::capFirst(type)
  if (type=="BevertonHolt") type <- "Typical"
  if (type=="vonBertalanffy") type <- "Original"
  methEV <- match.arg(methEV)
  meth0 <- match.arg(meth0)
  ## get the length and age vectors
  tmp <- FSA:::iHndlFormula(formula,data,expNumR=1,expNumE=1)
  len <- tmp$mf[,tmp$Rname[1]]
  age <- tmp$mf[,tmp$Enames[1]]
  ## Get starting values
  sv <- vbStarts(formula,data,type,ages2use,methEV,meth0)
  if (!requireNamespace("relax")) stop("'vbStarts' requires the 'relax' package to be installed to construct the dynamic plot.",call.=FALSE)
  else iVBStartsDynPlot(age,len,type,sv,ages2use)
}



##############################################################
# INTERNAL FUNCTIONS
##############################################################
#=============================================================
# Dynamics plots for finding starting values -- main function
#=============================================================
iVBStartsDynPlot <- function(age,len,type,sv,ages2use) {
  ## internal refresh function for the dialog box
  refresh <- function(...) {
    p1 <- relax::slider(no=1)
    p2 <- relax::slider(no=2)
    p3 <- relax::slider(no=3)
    iVBDynPlot(age,len,type,p1,p2,p3,ages2use)
  } # end internal refresh

  ## internal function to make minimum values for the sliders
  iMake.slMins <- function(sv) {
    svnms <- names(sv)
    tmp <- c(Linf=NA,K=NA,t0=NA,L0=NA,omega=NA,t50=NA,L1=NA,L2=NA,L3=NA)
    if ("Linf" %in% svnms) tmp["Linf"] <- 0.5*sv[["Linf"]]
    if ("K" %in% svnms) tmp["K"] <- 0.01
    if ("t0" %in% svnms) tmp["t0"] <- -5
    if ("L0" %in% svnms) tmp["L0"] <- 0
    if ("omega" %in% svnms) tmp["omega"] <- 0.5*sv[["omega"]]
    if ("t50" %in% svnms) tmp["t50"] <- 0.1*sv[["t50"]]
    if ("L1" %in% svnms) tmp["L1"] <- 0.5*sv[["L1"]]
    if ("L2" %in% svnms) tmp["L2"] <- 0.5*(sv[["L1"]]+sv[["L2"]])
    if ("L3" %in% svnms) tmp["L3"] <- ifelse("L2" %in% svnms,0.5*(sv[["L2"]]+sv[["L3"]]),0.5*(sv[["L1"]]+sv[["L3"]]))
    # reduce to only those in sv
    tmp <- tmp[which(names(tmp) %in% svnms)]
    # make sure they are in the same order as in sv
    tmp[svnms]
  }  # end iMake.slMins
  
  ## internal function to make maximum values for the sliders
  iMake.slMaxs <- function(sv,age) {
    svnms <- names(sv)
    tmp <- c(Linf=NA,K=NA,t0=NA,L0=NA,omega=NA,t50=NA,L1=NA,L2=NA,L3=NA)
    if ("Linf" %in% svnms) tmp["Linf"] <- 1.5*sv[["Linf"]]
    if ("K" %in% svnms) tmp["K"] <- 2*sv[["K"]]
    if ("t0" %in% svnms) tmp["t0"] <- 5
    if ("L0" %in% svnms) tmp["L0"] <- 2*sv[["L0"]]
    if ("omega" %in% svnms) tmp["omega"] <- 2*sv[["omega"]]
    if ("t50" %in% svnms) tmp["t50"] <- 0.6*max(age,na.rm=TRUE)
    if ("L1" %in% svnms) tmp["L1"] <- ifelse("L2" %in% svnms,0.5*(sv[["L1"]]+sv[["L2"]]),
                                             0.5*(sv[["L1"]]+sv[["L3"]]))
    if ("L2" %in% svnms) tmp["L2"] <- 0.5*(sv[["L2"]]+sv[["L3"]])
    if ("L3" %in% svnms) tmp["L3"] <- 1.5*sv[["L3"]]
    # reduce to only those in sv
    tmp <- tmp[which(names(tmp) %in% svnms)]
    # make sure they are in the same order as in sv
    tmp[svnms]
  } # end iMake.slMaxs
  
  ## internal function to make delta values for the sliders
  iMake.slDeltas <- function(sv) {
    svnms <- names(sv)
    tmp <- c(Linf=NA,K=NA,t0=NA,L0=NA,omega=NA,t50=NA,L1=NA,L2=NA,L3=NA)
    if ("Linf" %in% svnms) tmp["Linf"] <- 0.01*sv[["Linf"]]
    if ("K" %in% svnms) tmp["K"] <- 0.01
    if ("t0" %in% svnms) tmp["t0"] <- 0.01
    if ("L0" %in% svnms) tmp["L0"] <- 0.01*sv[["L0"]]
    if ("omega" %in% svnms) tmp["omega"] <- 0.01*sv[["omega"]]
    if ("t50" %in% svnms) tmp["t50"] <- 0.01
    if ("L1" %in% svnms) tmp["L1"] <- 0.01*sv[["L1"]]
    if ("L2" %in% svnms) tmp["L2"] <- 0.01*sv[["L2"]]
    if ("L3" %in% svnms) tmp["L3"] <- 0.01*sv[["L3"]]
    # reduce to only those in sv
    tmp <- tmp[which(names(tmp) %in% svnms)]
    # make sure they are in the same order as in sv
    tmp[svnms]
  } # end iMake.slDeltas
  
  ## Main function
  ## The default values for the sliders will be at the starting
  ## values as determined above.  Unlist first to make as a vector.
  sl.defaults <- unlist(sv)
  ## Grab names from the sv vector
  sl.names <- names(sl.defaults)
  ## Make minimum, maximum and delta values
  sl.mins <- iMake.slMins(sl.defaults)
  sl.maxs <- iMake.slMaxs(sl.defaults,age)
  sl.deltas <- iMake.slDeltas(sl.defaults)
  ## Make a title
  sl.ttl <- paste0("Von Bertalanffy (",type,")")
  ## Set up names that are specific to type and param
  relax::gslider(refresh,prompt=TRUE,hscale=1.5,pos.of.panel="left",
                 title=sl.ttl,sl.names=sl.names,
                 sl.mins=sl.mins,sl.maxs=sl.maxs,
                 sl.deltas=sl.deltas,sl.defaults=sl.defaults)
}

#=============================================================
# Constructs the actual plot in the dynamics plots for finding
# starting values
#=============================================================
iVBDynPlot <- function(age,len,type,p1,p2,p3,ages2use) {
  ## create a sequence of age values
  max.age <- max(age,na.rm=TRUE)
  x <- seq(0,max.age,length.out=20*max.age)
  ## create a sequence of matching recruit values according to
  ## the model type and param
  y <- iVBDynPlot_makeL(x,type,p1,p2,p3,ages2use)
  ## Construct the scatterplot with superimposed model
  opar <- graphics::par(mar=c(3.5,3.5,1.25,1.25), mgp=c(2,0.4,0), tcl=-0.2, pch=19)
  graphics::plot(age,len,xlab="Age",ylab="Mean Length")
  graphics::lines(x,y,lwd=2,lty=1,col="blue")
  graphics::par(opar)
}

#=============================================================
# Construct values for length given values of age, a Von B
# model and values of the parameters.  This is called by iVBDynPlot()
#=============================================================
iVBDynPlot_makeL <- function(x,type,p1,p2,p3,ages2use){
  switch(type,
         Typical=            { # p1=Linf, p2=K,  p3=to
                               y <- p1*(1-exp(-p2*(x-p3))) },
         Original=           { # p1=Linf, p2=L0, p3=K
                               y <- (p1-(p1-p2)*exp(-p3*x)) },
         GQ=, GallucciQuinn= { # p1=omega,p2=K,  p3=t0
                               y <- (p1/p2)*(1-exp(-p2*(x-p3))) },
         Weisberg=           { # p1=Linf, p2=t50,  p3=to
                               y <- p1*(1-exp(-(log(2)/(p2-p3))*(x-p3))) },
         Mooij=              { # p1=Linf, p2=L0, p3=omega
                               y <- p1-(p1-p2)*exp(-(p3/p1)*x) },
         Francis=            { # p1=L1, p2=L2, p3=L3
                               r <- (p3-p2)/(p2-p1)
                               t <- c(ages2use[1],mean(ages2use),ages2use[2])
                               y <- p1+(p3-p1)*((1-r^(2*((x-t[1])/(t[3]-t[1]))))/(1-r^2)) },
         Schnute=            { # p1=L1, p2=L3, p3=K
                               t <- ages2use
                               y <- p1+(p2-p1)*((1-exp(-p3*(x-t[1])))/(1-exp(-p3*(t[2]-t[1])))) }
  ) # end type switch
  y
}