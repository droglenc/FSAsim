#' @title Dynamic plots to explore typical fisheries stock-recruitment models.
#'
#' @description Plots modeled number of recruits versus stock size for common parameterizations of the Beverton-Holt and Rickerstock-recruit models. Slider bars are used to dynamically alter the parameters of each model.
#'
#' @details This function can be used to explore the dynamics of stock-recruitment models for various parameter choices. In these instances of model exploration the \code{S=} and \code{R=} arguments should be (left) set at \code{NULL}.
#' 
#' The \code{type=} argument is used to choose either the \code{"BevertonHolt"} or \code{"Ricker"} stock-recruitment models. Common parameterizations of the \code{"BevertonHolt"} and \code{"Ricker"} models can be chosen with \code{param=}. Four paramaterizations of the Beverton-Holt model and three parameterizations of the Ricker model are allowed. See \code{srFunShow} described in \code{\link[FSA]{stockRecruitment}} to see equations for each model.
#'
#' @param type A string that indicates the type of the stock-recruitment model. Must be one of \code{"BevertonHolt"}, \code{"Ricker"}, \code{"Shepherd"}, or \code{"SailaLorda"}.
#' @param param A numeric that indicates the parameterization of the stock-recruitment model type.
#' @param max.S A single numeric that indicates the maximum spawning stock to use for scaling the x-axis. Ignored if \code{S} is not NULL.
#' @param max.R A single numeric that indicates the maximum recruitment to use for scaling the y-axis. Ignored if \code{S} is not NULL.
#'
#' @return None. However a dynamic graphic connected to slider bar controls of the \bold{\sQuote{a}}, \bold{\sQuote{b}}, or \bold{\sQuote{Rp}} parameters specific to the chosen stock-recruit model.
#'
#' @author Derek H. Ogle, \email{dogle@@northland.edu}
#'
#' @seealso \code{srStartsDP} and \code{srFunShow} described in \code{\link[FSA]{stockRecruitment}}.
#'
#' @aliases srSim
#'
#' @keywords iplot
#'
#' @examples
#' ## ONLY RUN IN INTERACTIVE MODE
#' if (interactive()) {
#' # Beverton-Holt models
#' srSim()
#' srSim(param=2)
#' srSim(param=3)
#' srSim(param=4)
#'
#' # Ricker models
#' srSim(type="Ricker")
#' srSim(type="Ricker",param=2)
#' srSim(type="Ricker",param=3)
#' }  ## END INTERACTIVE MODE
#'
#' @rdname srSim
#' @export
srSim <- function(type=c("BevertonHolt","Ricker"),param=1,
                  max.S=500,max.R=1000) {
  p1 <- p2 <- NULL
  type <- match.arg(type)
  if (!iCheckRStudio()) FSA:::STOP("'srSim' only works in RStudio.")
  if (iChk4Namespace("manipulate")) {
    switch(type,
           BevertonHolt= {
             if (!param %in% 1:4)
               FSA:::STOP("'param' must be in 1:4 when type='BevertonHolt'.")
             if (param==1) {
               delta.a <- max.R/10000  # works OK for default max.S and max.R
               delta.b <- max.S/50000
               max.rds <- 500*delta.a
               manipulate::manipulate(
                 { isrSimPlot(type,param,max.S,max.R,max.rds,p1,p2,NULL) },
                 p1=manipulate::slider(10*delta.a,500*delta.a,step=delta.a,
                                       initial=500*delta.a,label="a"),
                 p2=manipulate::slider(5*delta.b,50*delta.b,step=delta.b,
                                       initial=5*delta.b,label="b"))
             } else if (param==2) {
               delta.a <- max.R/10000  # works OK for default max.S and max.R
               max.rds <- 1000*delta.a
               manipulate::manipulate(
                 { isrSimPlot(type,param,max.S,max.R,max.rds,p1,p2,NULL) },
                 p1=manipulate::slider(10*delta.a,1000*delta.a,step=delta.a,
                                       initial=500*delta.a,label="a"),
                 p2=manipulate::slider(0.02*max.R,max.R,step=0.02*max.R,
                                       initial=max.R,label="Rp"))
             } else if (param==3) {
               delta.a <- 1/max.R  # works OK for default max.S and max.R
               delta.b <- 0.0001*max.S/max.R
               max.rds <- 1/(10*delta.a)
               manipulate::manipulate(
                 { isrSimPlot(type,param,max.S,max.R,max.rds,p1,p2,NULL) },
                 p1=manipulate::slider(10*delta.a,200*delta.a,step=delta.a,
                                       initial=10*delta.a,label="a"),
                 p2=manipulate::slider(20*delta.b,120*delta.b,step=delta.b,
                                       initial=20*delta.b,label="b"))
             } else {
               delta.a <- 1/max.R  # works OK for default max.S and max.R
               max.rds <- 1/(10*delta.a)
               manipulate::manipulate(
                 { isrSimPlot(type,param,max.S,max.R,max.rds,p1,p2,NULL) },
                 p1=manipulate::slider(10*delta.a,200*delta.a,step=delta.a,
                                       initial=10*delta.a,label="a"),
                 p2=manipulate::slider(0.02*max.R,max.R,step=0.02*max.R,
                                       initial=max.R,label="Rp"))
             }
           },
           Ricker= {
             if (!param %in% 1:3) 
               FSA:::STOP("'param' must be in 1:3 when type='Ricker'.")
             if (param==1) {
               delta.a <- 10/max.R  # works OK for default max.S and max.R
               delta.b <- 0.002*max.S/max.R 
               max.rds <- 2000*delta.a
               manipulate::manipulate(
                 { isrSimPlot(type,param,max.S,max.R,max.rds,p1,p2,NULL) },
                 p1=manipulate::slider(40*delta.a,2000*delta.a,step=delta.a,
                                       initial=2000*delta.a,label="a"),
                 p2=manipulate::slider(8*delta.b,100*delta.b,step=delta.b,
                                       initial=8*delta.b,label="b"))
             } else if (param==2) {
               delta.a <- 10/max.R  # works OK for default max.S and max.R
               delta.b <- 0.002*max.S/max.R    
               max.rds <- exp(300*delta.a)
               manipulate::manipulate(
                 { isrSimPlot(type,param,max.S,max.R,max.rds,p1,p2,NULL) },
                 p1=manipulate::slider(delta.a,300*delta.a,step=delta.a,
                                       initial=300*delta.a,label="a"),
                 p2=manipulate::slider(8*delta.b,100*delta.b,step=delta.b,
                                       initial=8*delta.b,label="b"))
             } else {
               delta.a <- 10/max.R
               max.rds <- 6000*delta.a
               manipulate::manipulate(
                 { isrSimPlot(type,param,max.S,max.R,max.rds,p1,p2,NULL) },
                 p1=manipulate::slider(550*delta.a,6000*delta.a,step=delta.a,
                                       initial=6000*delta.a,label="a"),
                 p2=manipulate::slider(0.02*max.R,max.R,step=0.02*max.R,
                                       initial=max.R,label="Rp"))
             } 
          })
  }
}

## Internal plotting function
isrSimPlot <- function(type,param,max.S,max.R,max.rds,a,b,c) {
  # create a sequence of stock values
  x <- c(seq(0,0.2*max.S,length.out=500),seq(0.2*max.S,max.S,length.out=500)) 
  # create a sequence of recruit values based on model choice
  if (type=="BevertonHolt") {  
    switch(toString(param),
           "1"={Rp <- a/b               # calculated Rp
              y <- a*x/(1+b*x) },       
           "2"={Rp <- b                 # 2nd param is Rp in this parameteriz.
              y <- a*x/(1+a*(x/Rp)) },
           "3"={Rp <- 1/b               # calculated Rp
              y <- x/(a+b*x) },
           "4"={Rp <- b                 # 2nd param is Rp in this parameteriz.
              y <- x/(a+x/Rp)}
           ) # end B-H switch
  } else {
    switch(toString(param),
           "1"={Rp <- a/(b*exp(1))      # calculated Rp
              Sp <- 1/b
              y <- a*x*exp(-b*x)  },
           "2"={Rp <- exp(a)/(b*exp(1)) # calculated Rp
              Sp <- 1/b
              y <- x*exp(a-b*x) },
           "3"={Rp <- b                 # 2nd param is Rp in this parameteriz.
              Sp <- (Rp*exp(1))/a
              y <- a*x*exp(-a*x/(Rp*exp(1))) }
           ) # end Ricker switch
  }
  withr::local_par(mfrow=c(1,2),mar=c(3.5,3.5,1.25,1.25),
                   mgp=c(2,0.4,0),tcl=-0.2, pch=19)
  graphics::plot(x,y,xlab="Parental (Spawner) Stock",ylab="Recruits",
                 type="l",lwd=2,col="blue",ylim=c(0,max.R),xlim=c(0,max.S))
  graphics::abline(h=Rp,lwd=2,lty=3,col="red")
  graphics::axis(4,Rp,"Rp",col.ticks="red",col.axis="red",las=1)
  if (type=="Ricker") {
    graphics::abline(v=Sp,lwd=2,lty=3,col="red")
    graphics::axis(3,Sp,"Sp",col.ticks="red",col.axis="red")
  }
  graphics::plot(x,y/x,xlab="Parental (Spawner) Stock",ylab="Recruits/Spawner",
                 type="l",lwd=2,col="blue",ylim=c(0,max.rds),xlim=c(0,max.S))
} # end isrSimPlot internal function
