#' @title Dynamics plots to explore typical fisheries growth models.
#'
#' @description Plots hypothetical size-at-age for one of seven possible parameterizations of the von Bertalanffy, three possible parameterizations of the Gompertz, and the Schnute growth models.  Slider bars are used to dynamically alter the parameters of each model.
#'
#' @details This function can be used to explore the \dQuote{shape} of the growth models for various choices of the parameters.  In this usage the \code{x} and \code{y} arguments should be (left) set at \code{NULL}.  This function can also be used to visually \dQuote{fit} a growth model to a set of observed lengths and ages.  This usage may be used to provide reasonable starting values for the parameters when fitting the growth model to the data with non-linear least-squares.  The observed data are plotted by including a formula of the form \code{length~age} in \code{x} and a data frame from which to draw the variables in the formula in the \code{data} arguments.
#'
#' The \code{type} argument is used to choose a type of growth model and must be one of the following (the models can be seen with \code{\link[FSA]{vbModels}}, \code{\link[FSA]{GompertzModels}}, and \code{\link[FSA]{SchnuteModels}}):
#'
#' \tabular{ll}{
#' \code{"vbTypical"} \tab The "typical" Beverton-Holt parameterized von Bertalanffy model.\cr
#' \code{"vbOriginal"} \tab The original parameterization from von Bertalanffy.\cr
#' \code{"vbMooij"} \tab The Mooij et al (1999) paramaterization of the von Bertalanffy model.\cr
#' \code{"vbGQ"} \tab The Gallucci & Quinn (1979) parameterization of the von Bertalanffy model.\cr
#' \code{"vbWeisberg"} \tab The Weisberg et al. (2010) parameterization of the von Bertalanffy model.\cr
#' \code{"vbSchnute"} \tab The Schnute-like paramaterization of the von Bertalanffy model.\cr 
#' \code{"vbTypicalW"} \tab The "typical" Beverton-Holt parameterized von Bertalanffy model, but for weights rather than lengths (thus, includes one more parameter).\cr
#' \code{"vbOriginalW"} \tab The original parameterization from von Bertalanffy, but for weights rather than lengths (thus, includes one more parameter).\cr
#' \code{"Gompertz1"} \tab The "first" parameterization of the Gompertz model.\cr 
#' \code{"Gompertz2"} \tab The "second" parameterization of the Gompertz model.\cr 
#' \code{"Gompertz3"} \tab The "third" parameterization of the Gompertz model.\cr
#' \code{"schnute"} \tab The Schnute(1981) four-parameter general growth model.
#' }
#'
#' @param type A single character string that indicates which growth model to use.  See details.
#' @param max.len A single numeric that indicates the maximum length to use in the simulations.
#' @param max.wt A single numeric that indicates the maximum weight to use in the simulations (only used of \code{type=} \code{"vbTypicalW"} or \code{"vbOriginalW"}.
#'
#' @return None.  However a dynamic graphic connected to slider bar controls in which the user can change the maximum age over which the growth model is evaluated and change the parameters specific to the chosen growth model.
#'
#' @author Derek H. Ogle, \email{dogle@@northland.edu}
#'
#' @seealso See model equations with \code{\link[FSA]{vbModels}}, \code{\link[FSA]{GompertzModels}}, and \code{\link[FSA]{SchnuteModels}}.  See similar functionality for von Bertalanffy models in \code{\link[FSA]{vbStarts}}.
#'
#' @references Francis, R.I.C.C.  1988.  Are growth parameters estimated from tagging and age-length data comparable?  Canadian Journal of Fisheries and Aquatic Sciences, 45:936-942.
#'
#' Gallucci, V.F. and T.J. Quinn II. 1979.  Reparameterizing, fitting, and testing a simple growth model.  Transactions of the American Fisheries Society, 108:14-25.
#'
#' Mooij, W.M., J.M. Van Rooij, and S. Wijnhoven.  1999.  Analysis and comparison of fish growth from small samples of length-at-age data: Detection of sequal dimorphism in Eurasian perch as an example.  Transactions of the American Fisheries Society 128:483-490.
#'
#' Schnute, J.  1981.  A versatile growth model with statistically stable parameters. Canadian Journal of Fisheries & Aquatic Sciences, 38:1128-1140.
#'
#' Schnute, J. and D. Fournier. 1980.  A new approach to length-frequency analysis: Growth structure.  Canadian Journal of Fisheries and Aquatic Sciences, 37:1337-1351.
#'
#' Weisberg, S., G.R. Spangler, and L. S. Richmond. 2010. Mixed effects models for fish growth. Canadian Journal of Fisheries And Aquatic Sciences 67:269-277.
#'
#' @keywords iplot
#'
#' @examples
#' ## ONLY RUN IN INTERACTIVE MODE
#' if (interactive()) {
#'
#' # Explore growth models (no data) -- use the defaults
#' growthModelSim()
#'
#' # Schnute parameterization of the von Bertalanffy model
#' growthModelSim(type="vbSchnute")
#'
#' ## Explore growth models superimposed on length-at-age data
#' # get Smallmouth Bass data from FSA package
#' data(SMBassWB)
#'
#' # interactively "fit" the typical paramaterization of the von Bertalanffy model to the data
#' growthModelSim(lencap~agecap,data=SMBassWB)
#'
#' # interactively "fit" the second paramaterization of the Gompertz model to the data
#' growthModelSim(lencap~agecap,data=SMBassWB,type="Gompertz2")
#'
#' } ## END IF INTERACTIVE MODE
#'
#' @export growthModelSim
#' 
## Main function
growthModelSim <- function(type=c("vbTypical","vbOriginal","vbGQ","vbGallucciQuinn","vbMooij",
                                  "vbWeisberg","vbSchnute","vbTypicalW","vbOriginalW",
                                  "Gompertz1","Gompertz2","Gompertz3",
                                  "Schnute"),
                           max.len=500,max.wt=500) {
  # Trying to deal with no visible bindings problem
  p1 <- p2 <- p3 <- p4 <- t.max <- NULL
  type <- match.arg(type)
  if (!iCheckRStudio()) stop("'growthModelSim' only works in RStudio.",call.=FALSE)
  if (iChk4Namespace("manipulate")) {
    switch(type,
           vbOriginal= {
             manipulate::manipulate(
               {
                 iGrowthSimPlot(type,t.max,p1,p2,p3,NULL)
               },
               p1=manipulate::slider(100,800,step=10,initial=500,label="L_inf"),
               p2=manipulate::slider(0,90,step=10,initial=0,label="L_0"),
               p3=manipulate::slider(0,0.5,step=0.05,initial=0.3,label="K"),
               t.max=manipulate::slider(5,60,step=1,initial=10,label="Max Age")
             )},
           vbTypical= {
             manipulate::manipulate(
               {
                 iGrowthSimPlot(type,t.max,p1,p2,p3,NULL)
               },
               p1=manipulate::slider(100,800,step=10,initial=500,label="L_inf"),
               p2=manipulate::slider(0,0.5,step=0.05,initial=0.3,label="K"),
               p3=manipulate::slider(-4,4,step=0.5,initial=0,label="t_0"),
               t.max=manipulate::slider(5,60,step=1,initial=10,label="Max Age")
             )},
           vbOriginalW= {
             manipulate::manipulate(
               {
                 iGrowthSimPlot(type,t.max,p1,p2,p3,p4)
               },
               p1=manipulate::slider(100,1000,step=10,initial=500,label="W_inf"),
               p2=manipulate::slider(0,90,step=10,initial=0,label="W_0"),
               p3=manipulate::slider(0,0.5,step=0.05,initial=0.3,label="K"),
               p4=manipulate::slider(0.25,4,step=0.05,initial=3,label="b"),
               t.max=manipulate::slider(5,60,step=1,initial=10,label="Max Age")
             )},
           vbTypicalW= {
             manipulate::manipulate(
               {
                 iGrowthSimPlot(type,t.max,p1,p2,p3,p4)
               },
               p1=manipulate::slider(100,800,step=10,initial=500,label="W_inf"),
               p2=manipulate::slider(0,0.5,step=0.05,initial=0.3,label="K"),
               p3=manipulate::slider(-4,4,step=0.5,initial=0,label="t_0"),
               p4=manipulate::slider(0.25,4,step=0.05,initial=3,label="b"),
               t.max=manipulate::slider(5,60,step=1,initial=10,label="Max Age")
             )},
           vbGQ=, vbGallucciQuinn= {
             manipulate::manipulate(
               {
                 iGrowthSimPlot(type,t.max,p1,p2,p3,NULL)
               },
               p1=manipulate::slider(0,200,step=5,initial=75,label="omega"),
               p2=manipulate::slider(0,0.5,step=0.05,initial=0.3,label="K"),
               p3=manipulate::slider(-4,4,step=0.5,initial=0,label="t_0"),
               t.max=manipulate::slider(5,60,step=1,initial=10,label="Max Age")
             )},
           vbMooij= {
             manipulate::manipulate(
               {
                 iGrowthSimPlot(type,t.max,p1,p2,p3,NULL)
               },
               p1=manipulate::slider(100,800,step=10,initial=500,label="L_inf"),
               p2=manipulate::slider(0,90,step=10,initial=0,label="L_0"),
               p3=manipulate::slider(0,200,step=5,initial=75,label="omega"),
               t.max=manipulate::slider(5,60,step=1,initial=10,label="Max Age")
             )},
           vbWeisberg= {
             manipulate::manipulate(
               {
                 iGrowthSimPlot(type,t.max,p1,p2,p3,NULL)
               },
               p1=manipulate::slider(100,800,step=10,initial=500,label="L_inf"),
               p2=manipulate::slider(1,50,step=0.1,initial=5,label="t_50"),
               p3=manipulate::slider(-4,4,step=0.5,initial=0,label="t_0"),
               t.max=manipulate::slider(5,60,step=1,initial=10,label="Max Age")
             )},
           vbSchnute= {
             manipulate::manipulate(
               {
                 iGrowthSimPlot(type,t.max,p1,p2,p3,p4)
               },
               p1=manipulate::slider(0,90,step=10,initial=0,label="L_1"),
               p2=manipulate::slider(100,800,step=10,initial=500,label="L_2"),
               p3=manipulate::slider(0,0.5,step=0.05,initial=0.3,label="K"),
               p4=manipulate::slider(0.1,3,step=0.1,initial=1,label="b"),
               t.max=manipulate::slider(5,60,step=1,initial=10,label="Max Age")
             )},
           Gompertz1= {
             manipulate::manipulate(
               {
                 iGrowthSimPlot(type,t.max,p1,p2,p3,NULL)
               },
               p1=manipulate::slider(0,90,step=10,initial=0,label="L_0"),
               p2=manipulate::slider(0,2,step=0.01,initial=1.5,label="G"),
               p3=manipulate::slider(0,2,step=0.01,initial=0.5,label="g"),
               t.max=manipulate::slider(5,60,step=1,initial=10,label="Max Age")
             )},
           Gompertz2= {
             manipulate::manipulate(
               {
                 iGrowthSimPlot(type,t.max,p1,p2,p3,NULL)
               },
               p1=manipulate::slider(100,800,step=10,initial=500,label="L_inf"),
               p2=manipulate::slider(0,2,step=0.01,initial=0.5,label="g"),
               p3=manipulate::slider(0,10,step=0.05,initial=2,label="t*"),
               t.max=manipulate::slider(5,60,step=1,initial=10,label="Max Age")
             )},
           Gompertz3= {
             manipulate::manipulate(
               {
                 iGrowthSimPlot(type,t.max,p1,p2,p3,NULL)
               },
               p1=manipulate::slider(100,800,step=10,initial=500,label="L_inf"),
               p2=manipulate::slider(0,2,step=0.01,initial=0.5,label="g"),
               p3=manipulate::slider(-4,4,step=0.5,initial=0,label="t_0"),
               t.max=manipulate::slider(5,60,step=1,initial=10,label="Max Age")
             )},
           Schnute= {
             manipulate::manipulate(
               {
                 iGrowthSimPlot(type,t.max,p1,p2,p3,p4)
               },
               p1=manipulate::slider(0,90,step=10,initial=0,label="L_1"),
               p2=manipulate::slider(100,800,step=10,initial=500,label="L_2"),
               p3=manipulate::slider(-2,1,step=0.1,initial=0.1,label="c"),
               p4=manipulate::slider(0.1,2,step=0.1,initial=1,label="d"),
               t.max=manipulate::slider(5,60,step=1,initial=10,label="Max Age")
             )}
    ) # end switch
  }
}


## Internal function to predict length/weight for plotting -- used in iGrowthSimPlot
iPredLength <- function(type,t,p1,p2,p3,p4) {
  switch(type,
         # p1=Linf, p2=K,  p3=to, p4 not used
         vbTypical= {  sd <- p1*(1-exp(-p2*(t-p3))) },
         # p1=Winf, p2=K,  p3=to, p4=b
         vbTypicalW= { sd <- p1*(1-exp(-p2*(t-p3)))^p4 },         
         # p1=Linf, p2=L0, p3=K,  p4 not used
         vbOriginal= { sd <- (p1-(p1-p2)*exp(-p3*t)) },
         # p1=Winf, p2=W0, p3=K,  p4=b
         vbOriginalW= { sd <- (p1-(p1-p2)*exp(-p3*t))^p4 },
         # p1=omega,p2=K,  p3=t0, p4 not used
         vbGQ=, vbGallucciQuinn= { sd <- (p1/p2)*(1-exp(-p2*(t-p3))) },
         # p1=Linf, p2=t50,  p3=to, p4 not used
         vbWeisberg= {  sd <- p1*(1-exp(-(log(2)/(p2-p3))*(t-p3))) },
         # p1=Linf, p2=L0, p3=ome,p4 not used
         vbMooij= { sd <- p1-(p1-p2)*exp(-(p3/p1)*t) },
         # p1=L1,   p2=L2, p3=K,  p4=b
         vbSchnute= {
           sd <- ((p1^p4)+((p2^p4)-(p1^p4))*((1-exp(-p3*(t-min(t))))/(1-exp(-p3*(max(t)-min(t))))))^(1/p4)
         },
         # p1=Lo,   p2=G,  p3=g,  p4 not used
         Gompertz1= { sd <- p1*(exp(p2*(1-exp(-p3*t)))) },
         # p1=Linf, p2=g,  p3=t*, p4 not used
         Gompertz2= { sd <- p1*(exp(-exp(-p2*(t-p3)))) },
         # p1=Linf, p2=g,  p3=t0, p4 not used
         Gompertz3= { sd <- p1*(exp((-1/p2)*exp(-p2*(t-p3)))) },
         # p1=L1,   p2=L2, p3=c,  p4=d
         Schnute= {
           minage <- min(t)
           maxage <- max(t)
           diffage <- maxage-minage
           if (p3==0) {
             if (p4==0) {
               # Case 4
               sd <- p1*exp(log(p2/p1)*(t-minage)/diffage)
             } else { 
               # Case 3
               sd <- (p1^p4+(p2^p4-p1^p4)*(t-minage)/diffage)^(1/p4)
             }
           } else {
             if (p4==0) {
               # Case 2
               sd <- p1*exp(log(p2/p1)*(1-exp(-p3*(t-minage)))/(1-exp(-p3*diffage)))
             } else { 
               # Case 1
               sd <- (p1^p4+(p2^p4-p1^p4)*((1-exp(-p3*(t-minage)))/
                                             (1-exp(-p3*diffage))))^(1/p4) 
             }
           }
         } # end Schnute
  ) # end switch 
  sd
} ## end iPredLength internal function

## internal function for constructing the plot
iGrowthSimPlot <- function(type,t.max,p1,p2,p3,p4) {
  t <- seq(0,t.max,length.out=t.max*20)
  vals <- iPredLength(type,t,p1,p2,p3,p4)
  ylbl <- ifelse (type %in% c("vbTypicalW","vbOriginalW"),"Weight","Length")
  opar <- graphics::par(mar=c(3.5,3.5,1,1),mgp=c(2,0.4,0),tcl=-0.2,xaxs="i",yaxs="i")
  graphics::plot(t,vals,type="l",lwd=2,col="blue",
                 xlab="Age",xlim=c(0,t.max),
                 ylab=ylbl,ylim=c(0,800))
  graphics::par(opar)
} ## end iGrowthSimPlot internal function
