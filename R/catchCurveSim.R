#' @title Dynamic plots to explore catch-curve models.
#'
#' @description Constucts hypothetical catch-at-age data given choices for a hypothetical population (governed by the initial population size (No) and instantaneous mortality rate (Z)) and fishery (through age at recruitment to the gear and proportional vulnerabilities to the gear).  The natual log of catch versus age is plotted (i.e., a catch curve) and slider bars can be used to change population or fisheries parameter values.  This function an corresponding dynamics graphic can be used to explore the effects of a parameter on the catch curve and, thus, elucidate the effects of assumption violations.
#'
#' @details A catch curve plot that is dynamically linked to a set of slider controls.  The  user can then alter the population or fishery parameters to see the effect of those changes on the catch curve.  Each parameter that can be set or can be altered with the slider bars is described below.
#'
#' The \bold{\sQuote{Z mean}} controls the mean instantaneous mortality rate to use in the simulations.  If \bold{\sQuote{Z cv}} (see below) is set to 0 then this is the exact mortality rate used in the simulations.
#'
#' The \bold{\sQuote{Z cv}} controls the coefficient of variation (as a proportion, not a percentage) to use when modeling the instantaneous mortality rate.  If \bold{\sQuote{Z cv}} is set to 0 then no randomization is used.  If \bold{\sQuote{Z cv}} is greater than 0 then a random normal deviate centered on \bold{\sQuote{Z mean}} with a standard deviation of \bold{\sQuote{Z cv}}*\bold{\sQuote{Z mean}} is used to model the instantaneous mortality rate.
#'
#' A simulated change in instantaneous mortality rate is controlled with the \bold{\sQuote{Z delta}} slider bar and the \code{ZdeltaAge}, and \code{ZSteady} arguments.  The \bold{\sQuote{Z.delta}} slider controls a proportional adjustment in the mortality rate that begins at \code{ZdeltaAge}.  The adjustment is a constant proportion through time/ages if \code{ZSteady=TRUE}, but grows geometrically if \code{ZSteady=FALSE}.  For example, if \bold{\sQuote{Z mean}} is set to 0.40, \bold{\sQuote{Z cv}} is set to 0 (i.e., no variability), \bold{\sQuote{Z delta}} is set to 1.2, \code{ZdeltaAge=5}, and \code{ZSteady=TRUE} then the vector of mortalities used for ages 1 to 10 would be
#'
#' \code{c(0.4,0.4,0.4,0.4,0.4*1.2,0.4*1.2,0.4*1.2,0.4*1.2,0.4*1.2,0.4*1.2)}
#'
#' Alternatively, if all values were the same except that \code{ZSteady=FALSE}, then the vector of mortalities used for ages 1 to 10 would be 
#' 
#' \code{c(0.4,0.4,0.4,0.4,0.4*1.2,0.4*1.2^2,0.4*1.2^3,0.4*1.2^4,0.4*1.2^5,0.4*1.2^6)}
#'
#' The simulation code can consider two different tacks when modeling variability in the instantaneous mortality rate. The first tack considers that the mortality rate experienced throughout the existence of each year-class corresponding to ages in the cross-sectional sample older than some defined age was different than all previous year-classes.  Simulations using this tack can be used to assess the assumption of a constant mortality rate for all year-classes.  The second tack considers that the mortality rate for all ages prior to some defined age was the same for all year-classes but that a different mortality rate for all ages after some defined age for all year-classes should be used.  Simulations using this tack can be used to assess the assumption of a constant mortality rate for all ages.  Which track is used in a set of simulations is set by the \code{Zycl=} argument to the function (and not by a slider bar).  If \code{Zcyl=TRUE} then the first tack is used.
#'
#' Slider bars are used to control similar aspects of the initial population size (No).  \bold{\sQuote{No mean}} and \bold{\sQuote{No CV}} are the mean and coefficient of variation to be used in the simulations.  Changes in the initial population size are simulated with the \bold{\sQuote{No delta}} slider and the \code{NodeltaAge} and \code{NoSteady} arguments.  These adjustments are applied in the same manner as described for the instantaneous mortality rate.
#'
#' The vulnerabilties of each age of fish to the fishing gear can be simulated in one of two ways.  If a vector of \dQuote{1}s is sent in the \code{v} argument (this is the default) then the fish in all ages after the recruitment age (which is selected with a slider bar) are fully and equally vulnerable to the gear.  Fish in age-classes prior to the age of recruitment are modeled by a linear increase in proportional vulnerability from a value of 0.1 to 1 (for the age at recruitment).  A wider variety of changes in vulnerabilities can be modeled by sending various vectors of proportional vulnerabilities for each age in the simulation.  For example, if 10 ages are modeled, then the vector \code{c(0.1,0.3,0.6,1,1,1,0.8,0.6,0.4,0.2)} could be used to model incomplete vulnerabilities prior to age 4, complete vulnerablities from age 4 to age 6, and then incomplete vulnerabilities for age 7 and older.  If the vector of vulnerabilities sent does not contain all \dQuote{1}s then a recruitment age slider will not be available.
#'
#' @note The range of values allowed for each of the parameters was chosen to allow a wide variety of modeling opportunities.  However, it is highly likely that these ranges do not encompass every possible set of values that a user may wish to view.  Thus, this function is more useful as a learning tool rather than a research simulation tool.
#' 
#' @param max.age A single numeric integer that indicates the oldest age to use in the simulations.
#' @param v A single vector of length \code{max.age} of the proportional vulnerabilities to the gear.  See details.
#' @param Zycl A single logical that indicates whether changes in mortality rates are applied to and followed for each year-class or not.  See details.
#' @param ZdeltaAge A single numeric that indicates the age at which the mortality rate should change.  See details.
#' @param NodeltaAge A single numeric that indicates the age at which the initial population size should change.  See details.
#' @param ZSteady A single logical that indicates whether the change in mortality rate will be steady (\code{=TRUE}; default) or not.  See details.
#' @param NoSteady A single logical that indicates whether the change in initialpopulation size will be steady (\code{=TRUE}; default) or not.  See details.
#' @param recruit.age A single numeric that indicates what the age at full recruitment should be.  See detials.
#'
#' @return None.  A slide controlled dynamic graphic with the following description is produced.  The catch curve graphic is shown for as many as three lines.  The gray line is the catch-at-age data for the default values (the initial values when the simulation was begun).  This line is used simply as a basis for examining parameter changes.  The blue line is the catch-ag-age data for current choices of \bold{\sQuote{Z.mean}}, \bold{\sQuote{Z.cv}}, \bold{\sQuote{No.mean}}, and \bold{\sQuote{No.cv}}, but not for any of the other choices (i.e., no change in the \sQuote{delta} values).
#'
#' In other words, the blue line reflects the model for other than default parameter choices but WITHOUT any assumption violations.  This line serves as a basis for different parameter choices without assumption violations.  The red line is the catch-at-age data for all current choices of parameters in the slider box.  The lines are plotted in the order of \dQuote{gray}, \dQuote{red}, \dQuote{blue} so, if any two are equal then the color first plotted will not be seen.
#' 
#' @author Derek H. Ogle, \email{dogle@@northland.edu}
#' 
#' @seealso \code{\link{catchCurve}}
#' 
#' @keywords iplot
#' 
#' @examples
#' if (interactive()) {
#' ## default modeling
#' catchCurveSim()
#'
#' ## modeling with custom vulnerabilities
#' catchCurveSim(v=c(0.1,0.3,0.6,1,1,0.8,0.7,0.6,0.5,0.4))
#'
#' } # end if interactive
#' 
#'@export
#'
catchCurveSim <- function(max.age=15,v=rep(1,max.age),
                          Zycl=FALSE,ZdeltaAge=round(max.age/2,0),
                          NodeltaAge=ZdeltaAge,ZSteady=TRUE,NoSteady=ZSteady,
                          recruit.age=3) {
  type <- EM <- En <- NULL
  if (iChk4Namespace("manipulate")) {
    manipulate::manipulate(
      {
        set.seed(sample(1:100000))
        Z.param <- list(mean=Z.mean,cv=Z.cv,delta=Z.delta,
                        deltaAge=ZdeltaAge,dstdy=ZSteady)
        No.param <- list(mean=No.mean,cv=No.cv,delta=No.delta,
                         deltaAge=NodeltaAge,dstdy=NoSteady)    
        iCC.plot(max.age,v,Zycl,ZdeltaAge,NodeltaAge,recruit.age,Z.param,No.param)
      },
      Z.mean=manipulate::slider(0.05,1,step=0.05,initial=0.4,label="Z mean"),
      No.mean=manipulate::slider(100,2000,step=100,initial=1000,label="No mean"),
      Z.cv=manipulate::slider(0,0.2,step=0.05,initial=0,label="Z CV"),
      No.cv=manipulate::slider(0,0.5,step=0.05,initial=0,label="No CV"),
      Z.delta=manipulate::slider(0.3,2,step=0.1,initial=1,label="Z delta"),
      No.delta=manipulate::slider(0.3,2,step=0.1,initial=1,label="No delta"),
      steady=manipulate::checkbox(label="")
      rerand=manipulate::button("Rerandomize")
    )
  }
}  



## Internal function to create the plot
iCC.plot <- function(max.age,v,Zycl,ZdeltaAge,NodeltaAge,recruit.age,Z.param,No.param) {
  # create ages to model over
  Ages <- 1:max.age
  # Make catches for the current slider selections (include possible violations)
  C.ass <- iCC.makeCatch(Ages,v,Zycl,recruit.age,Z.param,No.param)
  # Make catches for the default settings
  Z.param.def <- list(mean=0.40,cv=0,delta=1,deltaAge=ZdeltaAge,dstdy=TRUE)
  No.param.def <- list(mean=1000,cv=0,delta=1,deltaAge=NodeltaAge,dstdy=TRUE)
  C.def <- iCC.makeCatch(Ages,v,Zycl,recruit.age,Z.param.def,No.param.def)
  # Make catches for the modified graph but wihout assumption violations
  # if no assumption violations, make same as actual so no extra line due to randomization
  if ((Z.param$delta==1) & (No.param$delta==1)) C.mod <- C.ass  
  else {
    Z.param.mod <- list(mean=Z.param$mean,cv=Z.param$cv,delta=1,
                        deltaAge=ZdeltaAge,dstdy=TRUE)
    No.param.mod <- list(mean=No.param$mean,cv=No.param$cv,delta=1,
                         deltaAge=NodeltaAge,dstdy=TRUE)    
    C.mod <- iCC.makeCatch(Ages,v,Zycl,recruit.age,Z.param.mod,No.param.mod)
  }
  # Make the plot
  op <- graphics::par(mar=c(3.5,3.5,1.5,1.5),mgp=c(2,0.4,0),tcl=-0.2,cex.main=0.75)
  graphics::plot(Ages,log(C.def),type="o",col="gray",pch=19,lwd=3,
                 ylab="log(Catch)",ylim=c(-3,8),main=iCC.title(Z.param,No.param))
  graphics::lines(Ages,log(C.ass),type="o",col="red",pch=19,lwd=3)
  graphics::lines(Ages,log(C.mod),type="o",col="blue",pch=19,lwd=3)
  graphics::legend("topright",c("default","assumptions not met","if assumptions met"),
                   lwd=2,pch=19,col=c("gray","red","blue"),cex=0.8)
  graphics::par(op)
} # end iCC.plot


## Internal function to make the catch vectors to plot
iCC.makeCatch <- function(Ages,v,Zycl,recruit.age,Z.param,No.param) {
  if (all(v==1)) { v <- iCC.makeV(Ages,recruit.age) }
  # Vector of initial year-class sizes
  No <- iCC.makeNoycl(Ages,No.param)
  if (Zycl) {
    # If Z will vary by year-class ...
    #   Vector of Zs for each year-class
    Z <- iCC.makeZycl(Ages,Z.param)
    #   Population sizes for each age
    N <- No*exp(-Z*Ages)
  } else {
    # or if Z varies by age (irrespective of year-class)
    #   Create Z vector with randomly select values from a Normal Distribution
    ifelse(Z.param$cv==0,Z <- rep(Z.param$mean,length(Ages)),
           Z <- rnorm(length(Ages),Z.param$mean,Z.param$mean*Z.param$cv))
    #   Adjust age vector for handling change in Zs by age
    Ages.adj <- iCC.makeZage(Ages,Z.param)
    #   Population sizes for each age
    N <- No*exp(-Z*Ages.adj)
  }
  v*N
} # end iCC.makeCatch


## Internal function to create a vector of vulnerabilities for each age
iCC.makeV <- function(Ages,recruit.age) {
  v <- rep(1,length(Ages))
  # slope for vulnerability for pre-recruitment ages
  slp <- .9/(recruit.age-1)
  for (i in 1:(recruit.age-1)) {
    # create linearly increasing vulnerabilities for pre-recruitment ages
    v[i] <- .1 + (i-1)*slp
  }
  v
} # end iCC.makeV

## Internal function to create a vector of initial year-class
##   sizes for each year-class
iCC.makeNoycl <- function(Ages,param) {
  #   Initiate No vector with randomly select values from a Normal Distribution
  ifelse(param$cv==0,No <- rep(param$mean,length(Ages)),
         No <- rnorm(length(Ages),param$mean,param$mean*param$cv))
  #   Initiate the multiplication vector with 1s
  md <- rep(1,length(Ages))
  if (param$delta!=1) {
    #   If a delta value is sent then must change the multiplication vector
    if (param$dstdy) {
      #   If a steady change, replace 1s with delta value
      md[param$deltaAge:length(Ages)] <- param$delta
    } else {
      #   If a non-steady change, ...
      for (i in (param$deltaAge:length(Ages))) {
        #   ... replace 1s with delta raised to power of number of ages past t1
        md[i] <- param$delta^(i-param$deltaAge+1)
      }
    }
  }
  #   Adjust initial No values by the multiplier vector
  No*md
} # end iCC.makeNoycl
  
## Internal function to create a vector of initial Zs for each
##   year-class (concept is same as make.No) 
iCC.makeZycl <- function(Ages,param) {
  #   Initiate Z vector with randomly select values from a Normal Distribution
  ifelse(param$cv==0,Z <- rep(param$mean,length(Ages)),
         Z <- rnorm(length(Ages),param$mean,param$mean*param$cv))
  #   Initiate the multiplication vector with 1s
  md <- rep(1,length(Ages))
  if (param$delta!=1) {
    #   If a delta value is sent then must change the multiplication vector
    if (param$dstdy) {
      #   If a steady change, ...
      #   ... replace 1s with delta value
      md[param$deltaAge:length(Ages)] <- param$delta
    } else {
      #   If a non-steady change, ...
      for (i in (param$deltaAge:length(Ages))) {
        #   ... replace 1s with delta raised to power of number of ages past t1
        md[i] <- param$delta^(i-param$deltaAge+1)
      }
    }
  }
  #   Adjust initial Z values by the multiplier vector
  Z*md
} # end iCC.makeZycl

## Internal function to create
iCC.makeZage <- function(Ages,param) {
  #   Initiate ages vector with the list of ages
  Ages.adj <- Ages
  if (param$delta!=1) {
    #   If a delta value is sent then must change the multiplication vector
    if (param$dstdy) {
      for (i in (param$deltaAge:length(Ages))) {
        Ages.adj[i] <- (1-param$delta)*Ages.adj[param$deltaAge]+param$delta*Ages.adj[i]
      }
    } else {
      for (i in (param$deltaAge:length(Ages))) {
        Ages.adj[i] <- Ages.adj[i-1]+param$delta^(i-param$deltaAge)
      }
    }
  }
  Ages.adj
} # end iCC.makeZage

## Internal function to create
iCC.title <- function(Z.param,No.param) {
  viol.Zdelta <- viol.Nodelta <- viol.rand <- FALSE
  lbl <- "VIOLATION:"
  if (Z.param$delta != 1) {
    viol.Zdelta <- TRUE
    if (Z.param$dstdy) {
      if (Z.param$delta < 1) { lbl <- paste(lbl,"Constant decreased mortality") }
      else { lbl <- paste(lbl,"Constant increased mortalilty") }        
    } else {
      if (Z.param$delta < 1) { lbl <- paste(lbl,"Constantly decreasing mortality") }
      else { lbl <- paste(lbl,"Constantly increasing mortalilty") }         
    }
  }
  if (No.param$delta != 1) {
    viol.Nodelta <- TRUE
    if (No.param$dstdy) {
      if (No.param$delta < 1) { lbl <- paste(lbl,"Constant decreased recruitment") }
      else { lbl <- paste(lbl,"Constant increased recruitment") }        
    } else {
      if (No.param$delta < 1) { lbl <- paste(lbl,"Constantly decreasing recruitment") }
      else { lbl <- paste(lbl,"Constantly increasing recruitment") }         
    }
  }    
  viol.sum <- sum(c(viol.Zdelta,viol.Nodelta))
  if(viol.sum==0) {
    lbl <- "No Assumptions Violations"
  } else if (viol.sum > 1) {
    lbl <- "VIOLATION: Multiple Assumptions"
  }
  if (Z.param$cv!=0 | No.param$cv!=0) {
    viol.rand <- TRUE
    lbl <- paste(lbl,", with randomness",sep="")
  }
  lbl
} # end iCC.title
