#' @title Dynamic plots to explore the Leslie depletion model.
#'
#' @description Constructs hypothetical catch and effort data given choices for number of removal events, initial population size, effort, catchability, survival, and recruitment for a hypothetical depletion fishery. Various plots are produced (see details) with corresponding slider bars that allow the user to modify parameter values to explore the effect of these modifications on the Leslie model dynamics or results. This function is used primarily to explore the effects of a parameter on the model and the effects of assumptions violations on modelydynamics.
#'
#' @details Three versions of the simulation are allowed. First, the model is purely \bold{\code{deterministic}} (i.e., without randomization) such that the shape of the Leslie model can be easily explored. Second, the model includes a \bold{\code{random}} component and one collection of data is plotted with the Leslie model and the estimates of the catchabitilty coefficient (q) and initial population size (No) superimposed. This version allows the user to explore the effect of variability on the model dynamics and how this variability and assumption violations effect the parameter estimates. The third version is like the second except that the results are repeated \code{rsmpls} times, the q and No are computed from each resample, and these results are plotted. This version, called the \bold{\code{resampling}} or \bold{\code{montecarlo}} version allows the user to examine the sampling distributions of the parameter estimates and compare those to the known parameters to explore the bias caused by assumption violations. The plots produced are described further in the \dQuote{returns} section below.
#'
#' In the \code{random} and \code{resampling} versions, randomness is included in the model by including binomial stochasticity in the catch and survival functions. Specifically, the number captured is, effectively, computed by assigning a uniform random number from between 0 and 1 to each individual in the population and then \dQuote{catching} those individuals where this value is less than q*E (where E is effort expended). A similar method is used for survivorship, but with q*E replaced with the user chosen probability of survival.
#
#' For each version, a plot or plots is produced that is linked to slider bars that allows the user to change model parameters or create assumption violations. A slider is created to control the initial population size (No), effort (E), and catchability coefficent (q). The remaining sliders allow for simulating specific violations to the assumption of a Leslie model. These sliders are described further below.
#'
#' The \bold{\sQuote{q factor}} value is a constant that modifies the catchability coefficient (q) for each subsequent sample. For example, if \bold{\sQuote{q.factor}} is set to 0.8 then the catchability decreases by a constant multiplier of 0.8 for each sample. In other words, the catchability set with the catchability slider is multiplied by the vector \code{c(1,0.8,0.8^2,0.8^3,...)} to determine a catchability for each removal event.
#'
#' The \bold{\sQuote{Survival}} value is a constant used as a proportion of fish alive at time t that survive to time t+1 or, in the the \code{random} and \code{resampling} versions, is the probability that a fish survives from time t to time t+1. The survival function is applied to the population after the catch at time t has already been removed from the population.
#'
#' The \bold{\sQuote{Recruitment}} value is a constant used to determine the number of \dQuote{new} fish to recruit to the population from time t to time t+1. The number to recruit is equal to the recruitment proportion of the extant number of fish alive at time t. For example, if 100 fish are alive at time t and the recruitment factor is 0.1 then 100*0.1=10 fish will be added to the population just before time t+1. The number of fish to recruit is computed after the catch at time t and any natural mortality at time t have been removed from the population.
#'
#' @note The range of values allowed for each of the parameters were chosen to allow a wide variety of model values. However, it is highly likely that these ranges do not encompass every possible set of values that a user may wish to view. Thus, this simulation should not be used for research-grade simulations.
#'
#' @param sim A single string that indicates the type or version of simulation that should be used. See the details.
#' @param removals A single numeric that indicates the number of removal events to simulate.
#' @param Ricker.mod A single logical value that indicates whether the modification proposed by Ricker should be used (\code{=TRUE}) or not (\code{=FALSE}, default).
#' @param rsmpls A single numeric for the number of simulations to run.
#'
#' @return None. An interactive graphic with corresponding slider bars, which differ depending on the version (as defined by \code{sim}) of simulation used, is produced.
#'
#' In the \code{deterministic} and \code{random} versions a plot of catch-per-unit-effort (CPE) against cumulative catch (i.e., the Leslie plot) is displayed. In the \code{deterministic} version, as many as three lines may be seen. The gray line is the Leslie model for the default values from the slider bars. This line is used simply as a basis for examining changes in parameters. The blue line is the Leslie model for current choices of \bold{\sQuote{Initial Size}}, \bold{\sQuote{Effort}}, and \bold{\sQuote{Catchability}}, but NOT for \bold{\sQuote{q factor}}, \bold{\sQuote{Survival}}, or \bold{\sQuote{Recruitment}}. In other words, the blue line reflects the model for other than default parameter choices but with NO assumption violations. This line serves as a basis for judging different parameter choices without any assumption violations. The red line is the Leslie model for all current choices of sliders. The lines are plotted in the order of \dQuote{gray}, \dQuote{red}, \dQuote{blue} so, if any two are equal then the color first plotted will not be seen.
#'
#' In the \code{random} version, the graphic is simply the traditional Leslie model graphic (see \code{\link[FSA]{depletion}}) with the \dQuote{random} catch-per-unit-effort values plotted against total catch with a best-fit linear regression line shown in blue. The current estimaes of q and No from the random data are also printed on the graph. A \bold{\sQuote{Rerandomize}} button is included with the sliders which can be used to evaluate the model again (with a different random seed) at the current slider choices.
#'
#' In the \code{resampling} version, a multipaned plot will be produced. The scatterplot is of the paired catchability and initial population size estimates with red lines showing the true values of the catchability and initial population size and blue lines at the means of the respective estimates. The top histogram is of the estimates of catchability from all resamples, whereas the right histogram is of the estimates of the initial population size (No) from all resamples. Both histograms will have a vertical red dashed line at the true value of the parameter (No or q, as provided by the user) and a vertical blue solid line at the mean value of the estimate from all resamples.
#'
#' @author Derek H. Ogle, \email{dogle@@northland.edu}
#'
#' @seealso \code{\link[FSA]{depletion}} in \pkg{FSA}.
#'
#' @keywords iplot
#'
#' @examples
#' if (interactive()) {
#'
#' ## Deterministic exploration of model dynamics
#' leslieSim()
#' 
#' ## Stochastic exploration of model dynamics -- Leslie model plot
#' leslieSim(type="random")
#'
#' ## Stochastic exploration of model dynamics -- sampling distribution plots
#' leslieSim(type="resampling")
#'
#'} # end if interactive
#'
#' @export
#'
leslieSim <- function(sim=c("deterministic","random","resampling","montecarlo"),
                      removals=8,Ricker.mod=FALSE,rsmpls=100) {
  sim <- match.arg(sim)
  switch(sim,
         deterministic={ iLeslie.DET(removals,Ricker.mod) } ,
         random={ iLeslie.RAND(removals,Ricker.mod) },
         resampling=,montecarlo={ iLeslie.MC(removals,Ricker.mod,rsmpls) }
         )
}


## Internal functions for DETerministic simulations
iLeslie.DET <- function(k,Ricker.mod) {
  # Trying to deal with no visible bindings problem
  No <- E <- q.adj.const <- p.surv <- r.prop <- NULL
  if (!iCheckRStudio()) FSA:::STOP("'leslieSim' only works in RStudio.")
  if (iChk4Namespace("manipulate")) {
    manipulate::manipulate(
      {
        # put default parameters for 8 removal events into a parameter list
        params.def <- list(k=k,Nt=rep(400,k),Et=rep(5,k),qt=rep(0.08,k),
                           p.surv=rep(1,k),r.prop=rep(0,k))
        # get values from sliders and put in list for no violations lines
        params.mod <- list(k=k,Nt=rep(No,k),Et=rep(E,k),qt=rep(q,k),
                           p.surv=rep(1,k),r.prop=rep(0,k))
        # get assumption violation values from sliders and put in a list
        params.ass <- list(k=k,Nt=rep(No,k),Et=rep(E,k),qt=q*q.adj.const^(0:(k-1)),
                           p.surv=rep(p.surv,k),r.prop=rep(r.prop,k))
        # Make the plot
        iLeslie.detPlotMain(params.def,params.mod,params.ass,q.adj.const,Ricker.mod)
      },
      No=manipulate::slider(100,800,step=50,initial=400,label="Initial Size (No)"),
      E=manipulate::slider(1,10,step=1,initial=5,label="Effort (E)"),
      q=manipulate::slider(0.01,0.1,step=0.01,initial=0.08,label="Catchability (q)"),
      q.adj.const=manipulate::slider(0.7,1.3,step=0.05,initial=1.0,label="q factor"),
      p.surv=manipulate::slider(0.7,1.0,step=0.05,initial=1.0,label="PR(Survival)"),
      r.prop=manipulate::slider(0,0.1,step=0.02,initial=0.0,label="Prop Recruitment")
    )
  }
}

iLeslie.detPlotMain <- function(params.def,params.mod,params.ass,
                                q.adj.const,Ricker.mod) {
  # Get results for different sets of parameter values
  res.def <- iLeslie.detRun(params.def)
  res.mod <- iLeslie.detRun(params.mod)
  res.ass <- iLeslie.detRun(params.ass)
  # Plot total results
  iLeslie.detPlot(res.def,res.mod,res.ass,Ricker.mod,
                  iLeslie.title(q.adj.const,params.ass$p.surv,params.ass$r.prop))
  graphics::legend("topright",
                   legend=c("assumptions not met","if assumptions met","default"),
                   lwd=2,col=c("red","blue","gray"),pch=19,bg="white",
                   box.col="white",inset=0.01)
} # end iLeslie.detPlotMain

iLeslie.detRun <- function(plist) {
  k <- plist$k
  Nt <- plist$Nt
  Et <- plist$Et
  qt <- plist$qt
  p.surv <- plist$p.surv
  r.prop <- plist$r.prop
  # Initialize catch, mortality, and recruitment vectors      
  Ct <- Mt <- Rt <- rep(0,k)
  for (i in 1:k) {
    Ct[i] <- qt[i]*Et[i]*Nt[i]
    # Update population size after catch for next pass
    if(i<k) {
      # Reduce by catch in current time
      Nt[i+1] <- Nt[i]-Ct[i]
      # Reduce by mortality, if any, in current time
      if (any(p.surv<1)) {
        Mt[i] <- (1-p.surv[i])*Nt[i+1]
        Nt[i+1] <- Nt[i+1]-Mt[i]
      }
      # Increase by recruitment, if any, in current time
      if (any(r.prop>0)) {
        Rt[i] <- round(Nt[i+1]*r.prop[i],0)
        Nt[i+1] <- Nt[i+1]+Rt[i]
      }
    }
  }
  data.frame(Nt,qt=round(qt,4),Mt,Rt,Et,Ct,Pt=round(Ct/Nt,3))
} # end iLeslie.detRun

iLeslie.detPlot <- function(resdef,res1,res2,Ricker.mod=FALSE,glbl) {
  resdef$cpe <- resdef$Ct/resdef$Et
  res1$cpe <- res1$Ct/res1$Et
  res2$cpe <- res2$Ct/res2$Et
  if (!Ricker.mod) {
    resdef$K <- cumsum(resdef$Ct)-resdef$Ct
    res1$K <- cumsum(res1$Ct)-res1$Ct
    res2$K <- cumsum(res2$Ct)-res2$Ct
  } else {
    resdef$K <- cumsum(resdef$Ct)-(resdef$Ct/2)
    res1$K <- cumsum(res1$Ct)-(res1$Ct/2)
    res2$K <- cumsum(res2$Ct)-(res2$Ct/2)
  }
  # Plot values for default parameters
  withr::local_par(mar=c(3.5,3.5,1.5,1.5),mgp=c(2,0.4,0),tcl=-0.2,
                      xaxs="i",yaxs="i")
  graphics::plot(cpe~K,data=resdef,type="o",pch=19,lwd=2,col="gray",main=glbl,
                 xlab="Cumulative Catch",xlim=c(0,1.75*max(resdef$K)),
                 ylab="CPE",ylim=c(0,1.75*max(resdef$cpe)))
  # Plot values for if assumptions are not met
  graphics::points(res2$K,res2$cpe,type="o",pch=19,lwd=2,col="red")
  # Plot values for modified parameters but if no assumptions violated
  graphics::points(res1$K,res1$cpe,type="o",pch=19,lwd=2,col="blue")
} # end iLeslie.detPlot





## Internal functions for RANDom simulations
iLeslie.RAND <- function(k,Ricker.mod) {
  # Trying to deal with no visible bindings problem
  No <- E <- q.adj.const <- p.surv <- r.prop <- NULL
  if (!iCheckRStudio()) FSA:::STOP("'leslieSim' only works in RStudio.")
  if (iChk4Namespace("manipulate")) {
    manipulate::manipulate(
      {
        # Make the plot
        params <- list(k=k,Nt=rep(No,k),Et=rep(E,k),qt=q*q.adj.const^(0:(k-1)),
                       p.surv=rep(p.surv,k),r.prop=rep(r.prop,k))
        set.seed(sample(1:100000))
        iLeslie.randPlotMain(params,Ricker.mod,q.adj.const)
      },
      No=manipulate::slider(100,800,step=50,initial=400,label="Initial Size (No)"),
      E=manipulate::slider(1,10,step=1,initial=5,label="Effort (E)"),
      q=manipulate::slider(0.01,0.1,step=0.01,initial=0.08,label="Catchability (q)"),
      q.adj.const=manipulate::slider(0.7,1.3,step=0.05,initial=1.0,label="q factor"),
      p.surv=manipulate::slider(0.7,1.0,step=0.05,initial=1.0,label="PR(Survival)"),
      r.prop=manipulate::slider(0,0.1,step=0.02,initial=0.0,label="Prop Recruitment"),
      rerand=manipulate::button("Rerandomize")
    )
  }
}

iLeslie.randPlotMain <- function(params,Ricker.mod,q.adj.const) {
  # Results for possible assumptions violations
  res.ass <- iLeslie.randRun(params)
  res <- FSA::depletion(res.ass$Ct,res.ass$Et,method="Leslie",Ricker.mod=Ricker.mod)
  withr::local_par(mar=c(3.5,3.5,1.5,1.5),mgp=c(2,0.4,0),tcl=-0.2)
  graphics::plot(res,xlim=c(0,800),ylim=c(0,50),
                 main=iLeslie.title(q.adj.const,params$p.surv,params$r.prop))
} # end iLeslie.randPlotMain

iLeslie.randRun <- function(params) {
  k <- params$k
  Nt <- params$Nt
  Et <- params$Et
  qt <- params$qt
  p.surv <- params$p.surv
  r.prop <- params$r.prop
  # Initialize catch, mortality, and recruitment vectors
  Ct <- Mt <- Rt <- rep(0,k)      
  for (i in 1:k) {
    # Determine whether caught or not -- capture probability is q*E*N/N or just q*E
    Ct[i] <- length(which(stats::runif(Nt[i])<qt[i]*Et[i]))
    # Update population size after catch for next pass
    if(i<k) {
      # Reduce by catch in current time
      Nt[i+1] <- Nt[i]-Ct[i]
      # Reduce by mortality, if any, in current time
      if (any(p.surv<1)) {
        # Mortality
        Mt[i] <- length(which(stats::runif(Nt[i+1])<(1-p.surv[i])))
        Nt[i+1] <- Nt[i+1]-Mt[i]
      }  
      # Increase by recruitment, if any, in current time
      if (any(r.prop>0)) {
        Rt[i] <- round(Nt[i+1]*r.prop[i],0)
        Nt[i+1] <- Nt[i+1]+Rt[i]
      }
    }
  }
  data.frame(Nt,qt=round(qt,4),Mt,Rt,Et,Ct,Pt=round(Ct/Nt,3))
} # end iLeslie.randRun


## Internal functions for Monte Carlo (resampling) simulations
iLeslie.MC <- function(k,Ricker.mod,rsmpls) {
  # Trying to deal with no visible bindings problem
  No <- E <- q.adj.const <- p.surv <- r.prop <- NULL
  if (!iCheckRStudio()) FSA:::STOP("'leslieSim' only works in RStudio.")
  if (iChk4Namespace("manipulate")) {
    manipulate::manipulate(
      {
        # Make the plot
        params <- list(k=k,Nt=rep(No,k),Et=rep(E,k),qt=q*q.adj.const^(0:(k-1)),
                       p.surv=rep(p.surv,k),r.prop=rep(r.prop,k))
        set.seed(sample(1:100000))
        iLeslie.mcPlotMain(params,Ricker.mod,rsmpls)
      },
      No=manipulate::slider(100,800,step=50,initial=400,label="Initial Size (No)"),
      E=manipulate::slider(1,10,step=1,initial=5,label="Effort (E)"),
      q=manipulate::slider(0.01,0.1,step=0.01,initial=0.08,label="Catchability (q)"),
      q.adj.const=manipulate::slider(0.7,1.3,step=0.05,initial=1.0,label="q factor"),
      p.surv=manipulate::slider(0.7,1.0,step=0.05,initial=1.0,label="PR(Survival)"),
      r.prop=manipulate::slider(0,0.1,step=0.02,initial=0.0,label="Prop Recruitment"),
      rerand=manipulate::button("Rerandomize")
    )
  }
}

iLeslie.mcPlotMain <- function(params,Ricker.mod,rsmpls) {
  res.qhat <- res.N0hat <- NULL
  for (i in 1:rsmpls) {
    # Results for possible assumptions violations
    res.ass <- iLeslie.randRun(params)
    if (!Ricker.mod) res <- FSA::depletion(res.ass$Ct,res.ass$Et,method="Leslie")
    else res <- FSA::depletion(res.ass$Ct,res.ass$Et,method="Leslie",Ricker.mod=Ricker.mod)
    res.qhat[i] <- res$est["q","Estimate"]
    res.N0hat[i] <- res$est["No","Estimate"]
  }
  ## histogram values
  Nohist <- graphics::hist(res.N0hat,right=FALSE,plot=FALSE)
  qhist <- graphics::hist(res.qhat,right=FALSE,plot=FALSE)
  ## axis ranges (make sure that truth is in range)
  Norng <- c(0.98,1.02)*range(c(params$Nt,res.N0hat,Nohist$breaks))
  qrng <- c(0.98,1.02)*range(c(params$qt[1],res.qhat,qhist$breaks))
  ## make the layout
  graphics::layout(matrix(c(2,4,1,3),nrow=2,byrow=TRUE),c(1.5,1),c(1,1.5),TRUE)
  ## Scatterplot of No vs q
  withr::local_par(mar=c(3,3,0,0),mgp=c(1.7,0.4,0),tcl=-0.2) 
  graphics::plot(res.qhat,res.N0hat,pch=16,cex=1.1,col=grDevices::rgb(0,0,0,0.8),
                 xlab="Catchability Estimate",xlim=qrng,
                 ylab="Estimated Population Size",ylim=Norng)
  # Add reference lines
  graphics::abline(h=c(params$Nt[1],mean(res.N0hat)),
                   col=c("red","blue"),lwd=c(3,2),lty=c(3,2))
  graphics::abline(v=c(params$qt[1],mean(res.qhat)),
                   col=c("red","blue"),lwd=c(3,2),lty=c(3,2))
  ## Histogram of catchability estimates
  graphics::par(mar=c(0,3,0.5,0),yaxs="i")
  graphics::plot(NULL,type="n",ylim=c(0,max(qhist$counts)),xlim=qrng,
                 axes=FALSE,xlab="",ylab="")
  graphics::rect(qhist$breaks[1:(length(qhist$breaks)-1)],0,
                 qhist$breaks[2:length(qhist$breaks)],qhist$counts,col="gray90")
  # Add reference lines
  graphics::abline(v=c(params$qt[1],mean(res.qhat)),
                   col=c("red","blue"),lwd=c(3,2),lty=c(3,2))
  graphics::abline(v=params$qt[1],col="red",lwd=3,lty=3)
  ## Histogram of population estimates
  graphics::par(mar=c(3,0,0,0.5),xaxs="i")
  graphics::plot(NULL,type="n",xlim=c(0,max(Nohist$counts)),ylim=Norng,
                 axes=FALSE,xlab="",ylab="")
  graphics::rect(0,Nohist$breaks[1:(length(Nohist$breaks)-1)],Nohist$counts,
                 Nohist$breaks[2:length(Nohist$breaks)],col="gray90")
  # Add reference lines
  graphics::abline(h=c(params$Nt[1],mean(res.N0hat)),
                   col=c("red","blue"),lwd=c(3,2),lty=c(3,2))
  ## Make legend
  # find % errors
  mn.N0.perr <- paste0("No % Error = ",round(mean(100*(res.N0hat-params$Nt[1])/params$Nt[1]),1))
  mn.q.perr <- paste0("q % Error = ",round(mean(100*(res.qhat-params$qt[1])/params$qt[1]),1))
  # place legend
  graphics::plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=FALSE)
  graphics::legend("topleft",c("Known","Mean",NA,mn.N0.perr,mn.q.perr),
                   cex=1.5,bty="n",col=c("red","blue",NA,NA,NA),
                   lwd=c(3,2,NA,NA,NA),lty=c(3,2,NA,NA,NA))
  ## Return to original plotting parameters
  graphics::layout(1)
} # end iLeslie.mcPlotMain


iLeslie.title <- function(q.adj.const,p.surv,r.prop) {
  viol.q <- viol.surv <- viol.recruits <- FALSE
  if (q.adj.const != 1) {
    viol.q <- TRUE
    if (q.adj.const < 1) lbl <- "VIOLATION: decreasing catchability"
    else lbl <- "VIOLATION: increasing catchability"
  }
  if (any(p.surv < 1)) {
    viol.surv <- TRUE
    lbl <- "VIOLATION: natural mortality is occurring"
  }
  if (any(r.prop > 0)) {
    viol.recruits <- TRUE
    lbl <- "VIOLATION: recruitment is occurring"
  }
  viol.sum <- sum(c(viol.q,viol.surv,viol.recruits))
  if(viol.sum==0) lbl <- "No Assumptions Violations"
  else if (viol.sum > 1) lbl <- "VIOLATION: Multiple Assumptions"
  lbl
} # end iLeslie.title internal function
