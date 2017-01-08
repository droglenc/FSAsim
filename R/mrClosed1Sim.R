#' @title Dynamic plots to explore single-census mark-recapture models.
#'
#' @description This simulates multiple results for a single census mark-recapture study.  The user can control aspects of the sampling including the number of marked animals (animals in the first sample) and the number of animals examined for marks (animals in the second sample) and aspects of the population including loss of marks on marked animals, survival rates for marked and unmarked fish, the addition of recruits, and differential probabilities of capture for marked and unmarked fish.  In addition, Petersen, Chapman, Ricker, and Bailey estimates of population size can be made (see \code{\link[FSA]{FSA}}).  \bold{Function can only be used with RStudio.}
#'
#' @details Two types of simulations can be conducted.  The first is primarily used to examine characteristics of the sampling distrubtion of the population estimates.  These simulations are selected with \code{sim="distribution"} and allow the user to dynamically select the expected number of fish to capture in the first (i.e., the marking or tagging) sample and the expected number of fish to be captured in the second (i.e., the recapture) sample.  The resulting plot (discussed more below) allows the user to examine the shape of the sampling distribution, see the width of chosen quantiles on the distribution (defaults to 95%), and the bias from the initial population size.
#' 
#' The second type of simulation is primarily used to examine the effect of assumption violations on the bias of the population estimates.  These simulations are selected with \code{sim="assumptions"} and allow the user to dynamically select the probablity that a marked fish loses the mark, the probability of survival for marked fish, the probability of survival for unmarked fish, a proportion of the initial population size to recruit to the population between the first and second samples, and the ratio of the probability of capture of marked fish to the probability of capture of unmarked fish in the second sample.
#' 
#' In both types of simulations, the user can dynamically select the type of population estimation method (Petersen, Chapman, Ricker, or Bailey; see \code{\link[FSA]{FSA}}) and rerun the simulation without changing any of the other items by selecting the \code{Rerandomize} button.  The user can also select the population size with the \code{N=} argument and the number of resamples use to construct the sampling distribution with the \code{rsmpls=} argument.  In the simulations to assess the assumptions the user can also set the the expected number of fish to capture in the first sample with the \code{EM} argument and the expected number of fish to be captured in the second sample with the \code{En} argument.
#'
#' In general, the simulation follows these steps:
#' \enumerate{
#'    \item a population of N fish is created;
#'    \item randomly select M fish in the first sample to be marked such that the average of all M values should be equal to the user-supplied \bold{\sQuote{Marked (M)}} (first type of simulations) or \code{EM} (second type of simulations) value;
#'    \item randomly select m fish to lose marks according to the user-supplied probability of mark loss (\bold{\sQuote{PR(Mark Loss)}});
#'    \item randomly select fish to die according to the user-supplied probabilities of survival for marked (\bold{\sQuote{PR(Surv Marked)}}) and unmarked fish (\bold{\sQuote{PR(Surv UNMarked)}});
#'    \item add recruits to the population in proportion to N and in accordance to the user-supplied proportion of recruits (\bold{\sQuote{Proportion Recruit}}) to add;
#'    \item identify the actual population size just before taking the final sample (N1);
#'    \item randomly select n fish in the second sample such that the average of all n values should be equal to the user-supplied \bold{\sQuote{Captured (n)}} (first type of simulations) or \code{En} (second type of simulations) value (see note below about how differential probabilties of capture are incorporated into the model);
#'    \item count the number of marked fish in the second sample; and 
#'    \item compute the population estimate using the chosen method (see below).
#' }
#'
#' The methods for estimating the population size are
#'
#' \tabular{ll}{
#' \code{type="Petersen"} \tab naive Petersen.\cr
#' \code{type="Chapman"} \tab Chapman(1951) modification of the Petersen.\cr
#' \code{type="Ricker"} \tab Ricker(1975) modification of the Chapman modification.\cr
#' \code{type="Bailey"} \tab Bailey(1951,1952) modification of the Petersen. 
#' }
#'
#' The effect of violating the assumption of mark loss is simulated by changing the probability of mark loss slider to a value greater than 0 but less than 1.  For example, setting \bold{\sQuote{PR(Mark Loss)}} to 0.1 is used to simulate a 10 percent probability of losing the mark.
#'
#' The effect of violating the assumption of no mortality from the first to second sample is simulated by changing one or both of the probabilities of survival for marked and unmarked fish to values less than 1 (but greater than 0).  For example, setting \bold{\sQuote{PR(Surv Marked)}} AND \bold{\sQuote{PR(Surv UNmarked)}} to 0.8 will simulate mortality between the first and final sample but NOT differential mortality between the marked and unmarked fish (i.e., the mortalities are the same for both groups of fish).  The effect of differential mortalities between marked and unmarked fish can be simulated by using different survival probabilities for marked and unmarked fish.
#'
#' The effect of violating the assumption of no recruitment from the first to second sample is simulated by changing the \bold{\sQuote{Proportion Recruit}} slider to a value greater than 0.  For example, setting \bold{\sQuote{Proportion Recruit}} to 0.1 will simulate 10 percent of N recruiting to the population just before the second sample.
#'
#' The effect of violating the assumption of equal catchabilities in the second sample for marked and unmarked fish is simulated by changing the \bold{\sQuote{PR(Capture) Ratio (M/U)}} slider to a value different than 1.  Values greater than 1 indicate that the catchability of marked fish is greater than the catchability of unmarked fish.  Values less than 1 indicate that the catchability of marked fish is less than that of unmarked fish.  For example, setting \bold{\sQuote{PR(Capture) Ratio (M/U)}} to 0.8 will simulate a situation where the capture probability of marked fish is 80 percent of the capture probablity of unmarked fish (i.e., simulates the situation where marked fish are less likely to be captured in the second sample than unmarked fish).
#'
#' The probability of capture in the final sample is equal to the expected number of fish to be collected in the final sample (set with \bold{\sQuote{Captured (n)}} or \code{En}) divided by the actual population size just prior to the final sample (N1) as long as \bold{\sQuote{PR(Capture) Ratio (M/U)}} is 1.  The probabilities of capture for the marked and unmarked fish are carefully adjusted if the \bold{\sQuote{PR(Capture) Ratio (M/U)}} value is different than 1.  Because of the different numbers of marked and unmarked individuals in the population, the probability of capture for marked and unmarked individuals must be computed by adjusting the overall probability of capture to assure that, on average, the user-provided expected number captured in the final sample is met.  This modification is found by solving the following system of equations for the probabilities of capture for the marked (PM) and unmarked (PU) fish, respectively,
#'
#' \deqn{PM*M + PU*(N-M) = P*N}
#'
#' \deqn{\frac{PM}{PU} = k}
#'
#' where M is the number of marked animals, and N-M is the number of unmarked animals, k is the \bold{\sQuote{PR(Capture) Ratio (M/U)}} value, and P is the overall probability of capture if there was no difference in catchability between marked and unmarked animals.  The solutions to this system, which are used in this function, are
#'
#' \deqn{PU = \frac{PN}{M+N}}
#'
#' \deqn{PM = PU*k}
#'
#' @param sim A single string that identifies the type of simulation to perform (see details).
#' @param N A single number that represents the known size of the simulated population just prior to the first sample.
#' @param rsmpls A single number that indicates the number of simulations to run.
#' @param EM A single number that indicates the expected number of fish to capture (and mark) in the first sample.  Only used if \code{sim="assumptions"}.
#' @param En A single number that indicates the expected number of fish to capture (and check for marks) in the second sample.  Only used if \code{sim="assumptions"}.
#' @param incl.final A single logical that indicates whether the mean final population size should be shown on the plot produce (only used if \code{sim="assumptions"}).
#' @param conf.level A single number that indicates the quantiles to compute and label on the x-axis of the histogram when \code{sim="distribution"}.
#'
#' @return None.  However, a dynamic graphic is produced that is controlled by slider bars as described in the details.  The dynamic graphic is a histogram of the population estimate from all resamples with a red vertical line at the initial population size (provided by the user), a blue vertical line at the population size just prior to the final sample (N1; if \code{incl.final=TRUE}), and a green vertical line at the mean population estimate.  The vertical blue and red lines may not be visible under some scenarios because of overprinting.  The x-axis will be labelled with only the quantiles computed with \code{conf.level} when \code{sim="distrib"}.
#'
#' @author Derek H. Ogle, \email{dogle@@northland.edu}
#'
#' @seealso \code{\link[FSA]{mrClosed}}.
#'
#' @keywords iplot
#'
#' @examples
#' if (interactive()) {
#'   mrClosed1Sim()
#'   mrClosed1Sim(type="C")  # use Chapman modification
#' } # end if interactive
#'
#' @export
#'
mrClosed1Sim <- function(sim=c("assumptions","distribution"),N=1000,rsmpls=5000,
                         EM=200,En=200,incl.final=TRUE,conf.level=0.95) {
  sim <- match.arg(sim)
  if (sim=="assumptions") iMRC1Assump(N,rsmpls,EM,En,incl.final)
    else iMRC1Dist(N,rsmpls,conf.level)
}


##############################################################

## Internal function for the assumption violations simulation
iMRC1Assump <- function(N,rsmpls,EM,En,incl.final) {
  # Trying to deal with no visible bindings problem
  type <- mark.loss <- surv.mark <- surv.unmark <- recruits <- cap.ratio <- NULL
  if (!iCheckRStudio()) stop("'mrClosed1Sim' only works in RStudio.",call.=FALSE)
  if (iChk4Namespace("manipulate")) {
    manipulate::manipulate(
      {
        set.seed(sample(1:100000))
        # Run population simulations
        mc.df <- iMRC1.genpopn(type,N,EM,En,mark.loss,surv.mark,surv.unmark,recruits,cap.ratio,rsmpls)
        # Find label for graph
        hlbl <- iMRC1.title(mark.loss,surv.mark,surv.unmark,recruits,cap.ratio)
        # Graph the results
        iMRC1.hist(mc.df,N,incl.final,hlbl,NULL,NA)  
        # Add legend
        iMRC1.legend("assumptions",mc.df,N,incl.final)
      },
      type=manipulate::picker("Chapman","Petersen","Ricker","Bailey",label="Method"),
      mark.loss=manipulate::slider(min=0,max=0.9,step=0.1,initial=0,label="PR(Mark Loss)"),
      surv.mark=manipulate::slider(min=0.1,max=1,step=0.1,initial=1,label="PR(Surv Marked)"),
      surv.unmark=manipulate::slider(min=0.1,max=1,step=0.1,initial=1,label="PR(Surv UNmarked)"),
      recruits=manipulate::slider(min=0,max=0.5,step=0.1,initial=0,label="Proportion Recruit"),
      cap.ratio=manipulate::slider(min=0.5,max=2,step=0.1,initial=1,label="PR(Capture) Ratio (M/U)"),
      rerand=manipulate::button("Rerandomize")
    )
  }
}

## Internal function for the distribution simulation
iMRC1Dist <- function(N,rsmpls,conf.level) {
  # Trying to deal with no visible bindings problem
  type <- EM <- En <- NULL
  if (!iCheckRStudio()) stop("'mrClosed1Sim' only works in RStudio.",call.=FALSE)
  if (iChk4Namespace("manipulate")) {
    N.1 <- round(0.01*N)
    manipulate::manipulate(
      {
        set.seed(sample(1:100000))
        # Run population simulations
        mc.df <- iMRC1.genpopn(type,N,EM,En,0,1,1,0,1,rsmpls)
        # Graph the results
        iMRC1.hist(mc.df,N,FALSE,"","intervals",conf.level)
        # Add legend
        iMRC1.legend("distribution",mc.df,N,FALSE)
      },
      type=manipulate::picker("Petersen","Chapman","Ricker","Bailey",label="Method"),
      EM=manipulate::slider(min=5*N.1,max=40*N.1,step=N.1,initial=20*N.1,label="Marked (M)"),
      En=manipulate::slider(min=5*N.1,max=40*N.1,step=N.1,initial=20*N.1,label="Captured (n)"),
      rerand=manipulate::button("Rerandomize")
    )
  }
}


## Internal function to create and run the population simulations
iMRC1.genpopn <- function(type,N,EM,En,mark.loss,surv.mark,surv.unmark,recruits,cap.ratio,rsmpls) {
  # Initialize vector sizes
  M <- n <- m <- N0 <- N1 <- rep(0,rsmpls)
  for (i in 1:rsmpls) {
    # Number marked in first sample (EM/N is probability of being marked)
    adj.M <- M[i] <- length(which(stats::runif(N)<(EM/N)))
    # Unmarked fish
    adj.U <- N-M[i]
    if (mark.loss>0) {
      # Apply mark loss probability
      lost.marks <- length(which(stats::runif(M[i])<mark.loss))
      # Removed lost mark fish from marked popn
      adj.M <- adj.M - lost.marks
      # Put lost mark fish back into unmarked popn
      adj.U <- adj.U + lost.marks
    }
    # Marked fish survival -- greater than survival prob is mortality
    if (!surv.mark==1)   adj.M <- adj.M - length(which(stats::runif(M[i])>surv.mark))
    # UnMarked fish survival
    if (!surv.unmark==1) adj.U <- adj.U - length(which(stats::runif(adj.U)>surv.unmark))
    # Add recruits to unmarked population
    if (!recruits==0) adj.U <- adj.U + N*recruits
    # Population size just before second capture
    N1[i] <- adj.M + adj.U
    # Probability of capturing unmarked fish in second sample
    p.unmark <- En/N1[i]*N1[i]/(M[i]+N1[i])
    # Probability of capturing marked fish in second sample
    p.mark <- p.unmark*cap.ratio
    # Number of marked fish recaptured in second sample
    m[i] <- length(which(stats::runif(adj.M)<p.mark))
    # Total number of fish captured in second sample
    n[i] <- m[i] + length(which(stats::runif(adj.U)<p.unmark))
    # Compute population estimates
    switch(type,
           P=,Petersen={ N0[i] <- round((M[i]*n[i])/m[i],0) },
           C=,Chapman={ N0[i] <- round((M[i]+1)*(n[i]+1)/(m[i]+1)-1,0) },
           R=,Ricker={ N0[i] <- round((M[i]+1)*(n[i]+1)/(m[i]+1),0) },
           B=,Bailey={ N0[i] <- round(M[i]*(n[i]+1)/(m[i]+1),0) }
    )
  }
  data.frame(M,n,m,N,N1,N0)
}

## Internal function to make a main label for the plot
iMRC1.title <- function(mark.loss,surv.mark,surv.unmark,recruits,cap.ratio) {
  # Initialize
  viol.markloss <- viol.diffmort <- viol.surv <- viol.recruits <- viol.capratio <- FALSE
  # Identify the violation type and make appropriate label
  if(mark.loss > 0) { viol.markloss <- TRUE; lbl <- "VIOLATION: Loss of Marks" }
  if(surv.mark < 1 || surv.unmark < 1) {
    if (surv.mark == surv.unmark) { viol.surv <- TRUE ; lbl <- "VIOLATION: Mortality" }
    else if (surv.unmark == 1) { viol.surv <- TRUE; lbl <- "VIOLATION: Mortality of Marked Fish" }
    else if (surv.mark ==1) { viol.surv <- TRUE; lbl <- "VIOLATION: Mortality of UNmarked Fish" }
    else { viol.surv <- TRUE; lbl <- "VIOLATION: Mortality (Differential)" }
  }
  if(recruits > 0) { viol.recruits <- TRUE; lbl <- "VIOLATION: Recruitment" }
  if(cap.ratio > 1) { viol.capratio <- TRUE; lbl <- "VIOLATION: Trap-Happy" } 
  if(cap.ratio < 1) { viol.capratio <- TRUE; lbl <- "VIOLATION: Trap-Shy" }
  viol.sum <- sum(c(viol.markloss,viol.diffmort,viol.surv,viol.recruits,viol.capratio))
  if(viol.sum==0) lbl <- "No Assumption Violations"
  else if (viol.sum > 1) lbl <- "VIOLATION: Multiple Assumptions"
  lbl
}

## Internal function to make the main histogram
iMRC1.hist <- function(df,N,incl.final,hlbl,xaxis,conf.level) {
  # Set the graphing parameters
  old.par <- graphics::par(mar=c(3.5,1.1,1.5,1.1),mgp=c(2,0.4,0),tcl=-0.2,yaxs="i")
  # Make the histogram
  h <- graphics::hist(df$N0,plot=FALSE,breaks=20,right=FALSE)
  graphics::hist(df$N0,breaks=20,right=FALSE,col="gray90",main=hlbl,
                 yaxt="n",ylab="",
                 xaxt="n",xlab="Population Estimate",
                 xlim=range(c(mean(df$N1),mean(df$N0),N,h$breaks)))
  # Handle x-axis
  if (is.null(xaxis)) graphics::axis(1)
  else graphics::axis(1,c(N,round(stats::quantile(df$N0,0.5+c(-1,1)*conf.level/2),0)),lwd=3)
  # Put vertical line for set initial pop
  graphics::abline(v=N,col="red",lwd=4,lty=2)
  # Put vertical line for mean estimate of initial pop
  graphics::abline(v=mean(df$N0),col="green",lwd=4)
  # Put vertical line for mean estimate of final pop (if requested by user)
  if (incl.final) graphics::abline(v=mean(df$N1),col="blue",lwd=4,lty=3)
  # Return to original graphing parameters
  graphics::par(old.par)
}

# Internal function to add the legend to the plot
iMRC1.legend <- function(sim,df,N,incl.final) {
  # Make legend labels
  mn.N0hat <- paste0("Mean Pop Est (=",round(mean(df$N0),0),")")
  init.N <- paste0("Initial Pop (=",N,")")
  mn.N0.perr <- paste0("  % Error = ",round(mean(100*(df$N0-N)/N),1))
  mn.N1 <- paste0("Mean Final Pop (=",round(mean(df$N1),0),")")
  mn.N1.perr <- paste0("  % Error = ",round(mean(100*(df$N0-df$N1)/df$N1),1))
  # Make legend values
  legs <- c(mn.N0hat,NA,init.N,mn.N0.perr,NA,mn.N1,mn.N1.perr)
  cols <- c("green",NA,"red",NA,NA,"blue",NA)
  ltys <- c(1,NA,2,NA,NA,3,NA)
  # Add legend
  if (incl.final) {
    graphics::legend("topright",legend=legs,col=cols,lty=ltys,lwd=2,box.col="white",bg="white")
  } else {
    graphics::legend("topright",legend=legs[-c(4:6)],col=cols[-c(4:6)],lty=ltys[-c(4:6)],
                     lwd=2,box.col="white",bg="white")
  }
}
