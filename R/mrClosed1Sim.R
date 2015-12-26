#' @title Dynamics plots to explore single-census mark-recapture models.
#'
#' @description This is used to simulate multiple results for single census mark-recapture study.  The user can control many aspects of the simulation including simulating the loss of marks, differential survival rates for tagged and untagged fish, the addition of recruits, and differential probabilities of capture for tagged and untagged fish.  In addition, Petersen, Chapman, and Bailey estimates of population size can be made.
#'
#' @details The user can use slider bars to choose values for the expected number of fish to capture in the first (i.e., the marking or tagging) sample, the expected number of fish to be captured in the second or final (i.e., the recapture) sample, a probablity that a tagged fish loses the tag, the probability of survival for tagged fish, the probability of survival for untagged fish, a proportion of the initial population size to recruit to the population between the first and final samples, and the ratio of the probability of capture of tagged fish to the probability of capture of untagged fish in the final sample.  The user selects the \dQuote{true} initial population size, the method of population estimator, and the number of resamples to use with arguments to the function.  The simulation can be re-run without changing any of the slider buttons by pressing the \bold{\sQuote{Re-Randomize}} button.
#'
#' In general, the simulation follows these steps:
#' \enumerate{
#'    \item a population of N fish is created;
#'    \item randomly select M fish in the first sample to be tagged such that the average of all M values should be equal to the user-supplied \bold{\sQuote{Tagged (M)}} value;
#'    \item randomly select m fish to lose tags according to the user-supplied probability of tag loss (\bold{\sQuote{PR(Tag Loss)}});
#'    \item randomly select fish to die according to the user-supplied probabilities of survival for tagged (\bold{\sQuote{PR(Surv Tagged)}}) and untagged fish (\bold{\sQuote{PR(Surv UNTagged)}});
#'    \item add recruits to the population in proportion to N and in accordance to the user-supplied proportion of recruits (\bold{\sQuote{Proportion Recruit}}) to add;
#'    \item identify the actual population size just before taking the final sample (N1);
#'    \item randomly select n fish in the second sample such that the average of all n values should be equal to the user-supplied \bold{\sQuote{Captured (n)}} value (see note below about how differential probabilties of capture are incorporated into the model);
#'    \item count the number of tagged fish in the second sample; and 
#'    \item compute the population estimate using the chosen method (see below).
#' }
#'
#' The methods for estimating the population size are
#'
#' \tabular{ll}{
#' \code{type="P"} \tab naive Petersen.\cr
#' \code{type="C"} \tab Chapman(1951) modification of the Petersen.\cr
#' \code{type="CR"} \tab Ricker(1975) modification of the Chapman modification.\cr
#' \code{type="B"} \tab Bailey(1951,1952) modification of the Petersen. 
#' }
#'
#' The effect of violating the assumption of tag loss is simulated by changing the probability of tag loss slider to a value greater than 0 but less than 1.  For example, setting \bold{\sQuote{PR(Tag Loss)}} to 0.1 is used to simulate a 10 percent probability of losing the tag.
#'
#' The effect of violating the assumption of no mortality from the first to second sample is simulated by changing one or both of the probabilities of survival for tagged and untagged fish to values less than 1 (but greater than 0).  For example, setting \bold{\sQuote{PR(Surv Tagged)}} AND \bold{\sQuote{PR(Surv UNTagged)}} to 0.8 will simulate mortality between the first and final sample but NOT differential mortality between the tagged and untagged fish (i.e., the mortalities are the same for both groups of fish).  The effect of differential mortalities between tagged and untagged fish can be simulated by using different survival probabilities for tagged and untagged fish.
#'
#' The effect of violating the assumption of no recruitment from the first to final sample is simulated by changing the Proportion Recruit slider to a value greater than 0.  For example, setting \bold{\sQuote{Proportion Recruit}} to 0.1 will simulate 10 percent of N recruiting to the population just before the final sample.
#'
#' The effect of violating the assumption of equal catchabilities in the final sample for tagged and untagged fish is simulated by changing the \bold{\sQuote{Ratio PR(Capture)}} slider to a value different than 1.  Values greater than 1 indicate that the catchability of tagged fish is greater than the catchability of untagged fish.  Values less than 1 indicate that the catchability of tagged fish is less than that of untagged fish.  For example, setting \bold{\sQuote{Ratio PR(Capture)}} to 0.8 will simulate a situation where the capture probability of tagged fish is 80 percent of the capture probablity of untagged fish (i.e., simulates the situation where marked fish are less likely to be captured).
#'
#' The probability of capture in the final sample is equal to the expected number of fish to be collected in the final sample (set with \bold{\sQuote{Captured (n)}}) divided by the actual population size just prior to the final sample (N1) as long as \bold{\sQuote{Ratio PR(Capture)}} is 1.  The probabilities of capture for the tagged and untagged fish are carefully adjusted if the \bold{\sQuote{Ratio PR(Capture)}} value is different than 1.  Because of the different numbers of tagged and untagged individuals in the population, the probability of capture for tagged and untagged individuals must be computed by adjusting the overall probability of capture to assure that, on average, the user-provided expected number captured in the final sample is met.  This modification is found by solving the following system of equations for the probabilities of capture for the tagged (PM) and untagged fish (PU), respectively,
#'
#' \deqn{PM*M + PU*(N-M) = P*N}
#'
#' \deqn{\frac{PM}{PU} = k}
#'
#' where M is the number of tagged animals, and N-M is the number of untagged animals, k is the \bold{\sQuote{Ratio PR(Capture)}} value, and P is the overall probability of capture if there was no difference in catchability between tagged and untagged animals.  The solutions to this system, which are used in this function, are
#'
#' \deqn{PU = \frac{PN}{M+N}}
#'
#' \deqn{PM = PU*k}
#'
#' @param type A single string that identifies the type of calculation method to use (see details).
#' @param N A single number that represents the \dQuote{known} size of the simulated population
#'just prior to the first sample.
#' @param rsmpls A single number that indicates the number of simulations to run.
#' @param incl.final A single logical that indicates whether the mean final population size and histogram of percent error from the final population size should be shown.
#'
#' @return None.  However, a dynamic graphic is produced that is controlled by slider bars as described in the details.  The dynamic graphic is a histogram of the population estimate from all resamples with a red vertical line at the initial population size (provided by the user), a blue vertical line at the population size just prior to the final sample (N1; if \code{incl.final=TRUE}), and a green vertical line at the mean population estimate.  The vertical blue and red lines may not be visible under some scenarios because of overprinting.
#'
#' @author Derek H. Ogle, \email{dogle@@northland.edu}
#'
#' @seealso \code{mrClosed} in \pkg{FSA}.
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
mrClosed1Sim <- function(sim=c("assumptions","distribution"),
                         N=1000,rsmpls=2000,EM=200,En=200,incl.final=TRUE) {
  if (!iCheckRStudio()) stop("'mrClosed1Sim()' only works in RStudio.",call.=FALSE)
  sim <- match.arg(sim)
  if (sim=="assumptions") iMRC1Assump(N,rsmpls,EM,En,incl.final)
    else iMRC1Dist(N,rsmpls,incl.final)
}


##############################################################

## Internal function for the assumption violations simulation
iMRC1Assump <- function(N,rsmpls,EM,En,incl.final) {
  if (iChk4Namespace("manipulate")) {
    manipulate::manipulate(
      {
        set.seed(sample(1:100000))
        # Run population simulations
        mc.df <- iMRC1.genpopn(type,N,EM,En,mark.loss,surv.mark,surv.unmark,recruits,cap.ratio,rsmpls)
        # Find label for graph
        hlbl <- iMRC1.title(mark.loss,surv.mark,surv.unmark,recruits,cap.ratio)
        # Graph the results
        iMRC1.hist(mc.df,N,incl.final,hlbl)  
        # Add legend
        iMRC1.legend("assumptions",mc.df,N,incl.final)
      },
      type=manipulate::picker("Petersen","Chapman","Ricker","Bailey",label="Method"),
      mark.loss=manipulate::slider(min=0,max=0.9,step=0.1,initial=0,label="PR(Tag Loss)"),
      surv.mark=manipulate::slider(min=0.1,max=1,step=0.1,initial=1,label="PR(Surv Tagged)"),
      surv.unmark=manipulate::slider(min=0.1,max=1,step=0.1,initial=1,label="PR(Surv UNtagged)"),
      recruits=manipulate::slider(min=0,max=0.5,step=0.1,initial=0,label="Proportion Recruit"),
      cap.ratio=manipulate::slider(min=0.5,max=2,step=0.1,initial=1,label="Ratio PR(Capture)"),
      rerand=manipulate::button("Rerandomize")
    )
  }
}

## Internal function for the distribution simulation
iMRC1Dist <- function(N,rsmpls,incl.final) {
  if (iChk4Namespace("manipulate")) {
    N.1 <- round(0.01*N)
    manipulate::manipulate(
      {
        set.seed(sample(1:100000))
        # Run population simulations
        mc.df <- iMRC1.genpopn(type,N,EM,En,0,1,1,0,1,rsmpls)
        # Graph the results
        iMRC1.hist(mc.df,N,incl.final,"")  
        # Add legend
        iMRC1.legend("distribution",mc.df,N,incl.final)
        
      },
      type=manipulate::picker("Petersen","Chapman","Ricker","Bailey",label="Method"),
      EM=manipulate::slider(min=5*N.1,max=40*N.1,step=N.1,initial=20*N.1,label="Tagged (M)"),
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
    adj.M <- M[i] <- length(which(runif(N)<(EM/N)))
    # Unmarked fish
    adj.U <- N-M[i]
    if (mark.loss>0) {
      # Apply tag loss probability
      lost.marks <- length(which(runif(M[i])<mark.loss))
      # Removed lost mark fish from marked popn
      adj.M <- adj.M - lost.marks
      # Put lost mark fish back into unmarked popn
      adj.U <- adj.U + lost.marks
    }
    # Marked fish survival -- greater than survival prob is mortality
    if (!surv.mark==1)   adj.M <- adj.M - length(which(runif(M[i])>surv.mark))
    # UnMarked fish survival
    if (!surv.unmark==1) adj.U <- adj.U - length(which(runif(adj.U)>surv.unmark))
    # Add recruits to unmarked population
    if (!recruits==0) adj.U <- adj.U + N*recruits
    # Population size just before second capture
    N1[i] <- adj.M + adj.U
    # Probability of capturing unmarked fish in second sample
    p.unmark <- En/N1[i]*N1[i]/(M[i]+N1[i])
    # Probability of capturing marked fish in second sample
    p.mark <- p.unmark*cap.ratio
    # Number of marked fish recaptured in second sample
    m[i] <- length(which(runif(adj.M)<p.mark))
    # Total number of fish captured in second sample
    n[i] <- m[i] + length(which(runif(adj.U)<p.unmark))
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
  if(mark.loss > 0) { viol.markloss <- TRUE; lbl <- "VIOLATION: Loss of Tags" }
  if(surv.mark < 1 || surv.unmark < 1) {
    if (surv.mark == surv.unmark) { viol.surv <- TRUE ; lbl <- "VIOLATION: Mortality" }
    else if (surv.unmark == 1) { viol.surv <- TRUE; lbl <- "VIOLATION: Mortality of Tagged Fish" }
    else if (surv.mark ==1) { viol.surv <- TRUE; lbl <- "VIOLATION: Mortality of UNtagged Fish" }
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
iMRC1.hist <- function(df,N,incl.final,hlbl) {
  # Set the graphing parameters
  old.par <- graphics::par(mar=c(3.5,1.1,1.5,1.1),mgp=c(2,0.4,0),tcl=-0.2,yaxs="i")
  # Make the histogram
  h <- graphics::hist(df$N0,plot=FALSE,right=FALSE)
  graphics::hist(df$N0,right=FALSE,main=hlbl,xlab="Population Estimate",ylab="",yaxt="n",
                 xlim=range(c(mean(df$N1),mean(df$N0),N,h$breaks)),col="gray90")
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
