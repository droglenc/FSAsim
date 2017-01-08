#'Simulate sampling fish based on length selectivity.
#'
#'Constructs a sample of fish based on a user-supplied length-based selectivity
#'curve.  The selectivity curve can be supplied via a parametric model (the
#'beta distribution) or manually for various length categories.
#'
#'NEED DETAIL HERE.
#'
#'@aliases simLenSelectP simLenSelectM
#'@param lens A vector containg the lengths of individual fish.
#'@param alpha A numeric shape parameter to the beta distribution.  See
#'\code{\link{dbeta}}.
#'@param beta A numeric shape parameter to the beta distribution.  See
#'\code{\link{dbeta}}.
#'@param max.height A numeric that controls the maximum height of the
#'probability distribution -- i.e., this will be the maximum probability of
#'capture.
#'@param show A logical indicating whether a graphic of the selectivity curve
#'should be shown (\code{=TRUE}) or not (\code{=TRUE}; default).
#'@param breaks A numeric vector of lower values for the break points of the
#'length categories.
#'@param probs A numeric vector of capture probabilities (i.e., selectivities)
#'for each length category.  Default is a vector containing all ones -- i.e.,
#'no selectivity by length category.
#'@param interact A logical indicating whether the capture probabilities (i.e.,
#'selectivities) should be chosen by the user interacting with a selectivity
#'plot.  See details.
#'@param digits A numeric indicating the number of digits that should be used
#'when selecting the capture probabilities.  Smaller values represent coarser
#'choices.
#'@return If \code{simLenSelectP} is used then a vector of logicals indicating
#'whether each fish was sampled (\code{TRUE}) or not.  If \code{simLenSelectM}
#'is used then a list that contains the following three items is returned:
#'\itemize{
#'\item smpld a vector of logicals indicating whether each fish was sampled
#'(\code{TRUE}) or not.
#'\item breaks the vector of length category breaks sent in \code{breaks}.
#'\item probs the vector of capture probabilities that corresponds to the
#'length categories in \code{breaks}.  This vector may not equal the supplied
#'\code{probs} vector if the user changed the capture probabilities with the
#'interactive graphic (i.e., using \code{interact=TRUE}).
#'}
#'@seealso \code{\link{simAges}}, \code{\link{simLenFromAge}}
#'@keywords misc
#'@examples
#'## Load FSA package for Summarize()
#'library(FSA)
#'
#'## set seed for repeatability
#'set.seed(5234734)
#'
#'## Simulated individual ages (random)
#'#    see simAges functions
#'bg.ages <- simAges(N0=500,A=0.35)
#'
#'## Simulated lengths, given the above ages
#'bg.lens <- simLenFromAge(bg.ages,228,0.206,0,8)
#'
#'## Simulated samplings given the above lengths and 
#'##   selectivities from a beta(3,5)
#'bg.smpl <- simLenSelectP(bg.lens,3,5)
#'# append ages, lengths, and samplings into one data frame
#'bg.df1 <- data.frame(age=bg.ages,len=bg.lens,smpld=bg.smpl)
#'# get only those that were sampled
#'bg.df1a <- Subset(bg.df1,smpld)
#'# Summaries
#'Summarize(len~age,data=bg.df1,digits=1)
#'Summarize(len~age,data=bg.df1a,digits=1)
#'
#'## Simulated samplings given the above lengths and user supplied selectivities
#'bg.brks <- seq(20,230,10)
#'bg.prbs1 <- c(0.0,0.0,0.0,0.1,0.3,0.6,0.9,0.9,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.9,
#'  0.8,0.6,0.3,0.2,0.1,0.0)
#'bg.smpl1 <- simLenSelectM(bg.lens,bg.brks,bg.prbs1,interact=FALSE)
#'bg.df2 <- data.frame(age=bg.ages,len=bg.lens,smpld=bg.smpl1$smpld)
#'bg.df2a <- Subset(bg.df2,smpld)
#'Summarize(len~age,data=bg.df2,digits=1)
#'Summarize(len~age,data=bg.df2a,digits=1)
#'
#'## Simulated samplings given the above lengths and selectivities from interactive choices
#'#  NOT RUN because of interactive choices
#'\dontrun{
#'bg.brks <- seq(20,230,10)
#'bg.smpl2 <- simLenSelectM(bg.lens,bg.brks)
#'bg.df3 <- data.frame(age=bg.ages,len=bg.lens,smpld=bg.smpl2$smpld)
#'bg.df3a <- Subset(bg.df3,smpld)
#'Summarize(len~age,data=bg.df3,digits=1)
#'Summarize(len~age,data=bg.df3a,digits=1)
#'}
#'
#'@rdname simLenSelect
#'@export simLenSelectP
simLenSelectP <- function(lens,alpha=1,beta=1,max.height=1,show=FALSE) {
 # some error catching
  if (length(alpha)>1) stop("'alpha' must contain one and only one value",call.=FALSE)
  if (alpha<0) stop("'alpha' must be a postive number",call.=FALSE)
  if (length(beta)>1) stop("'beta' must contain one and only one value",call.=FALSE)
  if (beta<0) stop("'beta' must be a postive number",call.=FALSE)
  if (length(max.height)>1) stop("'max.height' must contain one and only one value",call.=FALSE)
  if (max.height<=0 | max.height>1) stop("'max.height' must be a strictly postive number less than or equal to 1",call.=FALSE)

 # show plot
  if (show) {
    lens.show <- seq(min(lens),max(lens),length.out=200)
    vulns.show <- stats::dbeta(lens.show/max(lens.show),alpha,beta)
    vulns.show <- vulns.show/max(vulns.show)
    graphics::plot(vulns.show~lens.show,type="l",xlab="Length",ylab="Probability of Capture",lwd=2)
  }
 # generate results
 # get densities from beta distribution, with lengths put on 0-1 scale
  vulns <- stats::dbeta(lens/max(lens),alpha,beta)
 # scale vulnerabilities so that maximum vulnerability = max.height
  if (max.height > 1) {
    max.height <- 1
    warning("max.height was greater than 1; it has been rest to 1.",call.=FALSE)
  }
  vulns <- vulns/max(vulns)*max.height
 # generate random number between 0 and 1
  rands <- stats::runif(length(vulns))
 # if random number is less than vulnerability then fish is captured (TRUE)
  rands < vulns
}

#'@rdname simLenSelect
#'@export simLenSelectM
simLenSelectM <- function(lens,breaks,probs=rep(max.height,length(breaks)),max.height=1,interact=TRUE,digits=2) {
  if (interact) {
    old.par <- graphics::par(no.readonly=TRUE) 
    on.exit(graphics::par(old.par))
    options(locatorBell=FALSE)
    graphics::layout(matrix(c(1,2),nrow=2),heights=c(1,15))
    repeat {
     # setup the top panel
      graphics::par(mar=c(0.05,0.4,0.1,0.4),usr=c(0,1,0,1))
      graphics::frame()
      graphics::box(col="red",lwd=2)
      graphics::text(0.5,0.5,"Stop Interactive Choices & Perform Selections",col="red")
     # setup the plot
      graphics::par(mar=c(3.5,3.5,0.1,3.5),mgp=c(2,0.75,0))
      graphics::plot(breaks,probs,xlim=range(breaks),ylim=c(0,max.height),xlab='Length Categories',ylab='Probability of Capture',type="n")
      aty <- graphics::axTicks(2)
      graphics::axis(side=4,at=aty,labels=formatC(aty/aty[length(aty)]*100,format="f",digits=0))
      graphics::mtext("Percentage of Maximum Probability",side=4,line=2)
      graphics::abline(v=breaks,lty=3,lwd=1,col="gray90")
      graphics::abline(h=aty,lty=3,lwd=1,col="gray90")
      graphics::abline(h=c(0,max.height),lty=2,lwd=1,col="red")
      graphics::lines(breaks,probs,col="gray",lwd=2)
      graphics::points(breaks,probs,pch=16)
     # get point
      pnt <- graphics::locator(1)
      if (pnt$y > graphics::par('usr')[4]) { ## clicked in top panel
        break
      } else { ## clicked in bottom panel
       # move the point
        mov.i <- which.min((breaks-pnt$x)^2+(probs-pnt$y)^2)
        if (pnt$y > max.height) pnt$y <- max.height
        if (pnt$y < 0) pnt$y <- 0
        probs[mov.i] <- round(pnt$y,digits)
        next
      }
    } ## end repeat -- end selection of vulnerabilities
  } ## end interact
 # generate results
 # find length categories for each length
  df <- FSA::lencat(stats::as.formula("~lens"),data=data.frame(lens),breaks=breaks,as.fact=FALSE)
 # find probability of capture corresponding to each length category
  vulns <- probs[match(df$LCat,breaks)]
 # generate random number between 0 and 1
  rands <- stats::runif(length(lens))
 # if random number is less than vulnerability then fish is captured (TRUE)
  smpld <- rands < vulns
 # return results
  list(smpld=smpld,breaks=breaks,probs=probs)
}
