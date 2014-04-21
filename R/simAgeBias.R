#'Simulate an age bias and simulate the ages of fish using the age bias.
#'
#'Constructs an age-bias table interactively.  This age-bias table can then be
#'used to convert the true ages of a sample of fish to \sQuote{biased} ages.
#'
#'NEED DETAIL HERE.
#'
#'@aliases simAgeBias simApplyAgeBias
#'@param max.age A numeric indicating the maximum age to be modelled in the age
#'bias table.
#'@param show.props A logical indicating whether proportions of fish at x-axis
#'age or numbers of fish should be shown on the interactive plot (\code{=TRUE};
#'default) or not.
#'@param scale A logical indicating whether the plotted value of the
#'proportions or numbers should be scaled in relation to their numeric value
#'(\code{=TRUE}; default) or not.
#'@param ages A vector containg the \sQuote{true} ages of individual fish.
#'@param bias.table A table that contains, as columns, the proportion of fish
#'of a certain \sQuote{true} age in various \sQuote{biased} ages -- i.e., a
#'column-proportions table constructed from an age agreement table where the
#'\sQuote{true} ages correspond to columns.
#'@param agree.table A table that contains the age-agreement table where the
#'\sQuote{true} ages correspond to columns.
#'@return If \code{simApplyAgeBias} is used then a vector of \sQuote{biased}
#'ages is returned.  If \code{simAgeBias} is used then a list with the
#'following two items is returned:
#'\itemize{
#'\item agree a table containing the age agreement table resulting from the
#'interactive process.
#'\item bias a table containing the bias table resulting from the interactive
#'process.
#'}
#'@seealso \code{\link{simAges}}, \code{ageComp} in \pkg{FSA}
#'@keywords misc
#'@examples
#'## set seed for repeatability
#'set.seed(5234734)
#'
#'## Simulated individual ages (random)
#'#    see simAges functions
#'bg.ages <- simAges(N0=500,A=0.35)
#'summary(bg.ages)
#'
#'## Simulated ages given the above 'true' ages and age biases from interactive choices
#'#  NOT RUN because of interactive choices
#'\dontrun{
#'bg.ab <- simAgeBias(max.age=max(bg.ages))
#'bg.ages2 <- simApplyAgeBias(bg.ages,bg.ab$bias)
#'summary(bg.ages2)
#'}
#'
#'@rdname simAgeBias
#'@export simAgeBias
simAgeBias <- function(max.age=10,show.props=TRUE,scale=TRUE) {
  Freq <- NULL  # attempting to get by bindings warning in RCMD CHECK
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  options(locatorBell=FALSE)
  x <- y <- NULL
  mode <- "add"
  layout(matrix(c(2,1),nrow=1),widths=c(5,1))
  repeat { ## right panel
    par(mar=c(3.5,0,1,0.5),usr=c(0,1,0,1))
    frame()
    box()
    text(rep(0.5,4),c(0.8,0.625,0.425,0.225),lab=c("Stop\nInteraction","Add","Delete","Move") )
    lines(c(0.05,0.05,0.95,0.95,0.05),c(0.75,0.85,0.85,0.75,0.75)) 
    points(rep(0.5,3),c(0.575,0.375,0.175),
           pch=c(ifelse(mode=="add",16,1),
                 ifelse(mode=="del",16,1),
                 ifelse(mode=="mov",16,1)),cex=2.5 )
    ## left panel
    par(mar=c(3.5,3.5,1,1),mgp=c(2,0.75,0))
    plot(0,0,type="n",xlim=c(0,max.age),ylim=c(0,max.age),xlab="True Age",ylab="Biased Age")
    abline(a=0,b=1,lwd=1,col="gray50")
    abline(h=0:max.age,lty=3,lwd=1,col="gray90")
    abline(v=0:max.age,lty=3,lwd=1,col="gray90")  
    if (length(x)>0) {
      vals <- table(y,x)
      props <- prop.table(vals,margin=2)
      vals <- data.frame(vals)
      props <- data.frame(props)
      vals <- Subset(vals,Freq>0)
      props <- Subset(props,Freq>0)
      ifelse(scale,cxs <- 0.5*props$Freq+0.5,cxs <- 1)   # rescale props to be between 0.5 and 1 for plotting size
      if (show.props) {
        with(props,text(fact2num(x),fact2num(y),formatC(Freq,format="f",digits=2),cex=cxs))
      } else {
        with(vals,text(fact2num(x),fact2num(y),formatC(Freq,format="f",digits=0),cex=cxs))
      }
    }
    ns <- table(factor(x,levels=0:max.age))
    text(0:max.age,rep(0,max.age+1),ns,col="blue")
    # get point
    pnt <- locator(1)
    if (pnt$x > par('usr')[2]) { ## clicked in left panel
      pnt2 <- cnvrt.coords(pnt)$fig
      if (pnt2$y>0.7) { break }
      if (pnt2$y>0.5) { mode <- "add"
                        next  }
      if (pnt2$y>0.3) { mode <- "del"
                        next  }
      mode <- "mov"
      next
    } else { ## clicked in right panel
      if (mode=="add") {
        x <- c(x,round(pnt$x,0))
        y <- c(y,round(pnt$y,0))
        next
      }
      if (mode=="del") {
        min.i <- which.min((x-pnt$x)^2+(y-pnt$y)^2)
        x <- x[-min.i]
        y <- y[-min.i]
        next
      }
      if(mode=="mov") {
        mov.i <- which.min((x-pnt$x)^2+(y-pnt$y)^2)
        points(x[mov.i],y[mov.i],pch=16)
        pnt <- locator(1)
        x[mov.i] <- round(pnt$x,0)
        y[mov.i] <- round(pnt$y,0)
        next
      }
    }
  } ## end repeat
  if (is.null(x)) stop("No points were added to the plot.  Nothing is returned.",call.=FALSE)
  else {
    df <- data.frame(true=x,bias=y)
    ac <- ageComp(x,y,col.lab="Truth",row.lab="Biased")
    agree.raw <- ac$agree
    agree.prop <- prop.table(agree.raw,margin=2)
    list(agree=agree.raw,bias=agree.prop)
  }
}

#'@rdname simAgeBias
#'@export simApplyAgeBias
simApplyAgeBias <- function(ages,bias.table=NULL,agree.table=NULL) {
 # some checking
  if (length(ages) == 0) stop("'ages' must contain data.",call.=FALSE)
  if (is.null(bias.table) & is.null(agree.table)) {
    stop("One of 'bias.table' or 'agree.table' must be provided",call.=FALSE)
  }
  if (!is.null(bias.table) & !is.null(agree.table)) {
    warning("Both 'bias.table' and 'agree.table were provided.\n Only 'bias.table' will be used",call.=FALSE)
    agree.table <- NULL
  }
 # if only an agreement table given then conver to bias table 
  if (is.null(bias.table) & !is.null(agree.table)) {
    bias.table <- prop.table(agree.table,margin=2)
  }
 # find the ages present in the bias table
  ages.in.table <- fact2num(colnames(bias.table))
 # some more checking
  if (max(ages) > max(ages.in.table)) {
    stop("An observed age is greater than the maximum age in the bias table.",call.=FALSE)
  }
  if (min(ages) < min(ages.in.table)) {
    stop("An observed age is less than the minimum age in the bias table.",call.=FALSE)
  }
 # initiate the vector to hold the derived 'biased' ages 
  b.ages <- numeric(length(ages))
 # create the biased ages -- one individual at a time
  for (i in 1:length(ages)) {
   # find the column in the bias table that corresponds to the obseved true age
    col.ind <- match(ages[i],ages.in.table)
   # sample out of all ages with probabilities given by column in bias table 
    b.ages[i] <- sample(ages.in.table,1,prob=bias.table[,col.ind])
  }
 # return vector of derived 'biased' ages 
  b.ages
}
