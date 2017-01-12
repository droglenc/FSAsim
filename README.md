FSAsim
======

This is the github page for the FSAsim package for R developed by [Derek Ogle](http://derekogle.com/) at [Northland College](http://www.northland.edu/).  This is very much a *work-in-progress* so please use at your discretion.

## Installation
This package can be installed from github to your R with the following code

```r
if (!require('devtools')) install.packages('devtools'); require('devtools')
devtools::install_github('droglenc/FSAsim')
```

Descriptions of recent changes can be found in the [News.md file](https://github.com/droglenc/FSAsim/blob/master/NEWS.md)

## Note About Using Macs
**FSAsim** uses **TCL/TK** for some interactive plots.  Some Mac users report problems with using **TCL/TK**.  I do not have access to a Mac to test these problems. However, the CRAN page suggests that for recent versions of R (>3.0.0), [XQuartz](https://www.xquartz.org/) must be installed. In the past, some students have reported success installing the **TCL/TK** universal build [located here](http://cran.r-project.org/bin/macosx/tools/) (or [direct link to the file](http://cran.r-project.org/bin/macosx/tools/tcltk-8.5.5-x11.dmg)).

## Contact
Contact me with questions by sending a friendly e-mail to <dogle@northland.edu>.
