# Install the required packages by evaluating this chunk
# (only needs to be evaluated the first time)

# install.packages("devtools")
# install.packages("coda")
# install.packages("RColorBrewer")
# devtools::install_github("laduplessis/bdskytools")
# devtools::install_github("laduplessis/beastio")

# Load the required packages and set global options
library(coda)
library(bdskytools)
library(beastio)
library(RColorBrewer)

setwd("C:/Users/at991/Software/BEAST.v2.7.7.Windows/BEAST/bat")
knitr::opts_chunk$set(echo = TRUE, fig.path = "figs/",
                      dev = "png", fig.width = 7, fig.height = 5)

# Set up some colours
cols  <- list(blue   = RColorBrewer::brewer.pal(12,"Paired")[2],
              orange = RColorBrewer::brewer.pal(12,"Paired")[8])

set_alpha <- function(c, alpha=1.0) paste0(c,format(as.hexmode(round(alpha*255)), width=2))



# Load the trace file and check convergence

bdsky_trace   <- beastio::readLog("C:/Users/at991/Software/BEAST.v2.7.7.Windows/BEAST/bat/VA94_all_60_threshold_birthdeath.log", burnin=0.1)

summary(bdsky_trace)
varnames(bdsky_trace)


# We can use the `checkESS()` function to find which parameters have ESS < 200,

beastio::checkESS(bdsky_trace)


# or use the same function to plot the ESS values of all parameters.


beastio::checkESS(bdsky_trace,   cutoff=200, plot=TRUE, log='y', ylim=c(1,10000), title="All parameters", plot.grid=TRUE)


# Extract parameter estimates and HPDs

# Next we can extract the $R_e$ parameter values and their HPDs. 

Re_sky <- beastio::getLogFileSubset(bdsky_trace, "BDSKY_Serial")
Re_hpd <- t(beastio::getHPDMedian(Re_sky))
delta_hpd <- beastio::getHPDMedian(bdsky_trace[, "becomeUninfectiousRate_BDSKY_Serial"])


# Plotting non-gridded BDSKY estimates
#We can plot the raw $R_e$ HPD intervals. This is equivalent to the output in Tracer.


bdskytools::plotSkyline(1:14, Re_hpd, type='step', ylab="Re")


# Plotting a "smooth" skyline

tmrca_med  <- median(bdsky_trace[, "Tree.height"])
gridTimes  <- seq(0, median(tmrca_med), length.out=100)  

Re_gridded <- mcmc(bdskytools::gridSkyline(Re_sky, bdsky_trace[, "origin_BDSKY_Serial"], gridTimes))
Re_gridded_hpd <- t(getHPDMedian(Re_gridded))

times <- 2015 - gridTimes
plotSkyline(times, Re_gridded_hpd, xlab="Date", ylab="Re", type="smooth")

# Plotting multiple samples
layout(matrix(c(1:4), nrow=4))
plotSkyline(times, Re_gridded, type='steplines', traces=1, 
            col=cols$blue, ylims=c(0,3.5), xlab="Time", ylab="Re", main="1 random sample")
plotSkyline(times, Re_gridded, type='steplines', traces=10, 
            col=set_alpha(cols$blue,0.5), ylims=c(0,3.5), xlab="Time", ylab="Re", main="10 random samples")
plotSkyline(times, Re_gridded, type='steplines', traces=100, 
            col=set_alpha(cols$blue,0.5), ylims=c(0,3.5), xlab="Time", ylab="Re", main="100 random samples")
plotSkyline(times, Re_gridded, type='steplines', traces=1000, 
            col=set_alpha(cols$blue,0.1), ylims=c(0,3.5), xlab="Time", ylab="Re", main="1000 random samples")

# Combined plot
par(mar=c(4,4,1,4))
plot(1, type='n', xlim=c(1990,2015), ylim=c(0,1), 
     xlab='Year', ylab="", yaxt='n', xaxs='i', yaxs='i')

plotSkyline(range(times), as.matrix(delta_hpd), type='step', lwd=2, 
            xlab="", ylab="", add=TRUE, new=FALSE, axes=FALSE, fill=set_alpha(cols$blue, 0.5), col=cols$blue)
axis(4, las=1)
mtext(expression(delta), side=4, line=3)

plotSkyline(times, Re_gridded_hpd, lwd=2, xlims=c(1990,2015), xaxs='i', yaxs='i', 
            xlab="", ylab="", add=TRUE, new=TRUE, axes=FALSE, fill=set_alpha(cols$orange, 0.5), col=cols$orange)   
axis(2, las=1)
mtext(expression("R"[e]), side=2, line=3)

abline(h = 1, lty=2, col=cols$red, lwd=2)

legend("topright", legend=c(expression("R"[e]), expression(delta)), bty='n',
       fill=set_alpha(c(cols$orange, cols$blue), 0.5), border=c(cols$orange, cols$blue))