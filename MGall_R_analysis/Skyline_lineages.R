# Load the required packages and set global options
library(coda)
library(bdskytools)
library(beastio)
library(RColorBrewer)

setwd("C:/Users/at991/Software/BEAST.v2.7.7.Windows/BEAST/bat")

# Set up some colours
cols  <- list(blue   = RColorBrewer::brewer.pal(12,"Paired")[2],
              orange = RColorBrewer::brewer.pal(12,"Paired")[8])

set_alpha <- function(c, alpha=1.0) paste0(c,format(as.hexmode(round(alpha*255)), width=2))

# Load the trace files and check convergence
logfile1 <- "C:/Users/at991/Software/BEAST.v2.7.7.Windows/BEAST/bat/VA94_all_60threshold_lineage1.log"
logfile2 <- "C:/Users/at991/Software/BEAST.v2.7.7.Windows/BEAST/bat/VA94_all_60threshold_lineage2.log"

bdsky_trace1 <- beastio::readLog(logfile1, burnin=0.1)
bdsky_trace2 <- beastio::readLog(logfile2, burnin=0.1)

# Extract parameter estimates and HPDs for each lineage
Re_sky1 <- beastio::getLogFileSubset(bdsky_trace1, "BDSKY_Serial")
Re_sky2 <- beastio::getLogFileSubset(bdsky_trace2, "BDSKY_Serial")

Re_hpd1 <- t(beastio::getHPDMedian(Re_sky1))
Re_hpd2 <- t(beastio::getHPDMedian(Re_sky2))

# Plotting a "smooth" skyline

tmrca_med1  <- median(bdsky_trace1[, "Tree.height"])
gridTimes1  <- seq(0, median(tmrca_med1), length.out=100)

tmrca_med2  <- median(bdsky_trace2[, "Tree.height"])
gridTimes2  <- seq(0, median(tmrca_med2), length.out=100)

Re_gridded1 <- mcmc(bdskytools::gridSkyline(Re_sky1, bdsky_trace1[, "origin_BDSKY_Serial"], gridTimes1))
Re_gridded_hpd1 <- t(getHPDMedian(Re_gridded1))

Re_gridded2 <- mcmc(bdskytools::gridSkyline(Re_sky2, bdsky_trace2[, "origin_BDSKY_Serial"], gridTimes2))
Re_gridded_hpd2 <- t(getHPDMedian(Re_gridded2))

# Plotting combined Re values
times1 <- 2004 - gridTimes1
times2 <- 2015 - gridTimes2

# Combine the plots into one with log scale for y-axis
plotSkyline(times1, Re_gridded_hpd1, xlab="Date", ylab="Re", type="smooth", col=cols$blue, main="Combined Re values for two lineages", log="y")
plotSkyline(times2, Re_gridded_hpd2, xlab="Date", ylab="Re", type="smooth", col=cols$orange, add=TRUE, log="y")

legend("topright", legend=c("Lineage 1", "Lineage 2"), col=c(cols$blue, cols$orange), lty=1, bty='n')