# Load required packages
library(coda)
library(bdskytools)
library(beastio)
library(RColorBrewer)

setwd("C:/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Re_data")

# Define colors
cols <- list(
  blue = RColorBrewer::brewer.pal(12, "Paired")[2],
  orange = RColorBrewer::brewer.pal(12, "Paired")[8]
)

# Load trace files
logfile1 <- "C:/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/VA94_consensus_lineage1_trimmed_60threshold_50.log"
logfile2 <- "C:/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/VA94_consensus_lineage2_trimmed_60threshold_50.log"

bdsky_trace1 <- beastio::readLog(logfile1, burnin=0.1)
bdsky_trace2 <- beastio::readLog(logfile2, burnin=0.1)

# Extract Re values
Re_sky1 <- beastio::getLogFileSubset(bdsky_trace1, "BDSKY_Serial")
Re_sky2 <- beastio::getLogFileSubset(bdsky_trace2, "BDSKY_Serial")

# Compute TMRCA
tmrca_med1 <- median(bdsky_trace1[, "Tree.height"])
tmrca_med2 <- median(bdsky_trace2[, "Tree.height"])

# Define grid times separately for each lineage
gridTimes1 <- seq(0, tmrca_med1, length.out=59)
gridTimes2 <- seq(0, tmrca_med2, length.out=59)

# Interpolate both lineages separately
Re_gridded1 <- mcmc(bdskytools::gridSkyline(Re_sky1, bdsky_trace1[, "origin_BDSKY_Serial"], gridTimes1))
Re_gridded_hpd1 <- t(getHPDMedian(Re_gridded1))

Re_gridded2 <- mcmc(bdskytools::gridSkyline(Re_sky2, bdsky_trace2[, "origin_BDSKY_Serial"], gridTimes2))
Re_gridded_hpd2 <- t(getHPDMedian(Re_gridded2))

# Convert time to actual years (adjust based on reference year)
times1 <- 2004 - gridTimes1
times2 <- 2015 - gridTimes2

# Determine common y-axis range with 0 as the lower limit
ylim_range <- range(c(0, Re_gridded_hpd1, Re_gridded_hpd2), na.rm=TRUE)
ylim_range[1] <- 0

# Determine common x-axis range
xlim_range <- range(c(times1, times2))

# Save variables to an RData file for later use in Python
save(Re_gridded_hpd1, Re_gridded_hpd2, times1, times2, file="lineages_data.RData")

# Save the plot as a high-resolution PNG file
png("combined_lineages_plot.png", width=6000, height=3600, units="px", res=600)
# Plot first lineage
plotSkylinePretty(
  times1, Re_gridded_hpd1, type='smooth', axispadding=0.0,
  col=cols$blue, fill=adjustcolor(cols$blue, alpha.f=0.5),
  xlab="Date", ylab=expression("R"[e]),
  side=2, yline=2.5, xline=2, xgrid=TRUE, ygrid=TRUE,
  ylim=ylim_range, xlim=xlim_range
)

# Overlay second lineage with different x-values but same y-axis scale
par(new=TRUE)
plotSkylinePretty(
  times2, Re_gridded_hpd2, type='smooth', axispadding=0.0,
  col=cols$orange, fill=adjustcolor(cols$orange, alpha.f=0.5),
  xlab="", ylab="", side=2, yline=2.5, xline=2, xgrid=FALSE, ygrid=FALSE,
  ylim=ylim_range, xlim=xlim_range
)

# Add legend
legend("topright", legend=c("Lineage 1", "Lineage 2"),
       col=c(cols$blue, cols$orange), lty=1, bty='n')

# Close the PNG device
dev.off()
