##https://www.statlect.com/fundamentals-of-statistics/multivariate-normal-distribution-maximum-likelihood
##xsim <- mvrnorm(n = 1, mu=rep(mu.sim,n), Sigma=A*sig.sim, tol = 1e-6, empirical = FALSE)

library(ape)
#library(adephylo)
#library(matrixcalc)
library(MASS)

setwd('C:/Users/at991/PycharmProjects/MGall_Analysis/MGall_R_analysis/Brownian_Alberto')

logit <- function(z) log(z/(1-z)) 
#load('at_new.Rdata')
t <- read.table('meta.txt',sep='\t',header=T)
tre<-read.tree("Edited_VA94_consensus_all_trimmed_60threshold_50_combined.finaltree.newick")

tre$edge.length=tre$edge.length+(1/500000)

dop <- function(r)
{
	lam <- 0
	if(length(r$par)==3) lam <- r$par[3]
	cat('Mu = ',r$par[1],'; Sigma = ',r$par[2]^2,'; lambda = ',lam,' (lnL=',r$value,'; conv=',r$convergence,')\n',sep='')  
}
doopt <- function(x,tre)
{
	tre <- drop.tip(tre,names(x)[which(is.na(x))])
	
	A <- vcv(tre)
	invA <- solve(A)
	rtt <- diag(A)
	n <- dim(invA)[1]
	print(n)
	x <- x[match(colnames(A),names(x))]
			
	if(!all(rownames(invA)==names(x)) || !all(rownames(invA)==names(rtt)))
		stop('Problem with names')
	
	cat('Null\n')
	r <- optim(fn=lnL,par=runif(2),control=list(fnscale=-1),x=x,invA=invA,n=n,rtt=rtt)
	dop(r)
	lnLNull <- r$value
	r <- optim(fn=lnL,par=runif(2),control=list(fnscale=-1),x=x,invA=invA,n=n,rtt=rtt,method='SANN')
	dop(r)
	cat('Full\n')
	r <- optim(fn=lnL,par=runif(3),control=list(fnscale=-1),x=x,invA=invA,n=n,rtt=rtt)
	lnLFull <- r$value
	lambda <- r$par[3]
	anc.state <- r$par[1]
	dop(r)
	r <- optim(fn=lnL,par=runif(3),control=list(fnscale=-1),x=x,invA=invA,n=n,rtt=rtt,method='SANN')
	dop(r)
	pval <- 1-pchisq(2*(lnLFull-lnLNull),df=1)
	cat('\np=',pval,'\n')
	
	res <- list(p=pval,lam=lambda,anc=anc.state)
}

## Lineage-specific analysis will be performed below

#### The loglikelihood function... ####
lnL <- function(params,x,invA,n,rtt)
{
	mu <- params[1]
	sigma <- params[2]^2

	if(length(params)==2)
		lambda<- 0
	else
		lambda <- params[3]
	
	# Vector of means
	m <- mu + rtt*lambda
	v <- x-m
	lnL <- -(t(v) %*% (invA/sigma) %*% v)/2   
	lnL <- lnL -log(sigma)*n/2
	lnL
}


# Read lineage taxa lists
l1 <- trimws(readLines('leaves_L1.txt'))
l2 <- trimws(readLines('leaves_L2.txt'))

analyze_lineage <- function(leaves, label) {
	taxa <- intersect(leaves, tre$tip.label)
	taxa <- intersect(taxa, t$final.name)
	if (length(taxa) < 3) {
		cat('Skipping ', label, ' - insufficient taxa after intersection (', length(taxa), ')\n', sep='')
		return(NULL)
	}

	t_sub <- t[t$final.name %in% taxa, ]

	# Virulence index within lineage
	x <- sqrt(t_sub$peak_score) + sqrt(t_sub$mean_swel)
	names(x) <- t_sub$final.name
	mx <- mean(x, na.rm=TRUE); sdx <- sd(x, na.rm=TRUE)
	x <- (x - mx) / sdx

	# Dates for temporal rooting
	dates <- t_sub$year.sampling; names(dates) <- t_sub$final.name

	# Subset tree to lineage taxa
	tre_sub <- drop.tip(tre, setdiff(tre$tip.label, taxa))
	dvec <- dates[match(tre_sub$tip.label, names(dates))]

	# Root and scale by temporal signal
	tre_sub <- rtt(tre_sub, dvec)
	tre_sub$edge.length <- tre_sub$edge.length / max(tre_sub$edge.length)

	# Run optimization
	res <- doopt(x, tre_sub)

	# Return data for combined plotting
	list(res=res, x=x, dates=dates, label=label)
}

# Run analysis per lineage
res_L1 <- analyze_lineage(l1, 'L1')
res_L2 <- analyze_lineage(l2, 'L2')

# Combined plot for both lineages
if (!is.null(res_L1) || !is.null(res_L2)) {
	graphics.off()
	windows(width=8, height=6)

	# Prepare data
	d1 <- if (!is.null(res_L1)) res_L1$dates else NULL
	d2 <- if (!is.null(res_L2)) res_L2$dates else NULL
	x1 <- if (!is.null(res_L1)) res_L1$x[match(names(d1), names(res_L1$x))] else NULL
	x2 <- if (!is.null(res_L2)) res_L2$x[match(names(d2), names(res_L2$x))] else NULL

	d_all <- c(d1, d2)
	x_all <- c(x1, x2)

	# Create empty plot with combined ranges
	plot(range(d_all, na.rm=TRUE), range(x_all, na.rm=TRUE), type='n',
		 xlab='Sampling year', ylab='Normalized virulence',
		 main='Virulence vs sampling year (L1 blue, L2 orange)')

	# Plot points
	if (!is.null(d1)) points(d1, x1, pch=19, col='blue')
	if (!is.null(d2)) points(d2, x2, pch=19, col='orange')

	# Regression lines per lineage
	if (!is.null(d1)) {
		cf1 <- coef(lm(x1 ~ d1))
		lines(range(d1, na.rm=TRUE), cf1[1] + range(d1, na.rm=TRUE) * cf1[2], col='blue', lwd=2)
	}
	if (!is.null(d2)) {
		cf2 <- coef(lm(x2 ~ d2))
		lines(range(d2, na.rm=TRUE), cf2[1] + range(d2, na.rm=TRUE) * cf2[2], col='orange', lwd=2)
	}

	legend('topleft', legend=c('Lineage 1','Lineage 2'), col=c('blue','orange'),
		   pch=19, lwd=2, bty='n')

	# Summary text
	if (!is.null(res_L1) && !is.null(res_L2)) {
		mtext(side=3, paste('L1: lambda=', round(res_L1$res$lam,3), ', p=', round(res_L1$res$p,4),
							' | L2: lambda=', round(res_L2$res$lam,3), ', p=', round(res_L2$res$p,4), sep=''))
	} else if (!is.null(res_L1)) {
		mtext(side=3, paste('L1: lambda=', round(res_L1$res$lam,3), ', p=', round(res_L1$res$p,4), sep=''))
	} else if (!is.null(res_L2)) {
		mtext(side=3, paste('L2: lambda=', round(res_L2$res$lam,3), ', p=', round(res_L2$res$p,4), sep=''))
	}
}

# Console summary
cat('\nSummary:\n')
if (!is.null(res_L1)) cat('L1: lambda=', round(res_L1$res$lam,3), '; p=', round(res_L1$res$p,4), '\n', sep='')
if (!is.null(res_L2)) cat('L2: lambda=', round(res_L2$res$lam,3), '; p=', round(res_L2$res$p,4), '\n', sep='')
