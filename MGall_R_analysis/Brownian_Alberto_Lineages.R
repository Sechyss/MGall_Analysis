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

	n <- length(x)

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

	# Plot lineage-specific virulence vs dates
	graphics.off()
	windows(width=7, height=5)
	x_plot <- x[match(names(dates), names(x))]
	plot(dates, x_plot, pch=19, col='#00808080', main=paste('Lineage', label))
	cf <- coef(lm(x_plot ~ dates))
	lines(range(dates), cf[1] + range(dates) * cf[2])
	mtext(side=3, substitute(paste(ZZ,': (',hat(lambda),' = ',XX,'; ',italic(p),'=',YY,')'),
		  list(ZZ=paste('Lineage', label, 'virulence index'),
			   XX=round(res$lam, 3), YY=round(res$p, 4))))

	invisible(res)
}

# Run analysis per lineage
res_L1 <- analyze_lineage(l1, 'L1')
res_L2 <- analyze_lineage(l2, 'L2')

cat('\nSummary:\n')
if (!is.null(res_L1)) cat('L1: lambda=', round(res_L1$lam,3), '; p=', round(res_L1$p,4), '; anc=', round(res_L1$anc,3), '\n', sep='')
if (!is.null(res_L2)) cat('L2: lambda=', round(res_L2$lam,3), '; p=', round(res_L2$p,4), '; anc=', round(res_L2$anc,3), '\n', sep='')
