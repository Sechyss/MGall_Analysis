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

# Get the data together
#x <- logit(at12); names(x) <- rownames(core); xtit <- 'AT positions 1 and 2'
#x <- logit(at3); names(x) <- rownames(core); xtit <- 'AT 3rd positions'
#x <- sqrt(t$peak_score); names(x) <- t$Folder.name; xtit <- 'peak score'
#x <- sqrt(t$mean_swel); names(x) <- t$Folder.name; xtit <- 'mean swell'
x <- sqrt(t$peak_score)+sqrt(t$mean_swel); names(x) <- t$final.name; xtit <- 'JW virulence index'

mx <- mean(x,na.rm=T); sdx <- sd(x,na.rm=T)
x <- (x-mx)/sdx

# Use the temproal signal to root the tree
dates <- t$year.sampling; names(dates) <- t$final.name
dates <- dates[match(tre$tip.label,names(dates))]
tre <- rtt(tre,dates)
tre$edge.length <- tre$edge.length/max(tre$edge.length)

graphics.off()
layout(matrix(1:2,2,1))
plot(dates,x,pch=19,col='#00808080')

#plot(dates,meanswel,pch=19,col='#00808080')
#stop()

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
res <- doopt(x,tre)

graphics.off()
windows(width=7,height=5)
x<- x[match(names(dates),names(x))]
plot(dates,x,pch=19,col='#00808080')
cf <- coef(lm(x ~ dates))
lines(range(dates),cf[1]+range(dates)*cf[2])
#lines(range(dates),res$anc+range(dates)*res$lam,col='red')
mtext(side=3,substitute(paste(ZZ,': (',hat(lambda),' = ',XX,'; ',italic(p),'=',YY,')'),list(ZZ=xtit,XX=round(res$lam,3),YY=round(res$p,4))))
