# Testing Strategy for Curve Fitting

## Negative Binomial
The negative binomial distribution fitting algorithm was tested as follows.

We ran the following r script, which produced a series of negative binomial distributed
random variates for known input parameters.

library(gridExtra)
library(ggplot2)

numVariates <- 10000
variateList <- list()
plotList <- list()

for (i in 1:5) {
	print(i)
	r <- 2
	p <- i/20
	print(p)
	
	currVariates <- rnbinom(numVariates, size=r, prob = p)
	
	filename <- paste("negativeBinomialVariates.r.", r, ".p.", p, ".txt", sep="")
	writeLines(as.character(currVariates), con=filename)
	
	currHist <- ggplot(data.frame(x=currVariates), aes(x)) + geom_histogram(aes(y=..density..), bindwidth=1, fill="lightblue", color="black") + stat_function(fun=dnbinom, args = list(size = r, prob = p), color="red", linetype="dashed") + labs(title=paste("r=", r, ", p=",p), x="Values", y=" Frequency") + theme_classic() + xlim(0, 100)
	plotList[[i]] <- currHist
}

grid.arrange(grobs = plotList, ncol=2)

With these variates in hand, we then fed the variants into fit.NegativeBinomial and evaluated whether the sample estimated
parameters fell within 5% of the true parameters.


## Poisson

Poisson fitting testing followed a similar script as in Negative Binomial testing.

library(gridExtra)
library(ggplot2)

numVariates <- 10000
variateList <- list()
plotList <- list()

for (i in 1:5) {
	print(i)
	lambda <- i
	
	currVariates <- rpois(numVariates, lambda = lambda)
	
	filename <- paste("poissonVariates.lambda.", lambda, ".txt", sep="")
	writeLines(as.character(currVariates), con=filename)
	
	currHist <- ggplot(data.frame(x=currVariates), aes(x)) + geom_histogram(aes(y=..density..), bindwidth=1, fill="lightblue", color="black") + labs(title=paste("lambda=", lambda), x="Values", y=" Frequency") + theme_classic() + xlim(0, 15)
	plotList[[i]] <- currHist
}

grid.arrange(grobs = plotList, ncol=2)

