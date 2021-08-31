model <- function(){
	sigma ~ dlnorm(-1, pow(1, -2))
	beta0 ~ dnorm(0, pow(3, -2))
	for(j in 1:nHab){
		beta1[j] ~ dnorm(0, pow(3, -2))
	}
#-----------------------------------------------------------------------------
	# Link to habitat data
	for(i in 1:nPop){
		trendPred[i] <- beta0 + t(beta1[1:nHab]) %*% habPressures[i, 1:nHab]
	}

		for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2))
	}
} # end model

###############################################################################
# Simulate data
###############################################################################

set.seed(238)
beta0 <- -0.3
beta1 <- rnorm(10, mean = -0.1, sd = 1)
sigma <- 0.001

JAGSdat <- list(
	nPop = dim(popDat)[1],
	nHab = length(habNames),
	habPressures = habPressures
)

JAGSdat$trend <- rnorm(JAGSdat$nPop, mean = beta0 + t(as.matrix(beta1)) %*% t(JAGSdat$habPressures), sd = sigma)

plot(JAGSdat$habPressures[, 4], JAGSdat$trend)

###############################################################################
# Fitting
###############################################################################

# Parameters to track
parNames <- c(
	"beta0", # length nSpecies
	"beta1", # length nHab
	"sigma" # length 1
)


numChains <- 4
numIter <- 3000

cl <- makeCluster(numChains)

t.start<-proc.time()[3]

fit <- jags.parfit(
	cl,
	data = JAGSdat,
	params = parNames,
	model = model,
	n.iter = numIter,
	n.adapt = 5000,
	n.chains = numChains)

stopCluster(cl)
procTime <- round((proc.time()[3] - t.start)/(60), 1) # track process time
cat(paste("Process time for model =", procTime, "minutes"))

###############################################################################
# Looking at output
###############################################################################

S <- summary(fit)
parNames <- rownames(S[[1]])
out <- array(NA, dim = c(length(fit), dim(fit[[1]])))
for(i in 1:12){
	for(j in 1:numChains){
		out[j, ,i] <- fit[[j]][, i]
	}
}

par(mfrow = c(4,3), mar = c(1,1,2,1), oma = c(2,2,0,0))
for(i in 1:12){
	plot(1:numIter, out[1, , i], "n", ylim = quantile(out[, , i], c(0.025, 0.975)), xlab = "", ylab = "", main = parNames[i])
	for(j in 1:numChains){
		lines(1:numIter, out[j, , i], lty = 2, col = j)
	}
	
	abline(h = c(beta0, beta1, sigma)[i], lwd = 2)
}

# Simple model does pretty damn well