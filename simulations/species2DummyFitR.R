model <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	# Different for 6 marine ecotypes: 
	for(s in 1:nSpecies){
		beta0[s] ~ dnorm(0, pow(3, -2))
		for(j in 1:nHab){
			# Strength of popHab relationship
			beta1[s, j] ~ dnorm(0, pow(3, -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		trendPred[i] <- beta0[species[i]] + inprod(beta1[species[i], ], t(habPressures[i, ]))
		}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2))
	}
} # end model

###############################################################################
# Simulate data
###############################################################################

set.seed(238)
beta0 <- rnorm(5, mean = -0.1, sd = 0.1)
beta1 <- matrix(rnorm(10*5, mean = 0, sd = 0.03), nrow = 5, ncol = 10)
sigma <- 0.1

JAGSdat <- list(
	nPop = dim(popDat)[1],
	nFWeco = length(FWecoNames),
	nHab = length(habNames),
	nSpecies = length(speciesNames),
	FWeco = as.numeric(factor(popDat$FWecotype2, levels = FWecoNames)),
	species = as.numeric(factor(popDat$SPECIES, levels = speciesNames)),
	habPressures = habPressures
)

# Problem with how slope is calculated here in simulated data
# Need to sort out in model fitting too...

JAGSdat$trend <- rnorm(JAGSdat$nPop, mean = beta0[JAGSdat$species] + as.matrix(beta1[JAGSdat$species, 1:10]) %*% t(JAGSdat$habPressures), sd = sigma)

range(JAGSdat$trend)

plot(JAGSdat$habPressures[, 3], JAGSdat$trend, col = JAGSdat$species)

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
for(i in 1:dim(S[[1]])[1]){
	for(j in 1:numChains){
		out[j, ,i] <- fit[[j]][, i]
	}
}

gelman.diag(fit)
# Trace plots

# Beta0
par(mfrow = c(3,2), mar = c(1,1,2,1), oma = c(2,2,0,0))
for(i in grep("beta0", parNames)){
	plot(1:numIter, out[1, , i], "n", ylim = range(c(beta0[i], out[, , i])), xlab = "", ylab = "", main = parNames[i])
	for(j in 1:numChains){
		lines(1:numIter, out[j, , i], lty = 2, col = j)
	}
	
	abline(h = c(beta0)[i], lwd = 2)
}

# Beta1
for(s in 1:5){
	I <- grep(paste(s, ",", sep = ""), parNames)
	par(mfrow = c(4,3), mar = c(1,1,2,1), oma = c(2,2,0,0))
	for(i in 1:10){
		plot(1:numIter, out[1, , I[i]], "n", ylim = range(c(beta1[s, i], out[, , I[i]])), xlab = "", ylab = "", main = parNames[I[i]])
		for(j in 1:numChains){
			lines(1:numIter, out[j, , I[i]], lty = 2, col = j)
		}
	
	abline(h = beta1[s,i], lwd = 2)
}
}


plot(JAGSdat$habPressures[JAGSdat$species == 1, 1], JAGSdat$trend[JAGSdat$species == 1])
abline(a = S[[1]]['beta0[1]', 1], b = S[[1]]['beta1[1,1]', 1], col = 2)
abline(a = beta0[1], b = beta1[1,1], col = 4, lwd = 1.5)
