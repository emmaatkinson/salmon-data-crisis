###############################################################################
###############################################################################
#
# Testing for broad-scale relationships between freshwater habitat pressure 
# indicators and Pacific salmon population trends
# Stephanie Peacock <speacock@psf.ca>
#
# Model fitting approach:
# a) Use previously estimated trends that were calculated in 
#    dataCompilationNuSEDS_CU.R
#    using simple lm(log(smoothedMAX_ESTIMATE + 1) ~ ANALYSIS_YR)
#    for each population.
# b) Fit model relating population trends to habitat pressure in JAGS to allow
#    for more flexibility and diagnostics.
#
###############################################################################
###############################################################################

library(dclone)
source("functions.R")

###############################################################################
# Read in population data and transform to dat list for JAGS
###############################################################################

source("loadPopData.R")

#------------------------------------------------------------------------------
# Define list of variables for JAGS
#------------------------------------------------------------------------------

JAGSdat <- list(
	
	nPop = dim(popDat)[1],
  nFAZ = length(fazNames),
	nMAZ = length(mazNames),
	nSpawn = length(spawnNames),
	nRear = length(rearNames),
	nHab = length(habNames),
	trend = popDat$trend3,
 	streamOrder = popDat$StreamOrder - getmode(popDat$StreamOrder), # Centre stream order around zero
	faz = as.numeric(factor(popDat$FAZ, levels = fazNames)),
	maz = as.numeric(factor(popDat$MAZ, levels = mazNames)),
	rearEco = as.numeric(factor(popDat$rearEco, levels = rearNames)),
	spawnEco = as.numeric(factor(popDat$spawnEco, levels = spawnNames)),
	
	habPressures = habPressures,
	weight = 1/(1 + exp(popDat$QUALITY - 5))
)

###############################################################################
# Define model
###############################################################################

model <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
		}
		
		phi[j] ~ dnorm(0, pow(3, -2))
		
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		
		for(j in 1:nHab){
			slopeBeta[i, j] <- beta1[spawnEco[i], j] + phi[j] * streamOrder[i] + thetaFAZ[faz[i], j]
		}
		
		trendPred[i] <- beta0	+ thetaMAZ[maz[i], rearEco[i]] + inprod(slopeBeta[i, ], t(habPressures[i, ]))
		
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * weight[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * weight[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model



###############################################################################
# Fitting
###############################################################################

# Parameters to track
parNames <- c(
	"log_lik",
	"beta0", # length 1
	"beta1", # length nHab
  "phi",
	"sigmaMAZ", # length 1
	"thetaMAZ", # length nPFMA
  "sigmaFAZ", # length nHab
	"thetaFAZ", # matrix nPop x nHab
	"sigma" # length 1
)

numChains <- 3
numIter <- 50000
numUpdate <- 20000
numAdapt <- 20000
# Fit model to subsetted data
cl <- makeCluster(numChains)

t.start<-proc.time()[3]

fit <- jags.parfit(
	cl,
	data = JAGSdat,
	params = parNames,
	model = model,
	n.adapt = numAdapt,
	n.update = numUpdate,
	n.iter = numIter,
	n.chains = numChains)

stopCluster(cl)
procTime <- round((proc.time()[3] - t.start)/(60), 1) # track process time
cat(paste("Process time for model =", procTime, "minutes"))


# saveRDS(fit, file = "output/fit_coho3_20211029.rds")
