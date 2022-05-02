###############################################################################
# Fitting approach #1:
# a) Use previously estimated trends that were calculated in 
#    dataCompilationNuSEDS_CU.R
#    using simple lm(log(smoothedMAX_ESTIMATE + 1) ~ ANALYSIS_YR)
#    for each population.
# b) Fit model relating population trends to habitat pressure in JAGS to allow
#    for more flexibility and diagnostics.
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

source("models.R")

model <- model_10

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
