###############################################################################
###############################################################################
#
# Quantifying the relationship between population trends and freshwater habitat
# Stephanie Peacock <speacock@psf.ca>
#
# This code loads to .rds file of MCMC output from model fits and creates some
# lookup lists that make finding the output for certain parameters easier.
# It is sourced by "lookingAtOutput.R" and "vulnerability.R"
#
###############################################################################
###############################################################################

# Load packages
library(gplots)
library(dclone)
library(sjmisc) # for is_even() function

# Source functions like getmode()
source("functions.R")

# Colour palettes for plotting
library(PNWColors)
spawnCol <- pnw_palette("Bay", n = 5)[c(1,2,4,5)]
mazCol <- pnw_palette("Bay", n = 7)
fazCol <- pnw_palette("Bay", n = 22)

#--------------------------------------------------------------------
# Setup MCMC approach
#--------------------------------------------------------------------
numChains <- 3
numIter <- 50000
numUpdate <- 20000
numAdapt <- 20000

iter <- c(1:numIter)
#--------------------------------------------------------------------
# Setup data list
#--------------------------------------------------------------------
# Load population data used in fitting
source("loadPopData.R")

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

habNames2 <- c("%Agri", "%Urban", "RipDist", "LinDev", "ForestRds", "nonForestRds", "StreamXing", "ForestDist", "ECA", "MPB")

regionNames <- sort(unique(popDat$region))
speciesNames <- sort(unique(popDat$SPECIES))

nSpecies <- length(speciesNames)
nPop <- dim(popDat)[1]
nFAZ <- length(fazNames)
nMAZ <- length(mazNames)
nSpawn <- length(spawnNames)
nRear <- length(rearNames)
nHab <- length(habNames)


###############################################################################
# Calculate Rhat, Summary output, and delist MCMC chains
###############################################################################

# fit <- readRDS("output/fit_model10_centeredSO_20210512.rds") #This one actually convered better
# fit <- readRDS("output/fit_model10_centeredSO_20210609_3chains0.rds") # Three longer chains had total convergence; 0 is for the fits with negative pressures corrected to zero
fit <- readRDS("output/fit_coho3_20211029.rds") # Fit from Oct 2021 changing coho age-at-return to 3 years

rHat <- gelman.diag(fit) # Takes a long time on new Mac...?
cn <- colnames(fit[[1]])

S <- summary(fit)

#-----------------------------------------------------------------------------
# Delist results for easy summary of range etc. among chains
#-----------------------------------------------------------------------------
out <- array(NA, dim = c(length(fit), dim(fit[[1]])), dimnames = list(paste("chain", c(1:length(fit)), sep =""), NULL, cn))
for(i in 1:length(fit)){
	out[i, , ] <- fit[[i]]
}


# Create list of parameter indices for beta0 and beta1
outInd <- list(
	beta0 = which(cn == "beta0"),
	beta1 = array(NA, dim = c(nSpawn, nHab)),
	phi = grep("phi", cn),
	sigmaMAZ = which(cn == "sigmaMAZ"),
	sigma = which(cn == "sigma")
)

for(i in 1:nSpawn){
	for(j in 1:nHab){
		outInd$beta1[i, j] <- match(paste("beta1[", i, ",", j, "]", sep = ""), cn)
	}
}

# Random effect for region
outInd$thetaMAZ <- array(NA, dim = c(nMAZ, nRear))
for(i in 1:nMAZ){
	for(j in 1:nRear){
		outInd$thetaMAZ[i, j] <- match(paste("thetaMAZ[", i, ",", j, "]", sep = ""), cn)
	}
}

# Random effect for FAZ
outInd$sigmaFAZ <- numeric(nHab)
for(j in 1:nHab){ 
	outInd$sigmaFAZ[j] <- match(paste("sigmaFAZ[", j, "]", sep = ""), cn)
}

outInd$thetaFAZ <- array(NA, dim = c(nFAZ, nHab))
for(i in 1:nFAZ){
	for(j in 1:nHab){ 
		outInd$thetaFAZ[i, j] <- match(paste("thetaFAZ[", i, ",", j, "]", sep = ""), cn)
	}
}

rHat[[1]][outInd$beta1]

# Table of parameter estimates
# write.csv(S[[1]], "output/parEst.csv")
