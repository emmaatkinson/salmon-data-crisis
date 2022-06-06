###############################################################################
###############################################################################
#
# Quantifying the relationship between population trends and freshwater habitat
# Stephanie Peacock <speacock@psf.ca>
#
# This code takes the model object from `fitting.R` and samples from the 
# posterior and data to calculate the sensitivity, exposure, and vulnerability 
# of different salmon species and Freshwater Adaptive Zones (FAZs) throughout BC.
# Produces large dot plots for main report showing these metrics.
#
###############################################################################
###############################################################################

# Source code that loads MCMC output and de-lists/creates index
# to easier access that output.

source("loadResults.R")

# fit = raw MCMC output from jags.fit()
# out = de-listed fit, with dimension (3, 50000, 333) = (chains, iterations, 
# parameters)

#--------------------------------------------------------------------
# What are the pressure increases  we're interested in?
#--------------------------------------------------------------------

# Summary of indicators used in analysis
hab <- read.csv("data/pse_habitatpressurevalues_2018_disaggregated_grid14fixed.csv")
dd <- read.csv("data/forest_disturbance_dd.csv")

# Remove Data Deficient watersheds for forest disturbance
hab <- hab[which(hab$WTRSHD_FID %in% dd$wtrshd_fid == FALSE), ]


hab <- hab[, match(habNames, names(hab))]

# Change MPB NAs to zero
hab$MPB_pct[is.na(hab$MPB_pct)] <- 0

# Set hab indicators slightly less than zero to zero
hab[which(hab < 0, arr.ind = TRUE)] <- 0

# Summarize mean, min, max, and quantiles for each habitat indicator
habSummary <- data.frame(
	indicator = habNames,
	mean = apply(hab, 2, mean, na.rm = TRUE),
	min = apply(hab, 2, min, na.rm = TRUE),
	t(apply(hab, 2, quantile, c(0.025, 0.25, 0.50, 0.75, 0.975))),
	max = apply(hab, 2, max, na.rm = TRUE)
)

habSummaryTable <- data.frame(
	indicator = habNames,
	distribution = paste0(sprintf('%.2f', apply(hab, 2, mean, na.rm = TRUE)), " (",  sprintf('%.2f', apply(hab, 2, quantile, 0.5, na.rm = TRUE)), ", ", sprintf('%.2f', apply(hab, 2, quantile, 0.975, na.rm = TRUE)), ")")
)
# write.csv(habSummaryTable, file = "habSummaryTable.csv")

max_inc <- numeric(nHab); names(max_inc) <- habNames
for(i in 1:nHab) max_inc[i] <- quantile(hab[, i], 0.975)

#--------------
# Supplemental figure showing distributions of pressure values
#--------------
# quartz(width = 8, height = 5, pointsize = 10)
# par(mfrow = c(2,5), mar = c(4,2,1,1), oma =c(0,3,1,0))
# for(i in 1:nHab){
# 	if(i %in% den) breaks <- seq(0, max(hab[, i]), length.out = 100) else breaks = c(-1:100)
# 	hist(hab[which(hab[, i] > 0), i], main = habNames2[i], xlab = "", ylab = "", breaks = breaks, border =NA, yaxs = "i", xaxs = "i")
# 	# hist(hab[which(hab[, i] > 0), i], add = TRUE, col = "#FF000050", border = NA, breaks = breaks)
# 	u <- par('usr')
# 	arrows(max_inc[i], u[3] + 0.6*(u[4] - u[3]), max_inc[i], u[3], col = 2, length = 0.08)
# 	text(max_inc[i], u[3] + 0.6*(u[4] - u[3]), round(max_inc[i], 2), xpd = NA, pos = 4, col = 2, font = 2)
# 	
# 	arrows(habSummary$mean[i], u[3] + 0.8*(u[4] - u[3]), habSummary$mean[i], u[3], col = 4, length = 0.08)
# 	text(habSummary$mean[i], u[3] + 0.8*(u[4] - u[3]), round(habSummary$mean[i], 2), xpd = NA, pos = 4, col = 4, font = 2)
# 	
# 	arrows(habSummary$max[i], u[3] + 0.4*(u[4] - u[3]), habSummary$max[i], u[3], length = 0.08)
# 	text(habSummary$max[i], u[3] + 0.4*(u[4] - u[3]), round(habSummary$max[i], 2), xpd = NA, pos = 2, font = 2)
# 	}
# mtext(side = 1, outer = TRUE, "Pressure value", line = -1)
# mtext(side = 2, outer = TRUE, "Number of watersheds", line = 1)

# Index to be able to track which indicators are percentages and which
# are densities
pct <- c(1, 2, 3, 8, 9, 10)
den <- c(4, 5, 6, 7)

# v9: Increase all pressure indicators as a percentage of the 75th quantile
pressure_inc <- seq(0, 1, 0.02) %*% t(max_inc)
max_x <- numeric(length(habNames)); names(max_x) <- habNames
max_x[pct] <- 100
max_x[den] <- habSummary$max[den] + pressure_inc[den]

# For density, increase from current to 100% * most impacted (i.e., habSummary$max)

#--------------------------------------------------------------------
# Setup parameters for random selection of MCMC output
#--------------------------------------------------------------------

n <- 10000 # Number of draws

# Select MCMC draws from model output
# For each of the n draws, which chain and iteration are we using?
set.seed(4569)
indMod <- cbind(
	sample(1:numChains, size = n, replace = TRUE), 
	sample(1:numIter, size= n, replace = TRUE))

#--------------------------------------------------------------------
# How many populations in each category?
#--------------------------------------------------------------------

nSlope <- array(
	data = NA, 
	dim = c(4, 10, nFAZ), 
	dimnames = list(spawnNames, habNames, fazNames))

for(i in 1:nSpawn){
	for(h in 1:nHab){
		for(j in 1:nFAZ){
			nSlope[i, h, j] <- length(which(popDat$spawnEco == spawnNames[i] & popDat[, habNames[h]] > 0 & popDat$FAZ == fazNames[j]))
		}
	}
}

nSlope2 <- array(
	data = NA, 
	dim = c(5, 10, nFAZ), 
	dimnames = list(speciesNames, habNames, fazNames))

for(s in 1:nSpecies){
	for(h in 1:nHab){
		for(i in 1:nFAZ){
			nSlope2[s, h, i] <- length(which(popDat$SPECIES == speciesNames[s] & popDat[, habNames[h]] > 0 & popDat$FAZ == fazNames[i]))
		}
	}
}

nPop <- array(
	data = NA, 
	dim = c(5, nFAZ), 
	dimnames = list(speciesNames, fazNames))

for(s in 1:nSpecies){
	for(i in 1:nFAZ){
			nPop[s, i] <- length(which(popDat$SPECIES == speciesNames[s] & popDat$FAZ == fazNames[i]))
	}
}
#--------------------------------------------------------------------
# Calculate by FAZ and species or ecotype:
#   1) sensitivity - degree to which fish populations respond to a pressure (slope)
#                 (beta[1,s,j] + theta[FAZ,j] + phi[j] * o)
#   2) exposure - pressure values by species, FAZ, and habitat indicator
#              x[i,j]
#   3) threat - sensitivity * exposure
#            (beta[1,s,j] + theta[FAZ,j] + phi[j] * o) * x[i,j]
#   4) threatTotal - threat summed across pressures
#   5) status - baseline trend
#                      beta0 + theta[MAZ,r]
#   6) vulnerability - status + overallThreat
#--------------------------------------------------------------------

y_pred <- array(
	data = NA, 
	dim = c(2, nSpecies, nFAZ, nHab, n), 
	dimnames = list(c("current_pressure", "increase_pressure"), speciesNames, fazNames, habNames, NULL))


#--------------------------------------------------------------------
# Calculate sensitivity and exposure
#--------------------------------------------------------------------

for(s in 1:nSpecies){
	
	for(i in 1:nFAZ){
		
		# Select MCMC draws from data
		ind.is <- which(popDat$SPECIES == speciesNames[s] & popDat$FAZ == fazNames[i])
		if(length(ind.is) > 0){
			indDat <- sample(ind.is, size = n, replace = TRUE)
			
			status <- out[cbind(indMod[, 1], indMod[, 2], rep(outInd$beta0, n))] + out[cbind(indMod[, 1], indMod[, 2], outInd$thetaMAZ[cbind(JAGSdat$maz[indDat], JAGSdat$rearEco[indDat])])]
				
				for(h in 1:nHab){ # For each habitat indicator
					
					Beta <- out[cbind(indMod[, 1], indMod[, 2], outInd$beta1[cbind(JAGSdat$spawnEco[indDat], rep(h, n))])] + out[cbind(indMod[, 1], indMod[, 2], rep(outInd$phi[h], n))] * JAGSdat$streamOrder[indDat] + out[cbind(indMod[, 1], indMod[, 2], outInd$thetaFAZ[cbind(JAGSdat$faz[indDat], rep(h, n))])]
					
					# Calculate sum of other pressures
					BetaXOther <- matrix(0, nrow = n, ncol = nHab)
					for(H in which((c(1:nHab) %in% h) == FALSE)){
						
						BetaXOther[, H] <- (out[cbind(indMod[, 1], indMod[, 2], outInd$beta1[cbind(JAGSdat$spawnEco[indDat], rep(H, n))])] + out[cbind(indMod[, 1], indMod[, 2], rep(outInd$phi[H], n))] * JAGSdat$streamOrder[indDat] + out[cbind(indMod[, 1], indMod[, 2], outInd$thetaFAZ[cbind(JAGSdat$faz[indDat], rep(H, n))])]) * JAGSdat$habPressures[cbind(indDat, rep(H, n))]
					} # end H
					
					# Notes: 
					# Use JAGSdat stream order, which is re-cetnered around 4, 
					# rather than popDat stream order
					# Use JAGSdat variables for spawnEco and rearEco, which are numeric
					
					# Pressure value for focal indicator:
					x <- JAGSdat$habPressures[cbind(indDat, rep(h, n))]
					
					# Sim. current pressure
					y_pred[1, s, i, h, ] <- status + Beta * x + apply(BetaXOther, 1, sum)
					
					# Sim. increased pressure
					y_pred[2, s, i, h, ] <- status + Beta * max_inc[h] + apply(BetaXOther, 1, sum)
	
					
					# # Check plot
					# breaks <- seq(min(y_pred[, s, i, h, ]), max(y_pred[, s, i, h, ]), length.out = 100)
					# hist(y_pred[1, s, i, h, ], breaks = breaks)
					# hist(y_pred[2, s, i, h, ], add = TRUE, col = "#FF000040", border = "#00000040", breaks = breaks)

				} #end h
				
		} # end if
	} # end i FAZ
} # end s species

# saveRDS(y_pred, file = "output/predictedTrend_forRiskAssessment.rds")
y_pred <- readRDS("output/predictedTrend_forRiskAssessment.rds")

#--------------------------------------------------------------------
# Calculate % threatened and endangered
#--------------------------------------------------------------------

# Classification as threatened according to the COSEWIC criteria of 
# a 10% decline within 10 years or 
# y_i = log(S[t]/S[t-1]) = log(0.9) = -0.1053605
# or withing 3 generations, with generation length G
# y_i = log(S)

# Declines over 10 years or 3 generations, whichever is longer
G <- pmax(10, 3 * c(5, 4, 3, 2, 5))
names(G) <- speciesNames

# Cutoffs
threatened <- log(0.5)/G
endangered <- log(0.3)/G

pThreat <- array(
	data = NA, 
	dim = c(2, nSpecies, nFAZ, nHab), 
	dimnames = list(
		c("current_pressure", "increased_pressure"), 
		speciesNames, 
		fazNames, 
		habNames))

for(s in 1:nSpecies){
	for(i in 1:nFAZ){
		for(p in 1:2){
			for(h in 1:nHab){
				pThreat[p, s, i, h] <- length(which(y_pred[p, s, i, h, ] < threatened[s]))/n
				}}}}


#------------------------------------------------------------------------------
# Table: Which species, FAZ's and threats are sensitive
#------------------------------------------------------------------------------

# What is the *change* in probability of being threatened for a 50% increase in 
# the pressure value?

threatIncrease <- array(NA, dim = c(nSpecies, nFAZ, nHab), dimnames = list(speciesNames, fazNames, habNames))
for(s in 1:nSpecies){
	for(i in 1:nFAZ){
		for(h in 1:10){
			threatIncrease[s, i, h] <- pThreat[2, s, i, h] - pThreat[1, s, i, h]
		}}}

# More than >=50% increase in threatened with a 50% increase in pressure
# and >80% ofpopulations threated with a 50% increase in pressure
# If there is currently no pressure, mark with asterix

inds <- which(threatIncrease >= 0.5, arr.ind = TRUE)
inc <- list(
	bySpecies = tapply(inds[, 1], inds[, 1], length),
	byFAZ = tapply(inds[, 1], inds[, 2], length),
	byThreat = tapply(inds[, 1], inds[, 3], length))
names(inc$bySpecies) <- speciesNames
names(inc$byFAZ) <- fazNames[as.numeric(names(inc$byFAZ))]
names(inc$byThreat) <- habNames[as.numeric(names(inc$byThreat))]


hist(pThreat[cbind(rep(2, nrow(inds)), inds)])


inds2 <- inds[which(pThreat[cbind(rep(2, nrow(inds)), inds)] > 0.8), ]
dim(inds)

# dims: s, i, h
# Create table of just those with pThreat > 0.7 -> inds2
tableInd <- inds2

popInd_inc <- list()
currrentPressureValues <- array(NA, dim = c(nrow(tableInd), 3), dimnames = list(NULL, c("mean", "min", "max")))
for(i in 1:nrow(tableInd)){
	popDat_inc[[i]] <- which(popDat$SPECIES == speciesNames[tableInd[i,1]] & popDat$FAZ == fazNames[tableInd[i,2]])
	currrentPressureValues[i, 1] <- mean(JAGSdat$habPressures[popDat_inc[[i]], tableInd[i,3]])
	currrentPressureValues[i, 2] <- min(JAGSdat$habPressures[popDat_inc[[i]], tableInd[i,3]])
	currrentPressureValues[i, 3] <- max(JAGSdat$habPressures[popDat_inc[[i]], tableInd[i,3]])
	
}

# Put together
threatTable <- data.frame(
	species = speciesNames[tableInd[, 1]],
	FAZ = fazNames[tableInd[, 2]],
	indicator = habNames4[tableInd[, 3]],
	
	nP = nPop[tableInd[, c(1,2)]],
	nSlope = nSlope2[tableInd[,c(1,3,2)]],
	
	currrentPressureValue = paste0(sprintf('%.2f', currrentPressureValues[, 1]), " (", sprintf('%.2f', currrentPressureValues[, 2]), ", ", sprintf('%.2f', currrentPressureValues[, 3]), ")"),
	curppnThreatened = pThreat[cbind(rep(1, nrow(tableInd)), tableInd)],
	
	incPressureValye = sprintf('%.2f', max_inc[tableInd[, 3]]),
	incppnThreatened = pThreat[cbind(rep(2, nrow(tableInd)), tableInd)],
	
	change = threatIncrease[tableInd]
)
threatTable <- threatTable[order(threatTable$species, threatTable$change, decreasing = c(FALSE, TRUE), method = "radix"), ]

write.csv(threatTable, file = "output/threatTable.csv")

###############################################################################
# Plotting
###############################################################################

	
#------------------------------------------------------------------------------
# Maps, for each species and impact
#------------------------------------------------------------------------------

# Load spatial packages
library(PBSmapping)
gshhg <- "~/Documents/Mapping/gshhg-bin-2.3.7/"
xlim <- c(-135, -118) + 360
ylim <- c(48, 58)

# five resolutions: crude(c), low(l), intermediate(i), high(h), and full(f).
res <- "i"
land <- importGSHHS(paste0(gshhg,"gshhs_", res, ".b"), xlim = xlim, ylim = ylim, maxLevel = 2, useWest = TRUE)
rivers <- importGSHHS(paste0(gshhg,"wdb_rivers_", res, ".b"), xlim = xlim, ylim = ylim, useWest = TRUE)
borders <- importGSHHS(paste0(gshhg,"wdb_borders_", res, ".b"), xlim = xlim, ylim = ylim, useWest = TRUE, maxLevel = 1)

faz <- as.PolySet(read.csv("data/ignore/PSF/fazLL_thinned_inDat.csv"), projection = "LL")
cbind(tapply(faz$FAZ_Acrony, faz$PID, unique), fazNames)
fazIndex <- unique(faz$PID)[match(fazNames, tapply(faz$FAZ_Acrony, faz$PID, unique))]


# faz <- thinPolys(faz, tol = 1)
# Darker red = higher proportion threatened
colThreat <- data.frame(
	threatIncrease = seq(-0.99, 1, 0.02),
	colour = colorRampPalette(c("#0f85a0", "white", "#dd4124"))(n = 100))
	
habNames4 <- c("%Agriculture", "%Urban Devel.", "%Riparian Dist.", "Linear Devel.", "Forestry Roads", "Non-forestry Roads", "Stream Crossings", "%Forest Dist.", "%ECA", "%Pine Beetle")

pdf(file = "figures/threatenedMaps.pdf", width = 8, height = 11, pointsize = 10)
for(s in 1:nSpecies){
	par(mfrow = c(4,3), mar = c(0,0,0,0), oma = c(1,1,3,1))
	for(h in 1:10){

		# Plot map, no border...?
		plotMap(land, col = "white", bg = NA, las = 1, border = NA, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
		
		for(i in 1:nFAZ){
			if(nSlope2[s,h,i] > 0){
				addPolys(faz[faz$PID == fazIndex[i], ], border = NA, col = colThreat$colour[findInterval(threatIncrease[s, i, h], colThreat$threatIncrease)], lwd = 0.6)
			} else {
				addPolys(faz[faz$PID == fazIndex[i], ], border = grey(0.8), col = grey(0.8), density = 25, lwd = 0.6)
			}
		}
		
		addLines(rivers, col = grey(0.6), lwd = 0.6)
		addPolys(land, border = grey(0.6), lwd = 0.3)
		# addLines(borders)
			
		mtext(side = 3, line = 0.5, paste(letters[h], ") ", habNames4[h], sep = ""), cex = 1, adj = 0)
				
	}
	mtext(side = 3, outer= TRUE, speciesNames[s])
	
	# Add legend
	plot(1,1,"n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
	# polygon(x = c(0.5, 3, 3, 0.5), y = c(1, 1, 1.3, 1.3), xpd = NA)
	for(i in seq(1, 100, 2)){
		polygon(x = c(0.45, 0.5, 0.5, 0.45) + i*0.02, y = c(1, 1, 1.3, 1.3), xpd = NA, col = colThreat$colour[i], border = NA)
		
	}
	segments(x0 = 0.45 + c(seq(1, 100, 10), 101)*0.02, y0 = 0.95, x1 = 0.45 + c(seq(1, 100, 10), 101)*0.02, 1, xpd = NA)
		text(0.45 + c(seq(1, 100, 10), 101)*0.02, 0.9, round(c(colThreat$threatIncrease[seq(1, 100, 10)] - 0.01, 1), 2), xpd = NA)
		text(1.5, 1.4, "Change in proportion threatened with 50% increase in pressure",xpd = NA, cex = 1.2)
}
dev.off()

##############################################################################
# Focal populations
##############################################################################
fazCentroid <- calcCentroid(faz, rollup = 1)
fazCentroid$FAZ_Acrony <- faz$FAZ_Acrony[match(fazCentroid$PID, faz$PID)]

#------------------------------------------------------------------------------
# Focal populations: Coho and Stream crossings
#------------------------------------------------------------------------------
s <- 3
h <- 7

# pdf(file = "threatenedMap_cohoStreamX.pdf", width = 4, height = 5, pointsize = 10)
quartz(width = 8, height = 4, pointsize = 10)
layout(matrix(c(1,1,2,3,4, 1,1,5,6,7), nrow = 2, byrow = TRUE))
par(mar = c(5,3,8,1), oma = c(1, 0, 0, 1))
plotMap(land, col = "white", bg = NA, las = 1, border = NA, xlab = "", ylab = "")

for(i in 1:nFAZ){
	if(nSlope2[s,h,i] > 0){
		addPolys(faz[faz$PID == fazIndex[i], ], border = grey(0.8), col = colThreat$colour[findInterval(threatIncrease[s, i, h], colThreat$threatIncrease)], lwd = 0.6)
	} else {
		addPolys(faz[faz$PID == fazIndex[i], ], border = grey(0.8), col = grey(0.8), density = 25, lwd = 0.6)
	}
}
addLines(rivers, lwd = 0.6)
addPolys(land, lwd = 0.6)
addLines(borders, col = grey(0.6), lwd = 1.2)
mtext(side = 3, line = 1.5, "a) Change in the proportion\nof coho populations that are threatened\nwith a 50% increase in stream crossing density", adj = 0)
# dev.off()

# Focal FAZs
par(mar = c(5, 6, 8, 1))
# quartz(width = 3, height = 2, pointsize = 10)
# par(mfrow = c(1,1), mar = c(3,3,1,1)))
focalFAZ <- match(c("LFR", "LILL"), fazNames)
# i <- 16
j <- 0
for(i in focalFAZ){
	j <- j + 1
	ind.is <- which(popDat$SPECIES == speciesNames[s] & popDat$FAZ == fazNames[i])
	indDat <- sample(ind.is, size = n, replace = TRUE)
	Beta <- out[cbind(indMod[, 1], indMod[, 2], outInd$beta1[cbind(JAGSdat$spawnEco[indDat], rep(h, n))])] + out[cbind(indMod[, 1], indMod[, 2], rep(outInd$phi[h], n))] * JAGSdat$streamOrder[indDat] + out[cbind(indMod[, 1], indMod[, 2], outInd$thetaFAZ[cbind(JAGSdat$faz[indDat], rep(h, n))])]
	
	hist(Beta, freq = FALSE, col = NA, border = grey(0.8),  main = "", yaxs = "i", xlim = c(-0.15, 0.05), breaks = seq(-0.3, 0.1, 0.01))
	lines(density(Beta), lwd = 2)
	abline(v = 0)
	mtext(side = 3, line = 0.5, adj = 0, paste0(c("b", "e")[j], ")"))
	mtext(side = 3, line = 3, fazNames[i], adj = 0, font = 2)
	
	xR <- range(y_pred[c(1, which(pressure_inc == 0.50)), s, i, h, ])
	hist(y_pred[1, s, i, h, ], breaks = seq(xR[1], xR[2], length.out = 30), yaxs = "i", main = "", xlab = "Predicted trend", xlim = c(-2.5, 1))
	hist(y_pred[which(pressure_inc == 0.50), s, i, h, ], add = TRUE, col = "#dd412450", border = "#dd4124", breaks = seq(xR[1], xR[2], length.out = 30))
	abline(v = threatened[s], lty = 2)
	mtext(side = 3, line = 0.5, adj = 0, paste0(c("c", "f")[j], ")"))
	
	plot(pressure_inc, pThreat[1, , s, i, h], "l", xlab = "Pressure increase", ylab = "Proportion threatened", ylim = c(0,1), bty = "l", xaxs = "i")
	segments(x0 = 0.5, x1 = 0.5, y0 = 0, y1 = pThreat[1, which(pressure_inc == 0.5), s, i, h], lty = 2, col = "#dd4124")
	segments(x0 = 0.5, x1 = 0, y0 = pThreat[1, which(pressure_inc == 0.5), s, i, h], y1 = pThreat[1, which(pressure_inc == 0.5), s, i, h], lty = 2, col = "#dd4124")
	segments(x0 = 0, x1 = 0.5, y0 = pThreat[1, 1, s, i, h], y1 = pThreat[1, 1, s, i, h], lty = 2, col = grey(0.6))
	mtext(side = 3, line = 0.5, adj = 0, paste0(c("d", "g")[j], ")"))

}

segments(x0 = c(-3.781119, -3.683241), x1 = c(-2.700872, -2.732900), y0 = c(1.3386811, 0.8162466), y1 = c(2.9864691, 1.2579034), lwd = 1.2, xpd = NA)

#------------------------------------------------------------------------------
# Chinook and % ECA
#------------------------------------------------------------------------------
s <- 1
h <- 9

# pdf(file = "threatenedMap_cohoStreamX.pdf", width = 4, height = 5, pointsize = 10)
quartz(width = 8, height = 4, pointsize = 10)
layout(matrix(c(1,1,2,3,4, 1,1,5,6,7), nrow = 2, byrow = TRUE))
par(mar = c(5,3,8,1), oma = c(1, 0, 0, 1))
plotMap(land, col = "white", bg = NA, las = 1, border = NA, xlab = "", ylab = "")

for(i in 1:nFAZ){
	if(nSlope2[s,h,i] > 0){
		addPolys(faz[faz$PID == fazIndex[i], ], border = grey(0.8), col = colThreat$colour[findInterval(threatIncrease[s, i, h], colThreat$threatIncrease)], lwd = 0.6)
	} else {
		addPolys(faz[faz$PID == fazIndex[i], ], border = grey(0.8), col = grey(0.8), density = 25, lwd = 0.6)
	}
}
addLines(rivers, lwd = 0.6)
addPolys(land, lwd = 0.6)
addLines(borders, col = grey(0.6), lwd = 1.2)
mtext(side = 3, line = 1.5, "a) Change in the proportion\nof Chinook populations that are threatened\nwith a 50% increase in %ECA", adj = 0)
# dev.off()

# Focal FAZs
par(mar = c(5, 6, 8, 1))
focalFAZ <- match(c("LNR-P", "NC"), fazNames)
# i <- 16
j <- 0
for(i in focalFAZ){
	j <- j + 1
	ind.is <- which(popDat$SPECIES == speciesNames[s] & popDat$FAZ == fazNames[i])
	indDat <- sample(ind.is, size = n, replace = TRUE)
	Beta <- out[cbind(indMod[, 1], indMod[, 2], outInd$beta1[cbind(JAGSdat$spawnEco[indDat], rep(h, n))])] + out[cbind(indMod[, 1], indMod[, 2], rep(outInd$phi[h], n))] * JAGSdat$streamOrder[indDat] + out[cbind(indMod[, 1], indMod[, 2], outInd$thetaFAZ[cbind(JAGSdat$faz[indDat], rep(h, n))])]
	
	hist(Beta, freq = FALSE, col = NA, border = grey(0.8),  main = "", yaxs = "i", xlim = c(-0.15, 0.05), breaks = seq(-0.3, 0.1, 0.01))
	lines(density(Beta), lwd = 2)
	abline(v = 0)
	mtext(side = 3, line = 0.5, adj = 0, paste0(c("b", "e")[j], ")"))
	mtext(side = 3, line = 3, fazNames[i], adj = 0, font = 2)
	
	xR <- range(y_pred[c(1, which(pressure_inc == 0.50)), s, i, h, ])
	hist(y_pred[1, s, i, h, ], breaks = seq(xR[1], xR[2], length.out = 30), yaxs = "i", main = "", xlab = "Predicted trend", xlim = c(-2.5, 1))
	hist(y_pred[which(pressure_inc == 0.50), s, i, h, ], add = TRUE, col = "#dd412450", border = "#dd4124", breaks = seq(xR[1], xR[2], length.out = 30))
	abline(v = threatened[s], lty = 2)
	mtext(side = 3, line = 0.5, adj = 0, paste0(c("c", "f")[j], ")"))
	
	plot(pressure_inc, pThreat[1, , s, i, h], "l", xlab = "Pressure increase", ylab = "Proportion threatened", ylim = c(0,1), bty = "l", xaxs = "i")
	segments(x0 = 0.5, x1 = 0.5, y0 = 0, y1 = pThreat[1, which(pressure_inc == 0.5), s, i, h], lty = 2, col = "#dd4124")
	segments(x0 = 0.5, x1 = 0, y0 = pThreat[1, which(pressure_inc == 0.5), s, i, h], y1 = pThreat[1, which(pressure_inc == 0.5), s, i, h], lty = 2, col = "#dd4124")
	segments(x0 = 0, x1 = 0.5, y0 = pThreat[1, 1, s, i, h], y1 = pThreat[1, 1, s, i, h], lty = 2, col = grey(0.6))
	mtext(side = 3, line = 0.5, adj = 0, paste0(c("d", "g")[j], ")"))
	
}

segments(x0 = c(-3.781119, -3.683241), x1 = c(-2.700872, -2.732900), y0 = c(1.3386811, 0.8162466), y1 = c(2.9864691, 1.2579034), lwd = 1.2, xpd = NA)




#------------------------------------------------------------------------------
# Look at beta, and distribution of trends for certain combinations of 
# species, FAZ and indicator
#------------------------------------------------------------------------------
pdf(file = "figures/mostVulnerable36.pdf", width = 6.3, height = 2.8, pointsize = 10)

# for(j in 1:nrow(inds)){
# s <- inds[j, 1]
# i <- inds[j, 2]
# h <- inds[j, 3]

s <- which(speciesNames == "Chinook")
i <- which(fazNames == "EVI")
h <- which(habNames == "AgriculturePCT")

quartz(width = 6.3, height = 2.8, pointsize = 10)
par(mfrow = c(1,3), mar = c(4,4,2,1), oma = c(0,0,3,0))

ind.is <- which(popDat$SPECIES == speciesNames[s] & popDat$FAZ == fazNames[i])
indDat <- sample(ind.is, size = n, replace = TRUE)

Beta0 <- out[cbind(indMod[, 1], indMod[, 2], outInd$beta1[cbind(JAGSdat$spawnEco[indDat], rep(h, n))])] + out[cbind(indMod[, 1], indMod[, 2], rep(outInd$phi[h], n))] * JAGSdat$streamOrder[indDat]

Beta <- out[cbind(indMod[, 1], indMod[, 2], outInd$beta1[cbind(JAGSdat$spawnEco[indDat], rep(h, n))])] + out[cbind(indMod[, 1], indMod[, 2], rep(outInd$phi[h], n))] * JAGSdat$streamOrder[indDat] + out[cbind(indMod[, 1], indMod[, 2], outInd$thetaFAZ[cbind(JAGSdat$faz[indDat], rep(h, n))])]

d <- density(Beta)
d0 <- density(Beta0)
# hist(Beta, freq = FALSE, col = NA, border =NA,  main = "", yaxs = "i", xlim = c(-0.15, 0.05), breaks = seq(-0.3, 0.1, 0.01))
plot(d$x, d$y, "n", main = "", yaxs = "i", yaxt = "n", bty = "n", ylab = "", xlab = "Beta")
mtext(side = 2, line = 1, "Marginal density", cex = 0.8)
u <- par('usr')
abline(v = 0)
arrows(x0 = u[1], x1 = u[1], y0 = u[3], y1 = u[4], length = 0.08, xpd = NA)
polygon(x = c(d0$x, rev(d0$x)), y = c(d0$y, rep(0, length(d0$y))), col = "#00000010", border = grey(0.6))
polygon(x = c(d$x, rev(d$x)), y = c(d$y, rep(0, length(d$y))), col = "#dd412430", border = "#dd4124")
mtext(side = 3, line = 0.5, adj = 0, "d)")
legend("topleft", fill = c("#00000010", "#dd412430"), border = c(grey(0.6), "#dd4124"), bty = "n", c("Overall", fazNames[i]))

#--
# xR2 <- quantile(y_pred[c(1, which(pressure_inc == 0.50)), s, i, h, ], c(0.025, 0.975))
xR <- range(y_pred[c(1, which(pressure_inc == 0.50)), s, i, h, ])
hist(y_pred[1, s, i, h, ], breaks = seq(xR[1], xR[2], length.out = 30), yaxs = "i", main = "", xlab = "Predicted trend in spawner abundance", border = grey(0.6))
hist(y_pred[which(pressure_inc == 0.50), s, i, h, ], add = TRUE, col = "#00000060", border = 1, breaks = seq(xR[1], xR[2], length.out = 30))
abline(v = threatened[s], lty = 2)
mtext(side = 3, line = 0.5, adj = 0, "e)")
legend("topleft", fill = c(grey(0.8), "#00000060"), border = c(grey(0.6), 1), bty = "n", c("Current pressure", "50% increase"))

plot(pressure_inc, pThreat[1, , s, i, h], "l", xlab = "% increase in pressure value", ylab = "Proportion threatened", ylim = c(0,1), bty = "l", xaxs = "i")
segments(x0 = 0.5, x1 = 0.5, y0 = 0, y1 = pThreat[1, which(pressure_inc == 0.5), s, i, h], lty = 2)
segments(x0 = 0.5, x1 = 0, y0 = pThreat[1, which(pressure_inc == 0.5), s, i, h], y1 = pThreat[1, which(pressure_inc == 0.5), s, i, h], lty = 2)
# segments(x0 = 0, x1 = 0.5, y0 = pThreat[1, 1, s, i, h], y1 = pThreat[1, 1, s, i, h], lty = 2, col = grey(0.6))
axis(side =2, at = pThreat[1, 1, s, i, h], col = grey(0.6), las =1)

mtext(side = 3, line = 0.5, adj = 0, "f)")

mtext(side = 3, outer = TRUE, paste("Vulnerability of", speciesNames[s], "in", fazNames[i], "to", habNames4[h]), line = 1)

# }
# 
# dev.off()
#------------------------------------------------------------------------------
# Across all FAZs, weighted by the number of populations in each FAZ
#------------------------------------------------------------------------------

pThreat.avg <- array(NA, dim = c(nSpecies, nHab, length(pressure_inc)))
for(s in 1:nSpecies){
	for(h in 1:nHab){
		# Weighted by number of populations
		pThreat.avg[s, h, ] <- apply(pThreat[1, , s, , h], 1, weighted.mean, w = nSlope2[s, h, ])
		
	}
}

speciesCol <- pnw_palette("Bay", n = nSpecies)

h <- 1
par(mfrow = c(4,3))
for(h in 1:nHab){
	plot(pressure_inc, pThreat.avg[1, h, ], "n", ylim = c(-1, 1))
	for(s in 1:nSpecies) lines(pressure_inc, pThreat.avg[s, h, ], col = speciesCol[s])
	mtext(side = 3, line = 0.5, paste(letters[h], ") ", habNames4[h], sep = ""), cex = 1, adj = 0)
}

legend(1.5, 0, xpd = NA, col = speciesCol, lwd = 1, legend = speciesNames)

par(mfrow = c(2, 2), mar = c(4,4,2,1))
for(h in c(4, 5, 6, 7))


# What is the starting distribution of trends?
par(mfrow = c(1, nSpecies), mar = c(4,2,1,1))
for(s in 1:nSpecies){
	hist(popDat$trend3[popDat$SPECIES == speciesNames[s]], col = paste0(speciesCol[s], 30), border = speciesCol[s], breaks = seq(-0.5, 0.8, 0.1), main = speciesNames[s])
	abline(v = 0)
	abline(v = threatened[s], lty = 2)
}










###############################################################################
# Chat with Doug May 3, 2022
###############################################################################
s <- which(speciesNames == "Coho")
h <- which(habNames == "STRMXDEN")
i <- which(fazNames == "LFR")
# Lower Fraser Coho

hist(y_pred[1, s, i, h, ], yaxs = "i", main = "", xlab = "Predicted trend", breaks = seq(-0.2, 0.35, 0.01))
# hist(y_pred[which(pressure_inc == 0.50), s, i, h, ], add = TRUE, col = "#dd412450", border = "#dd4124", breaks = seq(xR[1], xR[2], length.out = 30))
abline(v = threatened[s], lty = 2, lwd = 2)
abline(v = 0, lwd = 2)
par(new = TRUE)
hist(popDat$trend3[which(popDat$FAZ == "LFR" & popDat$SPECIES == "Coho")], col = "#dd412450", border = "#dd4124", breaks = seq(-0.2, 0.35, 0.01), yaxs = "i", main = "", yaxt = "n")
axis(side = 4, col = 2)

legend("topright", fill = c(grey(0.8), "#dd4124"), legend =c("Predicted", "Observed"), bty ="n")

popDat[which(popDat$FAZ == "LFR" & popDat$SPECIES == "Coho"), ]
