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

# Different increases for indicators that are percent of watershed
# from those that are density

# Summary of indicators used in analysis
hab <- read.csv("data/pse_habitatpressurevalues_2018_disaggregated_grid14fixed.csv")
dd <- read.csv("data/forest_disturbance_dd.csv")
# Remove Data Deficient watersheds for forest disturbance
hab <- hab[which(hab$WTRSHD_FID %in% dd$wtrshd_fid == FALSE), ]
hab <- hab[, match(habNames, names(hab))]
hab$MPB_pct[is.na(hab$MPB_pct)] <- 0
hab[which(hab < 0, arr.ind = TRUE)] <- 0
habSummary <- data.frame(
	indicator = habNames, 
	mean = apply(hab, 2, mean, na.rm = TRUE), 
	min = apply(hab, 2, min, na.rm = TRUE),
	t(apply(hab, 2, quantile, c(0.025, 0.25, 0.50, 0.75, 0.975))),
	max = apply(hab, 2, max, na.rm = TRUE)
)

write.csv(habSummary, file = "habSummary.csv")
# For percent, increase from 0 to 100%
pct <- c(1, 2, 3, 8, 9, 10)
den <- c(4, 5, 6, 7)

pressure_inc <- seq(0, 1, 0.02) 
max_x <- c(rep(100, 3), 2*habSummary$max[den], rep(100, 3)); names(max_x) <- habNames

# For density, increase from current to 100% * most impacted (i.e., habSummary$max)

#--------------------------------------------------------------------
# Setup parameters for random selection of MCNC output
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
	dim = c(length(pressure_inc), nSpecies, nFAZ, nHab, n), 
	dimnames = list(NULL, speciesNames, fazNames, habNames, NULL))


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
					
					# Sim
					for(p in 1:length(pressure_inc)){
						
						if(h %in% den){ # If pressure value is a density
							# Add 0 - 100% of the most impacted watershed
							y_pred[p, s, i, h, ] <- status + Beta * pmin(x + (pressure_inc[p])*habSummary$max[h], max_x[h]) + apply(BetaXOther, 1, sum)
						
							} else { # If percent watershed, then add
							# Add percentages up to a max of 100%
								y_pred[p, s, i, h, ] <- status + Beta * pmin(x + 100*pressure_inc[p], max_x[h]) + apply(BetaXOther, 1, sum)
						}
					} # end p
					
					# # Check plot
					# par(mfrow = c(2,2), mar = c(4,4,2,1))
					# for(p in 1:4){
					# 	hist(y_pred[p, s, i, h, ], xlim = range(y_pred[1:3, s, i, h, ]), xlab = "", ylab = "", main = c("current", "+10%", "+20%", "+50%")[p])
					# 	abline(v = c(threatened, endangered), col = 2, lwd = 1.5, lty  = c(2,1))
					# }
					# 
					

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
	dim = c(2, length(pressure_inc), nSpecies, nFAZ, nHab), 
	dimnames = list(
		c("threatened", "endangered"), 
		NULL, 
		speciesNames, 
		fazNames, 
		habNames))

for(s in 1:nSpecies){
	for(i in 1:nFAZ){
		for(p in 1:length(pressure_inc)){
			for(h in 1:nHab){
				pThreat[1, p, s, i, h] <- length(which(y_pred[p, s, i, h, ] < threatened[s]))/n
				pThreat[2, p, s, i, h] <- length(which(y_pred[p, s, i, h, ] < endangered[s]))/n
			}}}}


###############################################################################
# Plotting
###############################################################################
habCol <- pnw_palette("Bay", nHab)

s <- 5
i <- 13
par(mfrow = c(1,2))

# Threatened
plot(pressure_inc, pThreat[1, , s, i, 1], "n", ylim = c(0,1))
for(h in 2:10) lines(pressure_inc, pThreat[1, , s, i, h], col = habCol[h])

# Endangered
plot(pressure_inc, pThreat[2, , s, i, 1], "n", ylim = c(0,1))
for(h in 2:10) lines(pressure_inc, pThreat[2, , s, i, h], col = habCol[h])

#------------------------------------------------------------------------------
# Plot each indicator separately 
#------------------------------------------------------------------------------

quartz(width = 7, height = 6, pointsize = 10)
par(mfrow = c(3,4), mar = c(3,3,0,0), oma = c(2,2,3,2))
for(h in 1:10){
	plot(pressure_inc, pThreat[1, , s, i, h], "l", ylim = c(0,1))
	lines(pressure_inc, pThreat[2, , s, i, h], col = 2)
	mtext(side = 3, line = -2, habNames2[h])
}
mtext(side = 3, outer=TRUE, paste0(speciesNames[s], ": ", fazNames[i]), line = 1, cex = 1.5)
plot(1,1,"n", bty = "n", xaxt = "n", yaxt = "n")
legend("center", col = c(1,2), pch = 1, lwd = 1, c("Threatened", "Endangered"), bty = "n")
mtext(side = 1, outer = TRUE, "Percent increase in pressure value")
mtext(side = 2, outer = TRUE, "Proportion of populations threatened or endangered")

#------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------
pdf(width = 7, height = 6, pointsize = 10)
par(mfrow = c(3,4), mar = c(3,3,0,0), oma = c(2,2,3,2))

for(s in 1:nSpecies){
	for(i in 1:nFAZ){
		for(h in 1:10){
			
			# Plot only those that have a change in proportion of > 5%
			if(tail(pThreat[1, , s, i, h], 1) - head(pThreat[1, , s, i, h], 1) >= 0.10 & tail(pThreat[1, , s, i, h], 1) > 0.6){
				plot(pressure_inc, pThreat[1, , s, i, h], "l", ylim = c(0,1))
				lines(pressure_inc, pThreat[2, , s, i, h], col = 2)
				mtext(side = 3, line = -1, paste0(speciesNames[s], ": ", fazNames[i], ": ", habNames2[h]), cex = 0.7)
			}
}}}
dev.off()
	
#------------------------------------------------------------------------------
# Maps, for each species and impact
#------------------------------------------------------------------------------

# What is the *change* in probability of being threatened for a 50% increase in 
# the pressure value?

threatIncrease <- array(NA, dim = c(nSpecies, nFAZ, nHab), dimnames = list(speciesNames, fazNames, habNames))
for(s in 1:nSpecies){
	for(i in 1:nFAZ){
		for(h in 1:10){
			threatIncrease[s, i, h] <- pThreat[1, which(pressure_inc == 0.50), s, i, h] - pThreat[1, which(pressure_inc == 0), s, i, h]
		}}}

# Load spatial packages
library(PBSmapping)
gshhg <- "~/Google Drive/Mapping/gshhg-bin-2.3.7/"
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
# Table: Which species, FAZ's and threats are sensitive
#------------------------------------------------------------------------------
# More than >=50% increase in threatened with a 50% increase in pressure
# and >80% ofpopulations threated with a 50% increase in pressure
# If there is currently no pressure, mark with asterix

inds <- which(threatIncrease >= 0.5, arr.ind = TRUE)
inc <- tapply(inds[, 1], inds[,1], length)
names(inc) <- speciesNames
inc

inds <- inds[which(pThreat[cbind(rep(1, nrow(inds)), rep(which(pressure_inc == 0.5), nrow(inds)), inds)] > 0.8), ]
# dims: s, i, h
dim(inds)

dim(nSlope2)
nSlope2[inds[,c(1,3,2)]]

threatTable <- data.frame(
	species = speciesNames[inds[, 1]],
	FAZ = fazNames[inds[, 2]],
	indicator = habNames4[inds[, 3]],
	currentppnThreatened = pThreat[cbind(rep(1, nrow(inds)), rep(which(pressure_inc == 0), nrow(inds)), inds)],
	ppnThreatened50 = pThreat[cbind(rep(1, nrow(inds)), rep(which(pressure_inc == 0.5), nrow(inds)), inds)],
	increase = threatIncrease[inds],
	nP = nPop[inds[, c(1,2)]],
	nSlope = nSlope2[inds[,c(1,3,2)]]
)
threatTable <- threatTable[order(threatTable$species, threatTable$increase, decreasing = c(FALSE, TRUE), method = "radix"), ]

write.csv(threatTable, file = "output/threatTable.csv")

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
