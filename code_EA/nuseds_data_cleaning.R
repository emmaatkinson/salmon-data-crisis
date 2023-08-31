##############################################################################
###############################################################################
#
# Assembling spawner abundance data to visualise trends in monitoring
# coverage and data quality for wild Pacific salmon in BC
# Emma Atkinson <emmamargaretatkinson@gmail.com>
#
# Compiling escapement data and performing quality control
#
# Data were downloaded from Open Data Canada
# https://open.canada.ca/data/en/dataset/c48669a3-045b-400d-b730-48aafe8c5ee6
#
##############################################################################
##############################################################################
rm(list=ls())

library(here)
setwd(here())

source("code_EA/functions.R")

#------------------------------------------------------------------------------
# Read in NuSEDS data
#------------------------------------------------------------------------------

# (0) Original NuSEDS dataset, without any CU information
nuseds <- read.csv("data_EA/ignore/NuSEDS_2023-June.csv")
pse_decoder <- read.csv("data_EA/pse-region_dfo-cu_decoder.csv")

# Clean up Estimate Classification field of nuseds
nuseds$ESTIMATE_CLASSIFICATION[nuseds$ESTIMATE_CLASSIFICATION == "PRESENCE/ABSENCE (TYPE-6)"] <- "PRESENCE-ABSENCE (TYPE-6)"
nuseds$ESTIMATE_CLASSIFICATION[nuseds$ESTIMATE_CLASSIFICATION =="NO SURVEY THIS YEAR"] <- "NO SURVEY"
nuseds$ESTIMATE_CLASSIFICATION[nuseds$ESTIMATE_CLASSIFICATION ==""] <- "UNKNOWN"

# Re-level ESTIMATE_CLASSIFICATION for easy numeric calculations
nuseds$ESTIMATE_CLASSIFICATION <- factor(
	nuseds$ESTIMATE_CLASSIFICATION,
	levels = c(
		"TRUE ABUNDANCE (TYPE-1)",
		"TRUE ABUNDANCE (TYPE-2)",
		"RELATIVE ABUNDANCE (TYPE-3)",
		"RELATIVE ABUNDANCE (TYPE-4)",
		"RELATIVE ABUNDANCE (TYPE-5)",
		"PRESENCE-ABSENCE (TYPE-6)",
		"RELATIVE: CONSTANT MULTI-YEAR METHODS",
		"RELATIVE: VARYING MULTI-YEAR METHODS",
		"NO SURVEY",
		"UNKNOWN"))

# (1) NuSEDS dataset organized by CU
# This dataset has fewer rows than the full NuSEDS, but more info on lat/lon
# and watershed, etc. (as well as CU!), and the pink salmon are split into
# even and odd year populations, so its preferable.

nusedsCU <- read.csv("data_EA/ignore/Conservation_Unit_Data_20220902_2023-June.csv")

# Create unique id for each "SPECIES QUALIFIED" population
# I.e., different for even and odd year pinks
nusedsCU$SQ_POP_ID <- paste(nusedsCU$SPECIES_QUALIFIED, nusedsCU$POP_ID, sep = "")
popCount()

# Order in a sensible way
nusedsCU <- nusedsCU[order(nusedsCU$SQ_POP_ID, nusedsCU$ANALYSIS_YR), ]

# Clean up Estimate Classification field of nuseds
nusedsCU$ESTIMATE_CLASSIFICATION[nusedsCU$ESTIMATE_CLASSIFICATION == "PRESENCE/ABSENCE (TYPE-6)"] <- "PRESENCE-ABSENCE (TYPE-6)"
nusedsCU$ESTIMATE_CLASSIFICATION[nusedsCU$ESTIMATE_CLASSIFICATION =="NO SURVEY THIS YEAR"] <- "NO SURVEY"
nusedsCU$ESTIMATE_CLASSIFICATION[nusedsCU$ESTIMATE_CLASSIFICATION ==""] <- "UNKNOWN"

# Re-level ESTIMATE_CLASSIFICATION for easy numeric calculations
nusedsCU$ESTIMATE_CLASSIFICATION <- factor(
	nusedsCU$ESTIMATE_CLASSIFICATION, 
	levels = c(
		"TRUE ABUNDANCE (TYPE-1)", 
		"TRUE ABUNDANCE (TYPE-2)", 
		"RELATIVE ABUNDANCE (TYPE-3)", 
		"RELATIVE ABUNDANCE (TYPE-4)", 
		"RELATIVE ABUNDANCE (TYPE-5)", 
		"PRESENCE-ABSENCE (TYPE-6)", 
		"RELATIVE: CONSTANT MULTI-YEAR METHODS", 
		"RELATIVE: VARYING MULTI-YEAR METHODS", 
		"NO SURVEY", 
		"UNKNOWN"))

nusedsCU0 <- nusedsCU # Kepp a version without any filtering for reference
nusedsCU <- nusedsCU0
popFiltered <- data.frame(Step = 0, Descr = "Raw NuSEDS CU", nRemoved = 0, nPop = popCount())

#------------------------------------------------------------------------------
# (1) Filter out entries with MAX_ESTIMATE == NA
#------------------------------------------------------------------------------

nusedsCU <- nusedsCU[ - which(is.na(nusedsCU$MAX_ESTIMATE)), ] # Removes 157045 records

popFiltered <- recordStep(popFiltered, newDescr = "Remove MAX_ESTIMATE == NA")

#------------------------------------------------------------------------------
# (2) Remove years prior to when reliable monitoring started
#------------------------------------------------------------------------------

#-----
# For all populations, truncate to 1950
nusedsCU <- nusedsCU[ - which(nusedsCU$ANALYSIS_YR < 1950), ]
popFiltered <- recordStep(popFiltered, newDescr = "Truncate to 1950")

#-----
# Fraser sockeye: population specific based on input from Keri Benner
# Southern BC Chinook: population specific input from Brown et al. (2020)

sy <- read.csv("data_EA/startYears.csv")
sy <- sy[!is.na(sy$FirstYear),]

# # Plot for supplement
# par(mfrow = c(1,1), mar = c(4,5,2,10))
# plot(sy$FirstYear_nuseds, sy$FirstYear, pch = c(1,2)[as.numeric(factor(sy$Source))], xlab = "First year in NuSEDS data", ylab = "", col = c(1,2)[as.numeric(sy$FirstYear < sy$FirstYear_nuseds) + 1], las = 1)
# abline(a =0, b = 1)
# mtext(side = 2, line = 4, "First year of reliable data")
# legend(2015, 2010, pch = c(1,2,1), col = c(1,1,2), c("Brown et al. (2020)", "K. Benner", "Not used"), xpd = NA, bty = "n")
# sy[which(sy$FirstYear_nuseds > sy$FirstYear), ]
# pp(45014)
# abline(v = sy$FirstYear[sy$POP_ID == 45014])
# pp(46328)
# abline(v = 1992.5)


for(i in 1:dim(sy)[1]){
	# If that population is in the dataset
	if(sy$POP_ID[i] %in% nusedsCU$POP_ID){
		rmYrs <- which(nusedsCU$ANALYSIS_YR < sy$FirstYear[i] & nusedsCU$POP_ID == sy$POP_ID[i])
		if(length(rmYrs) > 0){
			nusedsCU <- nusedsCU[ - rmYrs, ]		
		}
	}
}

popCount() #6881 # EA (Aug 2023): 6886 (?)

popFiltered <- recordStep(popFiltered, newDescr = "Truncate sockeye, Chinook based on additional info")

# Import cu to determine GFE_TYPE
cu <- read.csv("data_EA/ignore/conservation_unit_system_sites_2023-June.csv")

#------------------------------------------------------------------------------
# (4) Disaggregation
#------------------------------------------------------------------------------

# Are there duplicates for single population- year combination?
boop <- tapply(paste(nusedsCU$SQ_POP_ID, nusedsCU$ANALYSIS_YR, sep = ":"), paste(nusedsCU$SQ_POP_ID, nusedsCU$ANALYSIS_YR, sep = ":"), length)
unique(boop)
names(boop[which(boop == 2)])

#_____
# Remove extra years (see NuSEDS_CU_checks.R for in-depth look at what's going on here)
# Barriere River (Clearwater) Coho Run 1 had two estimates: (1) FENNELL CREEK AND SASKUM CREEK
# and (2) BARRIERE RIVER
# Remove the BARRIERE RIVER estimate for now
nusedsCU <- nusedsCU[ - which(nusedsCU$ANALYSIS_YR == 2013 & nusedsCU$SQ_POP_ID == "CO46602" & nusedsCU$STREAM_ID == 258), ]

#_____
## Add North Thomson Coho tributaries together then remove duplicates
# Includes NORTH THOMPSON RIVER, BIRCH ISLAND CHANNEL, PIG CHANNEL COMBINED
ny <- tapply(nusedsCU$MAX_ESTIMATE[nusedsCU$SQ_POP_ID == "CO46582"], nusedsCU$ANALYSIS_YR[nusedsCU$SQ_POP_ID == "CO46582"], length)
unique(nusedsCU$GEOGRAPHICAL_EXTNT_OF_ESTIMATE[nusedsCU$SQ_POP_ID == "CO46582"])
for(i in which(ny == 2)){
	y <- as.numeric(names(ny)[i])
	nusedsCU$MAX_ESTIMATE[nusedsCU$SQ_POP_ID == "CO46582" & nusedsCU$ANALYSIS_YR == y] <- sum(nusedsCU$MAX_ESTIMATE[nusedsCU$SQ_POP_ID == "CO46582" & nusedsCU$ANALYSIS_YR == y])
	ind.i <- which(nusedsCU$SQ_POP_ID == "CO46582" & nusedsCU$ANALYSIS_YR == y & nusedsCU$STREAM_ID == 719256264)
	nusedsCU <- nusedsCU[-ind.i, ]
	print(paste("removed row", ind.i, "\n"))
}

#_____
# Add tributaries (?) for Alouette River (Coquitlam) Chum CM 47925
ny <- tapply(nusedsCU$MAX_ESTIMATE[nusedsCU$SQ_POP_ID == "CM47925"], nusedsCU$ANALYSIS_YR[nusedsCU$SQ_POP_ID == "CM47925"], length)

for(i in which(ny == 2)){
	y <- as.numeric(names(ny)[i])
	nusedsCU$MAX_ESTIMATE[nusedsCU$SQ_POP_ID == "CM47925" & nusedsCU$ANALYSIS_YR == y] <- sum(nusedsCU$MAX_ESTIMATE[nusedsCU$SQ_POP_ID == "CM47925" & nusedsCU$ANALYSIS_YR == y])
	ind.i <- which(nusedsCU$SQ_POP_ID == "CM47925" & nusedsCU$ANALYSIS_YR == y & nusedsCU$STREAM_ID == 15)
	nusedsCU <- nusedsCU[-ind.i, ]
	print(paste("removed row", ind.i, "\n"))
}

#_____
# Add recent Widgeon Slough (Coquitlam) Summer Sockeye data from NuSEDS that is missing
# from the CU dataset
# Create dummy rows to modify based on NuSEDS info for 2015-2021
addYrs <- rbind(
	nusedsCU[nusedsCU$POP_ID == 47954 & nusedsCU$ANALYSIS_YR == 2013, ],
	nusedsCU[nusedsCU$POP_ID == 47954 & nusedsCU$ANALYSIS_YR == 2013, ],
	nusedsCU[nusedsCU$POP_ID == 47954 & nusedsCU$ANALYSIS_YR == 2013, ],
	nusedsCU[nusedsCU$POP_ID == 47954 & nusedsCU$ANALYSIS_YR == 2013, ],
	nusedsCU[nusedsCU$POP_ID == 47954 & nusedsCU$ANALYSIS_YR == 2013, ],
	nusedsCU[nusedsCU$POP_ID == 47954 & nusedsCU$ANALYSIS_YR == 2013, ],
	nusedsCU[nusedsCU$POP_ID == 47954 & nusedsCU$ANALYSIS_YR == 2013, ])

# Update ANALYSIS_YR field to reflect missing years
addYrs$ANALYSIS_YR <- c(2015:2021)

# Columns to match in nuseds
colNames <- names(addYrs)[which(!is.na(match(names(addYrs), names(nuseds))))]
for(y in 2015:2021){
	addYrs[which(addYrs$ANALYSIS_YR == y), match(colNames, names(addYrs))] <- nuseds[which(nuseds$ANALYSIS_YR == y & nuseds$POP_ID == 47954), match(colNames, names(nuseds))]
}

# Add extra years to NuSEDS CU data
nusedsCU <- rbind(nusedsCU, addYrs)

#------
# Check that all duplicates have been dealt with:
unique(tapply(paste(nusedsCU$SQ_POP_ID, nusedsCU$ANALYSIS_YR, sep = ":"), paste(nusedsCU$SQ_POP_ID, nusedsCU$ANALYSIS_YR, sep = ":"), length)) # SHould be one.

popFiltered <- recordStep(popFiltered, newDescr = "Removed duplicate records, dis/aggregated populations")

###############################################################################
# Data QUALITY screening
###############################################################################

#------------------------------------------------------------------------------
# What does the data quality look like over time?
estClassAll <- array(NA, dim = c(length(1950:2021), length(levels(nusedsCU$ESTIMATE_CLASSIFICATION))), dimnames = list(c(1950:2021), levels(nusedsCU$ESTIMATE_CLASSIFICATION)))
for(i in 1:length(1950:2021)){
	dum0 <- nusedsCU$ESTIMATE_CLASSIFICATION[nusedsCU$ANALYSIS_YR == c(1950:2021)[i]]
	estClassAll[i, ] <- tapply(dum0, dum0, length)
}
estClassAll[which(is.na(estClassAll) == TRUE, arr.ind = TRUE)] <- 0

pdf("figures/EstimateClassificationTypes_Aug2023.pdf", width = 7, height = 4.2, pointsize = 8)
#windows(width = 6.3, height = 4, pointsize = 8)
par(mar = c(5,5,5,1))
bp <- barplot(t(estClassAll), beside = FALSE, las = 1, col = estCol, border = NA, xaxt = "n", ylab = "Number of populations", ylim=c(0,4500))
axis(side = 1, at = bp[seq(1, length(bp), 2)], labels = FALSE, tck = -0.01)
axis(side = 1, at = bp[seq(1, length(bp), 10)], labels = c(1950:2021)[seq(1, length(bp), 10)], tck = -0.01, tck = -0.03)

# plot(1,1, "n", xaxt="n", yaxt = "n", xlab = "", ylab = "", bty = "n")
legend(-1, 4500, ncol = 2, fill = estCol, legend = levels(nusedsCU$ESTIMATE_CLASSIFICATION), bg = "white", border = NA, xpd = NA)
dev.off()

#------------------------------------------------------------------------------


nusedsCU$ESTIMATE_CLASSIFICATION2 <- nusedsCU$ESTIMATE_CLASSIFICATION
nusedsCU$ESTIMATE_CLASSIFICATION2[which(as.numeric(nusedsCU$ESTIMATE_CLASSIFICATION) > 6)] <- NA

x <- tapply(as.numeric(nusedsCU$ESTIMATE_CLASSIFICATION2), nusedsCU$SQ_POP_ID, mean, na.rm = TRUE)
pdf("figures/EstimateClassificationHistogram_Aug2023.pdf", width = 7, height = 4.2, pointsize = 8)
hist(x, breaks = seq(1, 10, 0.2), col = colorRampPalette(c("forest green", "gold", "red", "black"))(n = 46), main = "", las = 1, xlab = "mean estimate classification type", ylab = "Number of populations")
dev.off()

# Create numeric score with NU SURVEY = UNKNOWN = 9
nusedsCU$ESTIMATE_CLASSIFICATION3 <- as.numeric(nusedsCU$ESTIMATE_CLASSIFICATION)
nusedsCU$ESTIMATE_CLASSIFICATION3[which(nusedsCU$ESTIMATE_CLASSIFICATION3 == 7)] <- 5
nusedsCU$ESTIMATE_CLASSIFICATION3[which(nusedsCU$ESTIMATE_CLASSIFICATION3 == 8)] <- 6
nusedsCU$ESTIMATE_CLASSIFICATION3[which(nusedsCU$ESTIMATE_CLASSIFICATION3 >= 9)] <- 7

# How to weight the data quality based on estimate classification?
pdf("figures/Weighting_Aug2023.pdf", width = 3.2, height = 2.3, pointsize = 8)
par(mar = c(4,4,1,1))
plot(seq(1, 9, 0.1), 1/(1 + exp((seq(1, 9, 0.1) - 5))), "l", xlab = "Estimate classification type", ylab = "Weights (w)", las = 1, ylim = c(0, 1), xlim = c(1, 7))
points(1:7, 1/(1+exp((1:7 - 5))), pch = 19, cex = 2, col = estCol[c(1:6, 10)])
abline(h = 0.5, lty = 2)
dev.off()

nusedsCU$weight <- 1/(1 + exp(nusedsCU$ESTIMATE_CLASSIFICATION3 - 5))

# What is the mean estimate classification since 2000?
# If there is no estimate since 2000, then use 1990-2000. 
# If no estimate then, then go back another decade.
nusedsCU$QUAL <- numeric(dim(nusedsCU)[1])
nPop <- length(unique(nusedsCU$SQ_POP_ID))
spid <- unique(nusedsCU$SQ_POP_ID)
for(i in 1:nPop){
	p.i <- which(nusedsCU$SQ_POP_ID == spid[i])
	nusedsCU.i <- nusedsCU[p.i, ]
	m <- max(nusedsCU.i$ANALYSIS_YR)
	if(m >= 2000){
		nusedsCU$QUAL[p.i] <- mean(nusedsCU.i$ESTIMATE_CLASSIFICATION3[nusedsCU.i$ANALYSIS_YR >= 2000])
	} else if(m >= 1990){
		nusedsCU$QUAL[p.i] <- mean(nusedsCU.i$ESTIMATE_CLASSIFICATION3[nusedsCU.i$ANALYSIS_YR >= 1990])
	} else if(m >= 1980){
		nusedsCU$QUAL[p.i] <- mean(nusedsCU.i$ESTIMATE_CLASSIFICATION3[nusedsCU.i$ANALYSIS_YR >= 1980])
	} else if(m >= 1970){
		nusedsCU$QUAL[p.i] <- mean(nusedsCU.i$ESTIMATE_CLASSIFICATION3[nusedsCU.i$ANALYSIS_YR >= 1970])
	} else {
		print(paste("Warning: No estimates for SQ_POP_ID", spid[i]))
		nusedsCU$QUAL[p.i] <- 7 # Changed from 9 to be the max est class type
	}
}
# THere are no estimates for a number of populations, but these all have unknown estimate
# classification so fine to assign to 10
sum(is.na(nusedsCU$QUAL))

###############################################################################
# Calculate smoothed spawner abundance
###############################################################################

nPop <- popCount()
spid <- unique(nusedsCU$SQ_POP_ID)
popInd <- numeric(nPop)
for(i in 1:nPop)popInd[i] <- which(nusedsCU$SQ_POP_ID == spid[i])[1]
pid <- nusedsCU$POP_ID[popInd]
# 
# # Import generation length information from PSE
# genL <- read.csv("data/All_regions_CU_decoder.csv")
# genL$gen_length[genL$species == "CO"] <- 3 # Changed from 4 to 3 in Oct 2021
# 
# # If missing generation length by CU, use the following
# gen.all <- c(Chum = 4, Chinook = 5, Coho = 3, Pink = 2, Sockeye = 5)
# 
# # Smooth MAX_ESTIMATE using geometric mean over generation length
# nusedsCU$smoothedMAX_ESTIMATE <- rep(NA, dim(nusedsCU)[1])
# nusedsCU$smoothedweightedMAX_ESTIMATE <- rep(NA, dim(nusedsCU)[1])
# 
# ptm <- proc.time() # this can take a few minutes
# for(i in 1:nPop){ # For each population
# 	
# 	ind.i <- which(nusedsCU$SQ_POP_ID == spid[i])
# 	# Subset MAX_ESTIMATE and ANALYSIS_YR (and weights)
# 	n <- nusedsCU$MAX_ESTIMATE[ind.i]
# 	y <- nusedsCU$ANALYSIS_YR[ind.i]
# 	w <- nusedsCU$weight[ind.i]
# 	
# 	# Set generation length
# 	cu_in <- cu$FULL_CU_IN[which(cu$POP_ID == pid[i] & cu$SPECIES_QUALIFIED == nusedsCU$SPECIES_QUALIFIED[popInd][i])]
# 	g <- genL$gen_length[which(genL$FULL_CU_IN == cu_in)]
# 	if(length(g) == 0){
# 		g <- gen.all[nusedsCU$SPECIES[nusedsCU$SQ_POP_ID == spid[i]][1]]
# 	} else if(is.na(g)){
# 		g <- gen.all[nusedsCU$SPECIES[nusedsCU$SQ_POP_ID == spid[i]][1]]
# 	}
# 	
# 	# Calculate geometric average
# 	nSmoothed <- rep(NA, length(n))
# 	for(k in 1:length(y)){ # For each year
# 		smoothYrs <- c(max(y[1], y[k] - g + 1):y[k]) # define previous years over which to smooth
# 		
# 		ind.ik <- which(nusedsCU$SQ_POP_ID == spid[i] & nusedsCU$ANALYSIS_YR == y[k])
# 		
# 		# Unweighted geometric mean
# 		nusedsCU$smoothedMAX_ESTIMATE[ind.ik] <- prod(n[y %in% smoothYrs]) ^ (1/sum(y %in% smoothYrs))
# 		
# 		# Weighted geometric mean
# 		nusedsCU$smoothedweightedMAX_ESTIMATE[ind.ik] <- prod((n^w)[y %in% smoothYrs]) ^ (1 / sum(w[y %in% smoothYrs]))
# 		
# 	}
# } # End all pops smooth spawners
# proc.time() - ptm # 3.5 mins
# 
#  plot(nusedsCU$smoothedMAX_ESTIMATE, nusedsCU$smoothedweightedMAX_ESTIMATE, ylim = c(0, 1e6), xlim = c(0, 1e6))
#  abline(a = 0, b = 1, lty = 2)

 
###############################################################################
# Write .csv of nusedsCU with these QA/QC and calculations completed
###############################################################################

write.csv(nusedsCU, file = "data_EA/Conservation_Unit_Data_filtered_20230831.csv")
write.csv(popFiltered, file = "data_EA/Conservation_Unit_Data_QAQC_20230831.csv")

###############################################################################
###############################################################################
# Create population data list
###############################################################################
###############################################################################

#------------------------------------------------------------------------------
# Freshwater Adaptive zones for habitat
#------------------------------------------------------------------------------

faz <- cu$FAZ_ACRO[match(pid, cu$POP_ID)]

pid[which(is.na(faz))]
pid[which(is.na(faz))] %in% cu$POP_ID
# Four POP_IDs that we don't have info on because they are not in the CU site data
# EA update: now there are 11
unique(nusedsCU$POPULATION[nusedsCU$POP_ID %in% pid[which(is.na(faz))]])

cuAdd <- data.frame(
	MAP_LABEL = rep(NA, 4),
	GFE_ID = c("Stream", rep(NA, 3)),
	SYSTEM_SITE = c("CHILCOTIN RIVER", "KENNEDY LAKE", "KENNEDY LAKE", "CLAYOQUOT ARM BEACHES"),
	GFE_TYPE = rep(NA, 4),
	SPECIES_QUALIFIED = c("CK", "CO", "SEL", "SEL"),
	Y_LAT = c(51.7404, rep(49.042106, 2), 49.114138),
	X_LONGT = c(-122.4023, rep(-125.593115, 2), -125.560631),
	FAZ_ACRO = c("MFR", rep("WVI", 3)),
	MAZ_ACRO = c("GStr", rep("WVI", 3)),
	JAZ_ACRO = c("MFR+GStr", rep(NA, 3)),
	CU_NAME = c("MIDDLE FRASER RIVER_SP_1.3", "CLAYOQUOT", rep("KENNEDY", 2)),
	CU_ACRO = c("MFR-spring", "CLAY", rep("Kennedy", 2)),
	CU_LAT = c(52.72282, 49.26781, rep(49.11552, 2)), 
	CU_LONGT = c(-123.0937 , -125.8691, rep(-125.5372, 2)),
	CU_TYPE = rep("Current", 4),
	CU_INDEX = c("10", "18", rep("L-13-14", 2)),
	FULL_CU_IN = c("CK-10", "CO-18", rep("SEL-13-14", 2)),
	SBJ_ID = rep(NA, 4),
	POP_ID = c(46841, 52082, 52080, 52120),
	IS_INDICATOR = rep("N", 4),
	CMNTS = rep(NA, 4),
	EFFECTIVE_DT = rep(NA, 4),
	WATERSHED_CDE = rep(NA, 4),
	# NEWWATERSHEDCDE = rep(NA, 4),
	FWA_WATERSHED_CDE = rep(NA, 4)
)

cu <- rbind(cu, cuAdd)
faz <- cu$FAZ_ACRO[match(pid, cu$POP_ID)]

#------------------------------------------------------------------------------
# Marine Adaptive zones for habitat
#------------------------------------------------------------------------------

maz <- cu$MAZ_ACRO[match(pid, cu$POP_ID)]
pid[which(is.na(maz))]
# All have MAZ info (EA update: 8 missing info)

#------------------------------------------------------------------------------#
# At population level (fit trend only)
#------------------------------------------------------------------------------#
# Refresh population list in case of changes
if(nPop != popCount()){
	nPop <- popCount()
	spid <- unique(nusedsCU$SQ_POP_ID)
	popInd <- numeric(nPop)
	for(i in 1:nPop)popInd[i] <- which(nusedsCU$SQ_POP_ID == spid[i])[1]
	pid <- nusedsCU$POP_ID[popInd]
}

popDat <- data.frame(
	POPULATION = nusedsCU$POPULATION[popInd],
	CU_NAME = nusedsCU$CU_NAME[popInd],
	FULL_CU_IN = cu$FULL_CU_IN[match(pid, cu$POP_ID)],
	POP_ID = pid,
	SPECIES_QUALIFIED = nusedsCU$SPECIES_QUALIFIED[popInd],
	SPECIES = nusedsCU$SPECIES[popInd],
	AREA = nusedsCU$AREA[popInd],
	FAZ = cu$FAZ_ACRO[match(pid, cu$POP_ID)],
	MAZ = cu$MAZ_ACRO[match(pid, cu$POP_ID)],
	#StreamOrder = wsd$StreamOrder[match(pid, wsd$POP_ID)],
	QUALITY = nusedsCU$QUAL[popInd],
	weight = nusedsCU$weight[popInd],
	nYears = tapply(nusedsCU$ANALYSIS_YR, nusedsCU$SQ_POP_ID, length)#,
	#forestDD = as.numeric(pid %in% popForestDD)
)

###############################################################################
# Final CU corrections based on Eric Hertz input & adding column for PSE regions
############################################################################### 

popDat$CU_findex[popDat$POP_ID==3119] <- "CM-16"
popDat$CU_findex[popDat$POP_ID==51771] <- "CM-16"
popDat$CU_findex[popDat$POP_ID==51772] <- "CM-16"
popDat$CU_findex[popDat$POP_ID==3143] <- "CM-16"
popDat$CU_findex[popDat$POP_ID==3122] <- "CM-16"
popDat$CU_findex[popDat$POP_ID==3125] <- "CM-16"
popDat$CU_findex[popDat$POP_ID==3138] <- "CM-16"
popDat$CU_findex[popDat$POP_ID==3128] <- "CM-16"
popDat$CU_findex[popDat$POP_ID==51778] <- "CM-16" 

# Adding column for PSE region
popDat$PSE_region = NA

pse_dec = pse_decoder[which(!(is.na(pse_decoder$cu_index))),]
missing_from_decoder = NA


for (r in unique(popDat$FULL_CU_IN)){
	if (r %in% unique(pse_dec$cu_index)){
	popDat$PSE_region[popDat$FULL_CU_IN==r] = pse_dec$region[pse_dec$cu_index==r]
	} else {
			missing_from_decoder = c(missing_from_decoder, r)
		}
}

popDat[is.na(popDat$FULL_CU_IN),]
missing_from_decoder

#------------------------------------------------------------------------------#
# Calculate trend for each population
#------------------------------------------------------------------------------#

# popDat$trend1 <- rep(NA, nPop)
# popDat$trend2 <- rep(NA, nPop)
# popDat$trend3 <- rep(NA, nPop)
# popDat$trend4 <- rep(NA, nPop)
# popDat$trend5 <- rep(NA, nPop)
# 
# for(i in 1:nPop){
# 	ind <- which(nusedsCU$SQ_POP_ID == spid[i])
# 	
# 	y1 <- log(nusedsCU$smoothedMAX_ESTIMATE[ind] + 1)
# 	y3 <- log(nusedsCU$smoothedweightedMAX_ESTIMATE[ind] + 1)
# 	x <- nusedsCU$ANALYSIS_YR[ind]
# 	w <- nusedsCU$weight[ind]
# 	
# 	popDat$trend1[i] <- coefficients(lm(y1 ~ x))[2]
# 	popDat$trend2[i] <- coefficients(lm(y3 ~ x))[2]
# 	popDat$trend3[i] <- coefficients(lm(y3 ~ x, weights = w))[2]
# 	
# 	
# 	ind4 <- which(x > 1984)
# 	if(length(ind4) >= 10)	popDat$trend4[i] <- coefficients(lm(y3[ind4] ~ x[ind4]))[2]
# 	
# 	ind5 <- which(x > 1999)
# 	if(length(ind5) >= 10)	popDat$trend5[i] <- coefficients(lm(y3[ind5] ~ x[ind5]))[2]
# 	
# }
# 
# #------
# # Plots to look at the difference among trend methods
# plot(popDat$trend1, popDat$trend3); abline(a = 0, b = 1, lty = 2)
# 
# i <- 3457
# i <- 536
# pp(popDat$POP_ID[i], cex = 1.5, pch = 1)
# points(nusedsCU$ANALYSIS_YR[nusedsCU$POP_ID == 45043], nusedsCU$smoothedMAX_ESTIMATE[nusedsCU$POP_ID == 45043], "o", pch = 19, cex =0.8, col = grey(0.8))
# points(nusedsCU$ANALYSIS_YR[nusedsCU$POP_ID == 45043], nusedsCU$smoothedweightedMAX_ESTIMATE[nusedsCU$POP_ID == 45043], "o", pch = 17, cex = 0.6, col = grey(0.6))
# 
# colp <- estCol[as.numeric(nusedsCU$ESTIMATE_CLASSIFICATION3[nusedsCU$POP_ID == popDat$POP_ID[i]])]
# 
# plot(x, y1, col = colp, xlab = "Year", ylab = "log(S + 1)", pch = 19, xlim = c(1980, 2019))
# points(x, y3, col = colp, cex = 1.5)
# abline(coefficients(lm(y1 ~ x)))
# abline(coefficients(lm(y3 ~ x)), lty = 2)
# abline(coefficients(lm(y3 ~ x, weights = w)), lty = 3)
# abline(coefficients(lm(y3[ind4] ~ x[ind4])), lty = 4)
# abline(coefficients(lm(y3[ind5] ~ x[ind5])), lty = 5)
# abline(v = 1984.5, lty = 4)
# abline(v = 1999.5, lty = 5)
# legend("bottomleft", lty = c(1:5), legend = paste("Trend", 1:5))
# 
# hist(popDat$trend1, border = NA)
# hist(popDat$trend3, border = NA, col = "#FF000030", add = TRUE)


#------------------------------------------------------------------------------
# Region rather than PFMA for marine RE
#------------------------------------------------------------------------------
# popDat$region <- rep(NA, nPop)
# popDat$region[popDat$AREA %in% c("1", "2E", "2W")] <- "HaidaGwaii"
# popDat$region[popDat$AREA %in% c("3A", "3B", "4A", "4B", "4C", "4D")] <- "NorthCoast"
# popDat$region[popDat$AREA %in% c("5", "6", "7", "8", "9", "10", "11")] <- "CentralCoast"
# popDat$region[popDat$AREA %in% c("12", "13")] <- "InnerSouthCoast"
# popDat$region[popDat$AREA %in% c("14N", "14S", "15", "16", "17", "18", "19", "20", "28A", "28B", "29B", "29C", "29D", "29E", "29F", "29G", "29H", "29I", "29J", "29K")] <- "SoG"
# popDat$region[popDat$AREA %in% c("21",  "22",  "23",  "24", "25",  "26", "27")] <- "WCVI"

# Use regions from the PSE
popDat$region_2 <- hab$Region[match(wsd$WTRSHD_FID[match(nusedsCU$POP_ID[popInd], wsd$POP_ID)], hab$WTRSHD_FID)]

sum(is.na(popDat$region))


#------------------------------------------------------------------------------
# Spawning ecotypes related to spawning habitat only
#------------------------------------------------------------------------------
popDat$spawnEco <- popDat$SPECIES

# Pink/chum
popDat$spawnEco[popDat$spawnEco == "Pink"] <- "Pink/Chum"
popDat$spawnEco[popDat$spawnEco == "Chum"] <- "Pink/Chum"

#------------------------------------------------------------------------------
# Marine/rearing ecotypes
#------------------------------------------------------------------------------
popDat$rearEco <- popDat$SPECIES

# Differentiate stream and ocean type Chinook
popDat$rearEco[popDat$SPECIES == "Chinook" & grepl("0.", popDat$CU_NAME)] <- "ChinookOcean"
popDat$rearEco[popDat$SPECIES == "Chinook" & grepl("1.", popDat$CU_NAME)] <- "ChinookStream"

popDat$rearEco[which(popDat$SPECIES == "Chinook" & (popDat$FULL_CU_IN %in% c(paste("CK", c(34, 36, 38, 39, 41, 46, 82), sep = "-"))))] <- "ChinookOcean"

popDat$rearEco[which(popDat$SPECIES == "Chinook" & (popDat$FULL_CU_IN %in% c(paste("CK", c(37, 42, 43,48:51,53:58,80), sep = "-"))))] <- "ChinookStream"

# Dean River Chinook
popDat$rearEco[popDat$POP_ID == 51809] <-  "ChinookStream"

# Differentiate lake-spawning and river-spawning sockeye
popDat$rearEco[popDat$SPECIES_QUALIFIED == "SEL"] <- "SockeyeLake"
popDat$rearEco[popDat$SPECIES_QUALIFIED == "SER"] <- "SockeyeRiver"

unique(popDat$rearEco)

#------------------------------------------------------------------------------
# Habitat pressure values
#------------------------------------------------------------------------------
popDat <- cbind(popDat, hab[match(wsd$WTRSHD_FID[match(nusedsCU$POP_ID[popInd], wsd$POP_ID)], hab$WTRSHD_FID), c('AgriculturePCT', 'UrbanPCT', 'RIPDSTPCT', 'Linear_Dev_noRoads', 'ForestRoadsDEN_km_km2', 'NonForestRoadsDEN_km_km2', 'STRMXDEN', 'FORDST_PCT', 'ECAPCT', 'MPB_pct', 'WATER_LIC', 'WWD_count')])

write.csv(popDat, file = "data_EA/popDat_20230831.csv")
