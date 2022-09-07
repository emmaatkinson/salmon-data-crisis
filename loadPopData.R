###############################################################################
###############################################################################
#
# Testing for broad-scale relationships between freshwater habitat pressure 
# indicators and Pacific salmon population trends
# Stephanie Peacock <speacock@psf.ca>
#
# Load population data from output of dataCompilationNuSEDS_CU.R
# For use in model fitting and plotting.
###############################################################################
###############################################################################

###############################################################################
# Read in population data and transform to dat list for JAGS
###############################################################################

popDat <- read.csv("data/popDat.csv")

#------------------------------------------------------------------------------
# Remove populations for which forest disturbance is data deficient
# (i.e., >50% of the corresponding spawning watershed is privately owned)
#------------------------------------------------------------------------------
popDat <- popDat[which(popDat$forestDD == 0), ]

#------------------------------------------------------------------------------
# Replace NAs in mtn pine beetle with 0s
#------------------------------------------------------------------------------
popDat$MPB_pct[is.na(popDat$MPB_pct)] <- 0

#------------------------------------------------------------------------------
# Replace NonForestRoadsDEN_km_km2 < 0
#------------------------------------------------------------------------------
popDat$NonForestRoadsDEN_km_km2[which(popDat$NonForestRoadsDEN_km_km2 < 0)] <- 0

#------------------------------------------------------------------------------
# Replace outdated MAZ names
#------------------------------------------------------------------------------
popDat$MAZ[popDat$MAZ == "NQCI"] <- "NHG"
popDat$MAZ[popDat$MAZ == "WQCI"] <- "WHG"

#------------------------------------------------------------------------------
# Define levels for categorical factors
#------------------------------------------------------------------------------

habNames <- c("AgriculturePCT","UrbanPCT","RIPDSTPCT","Linear_Dev_noRoads", "ForestRoadsDEN_km_km2", "NonForestRoadsDEN_km_km2","STRMXDEN","FORDST_PCT","ECAPCT","MPB_pct")

rearNames <- sort(unique(popDat$rearEco))
spawnNames <- sort(unique(popDat$spawnEco))

# For FAZ and MAZ, sort by latitude
fazNames <- c("UNR", "USK", "MSK", "LNR-P", "LSK", "HG", "NC", "BCD", "HecLow", "RSI", "SC", "HK", "LILL", "EVI", "WVI", "LFR", "FRCany", "MFR", "UFR", "LTh", "STh", "NTh")
# # Check they're all in; should be zero
# sum(fazNames %in% unique(popDat$FAZ) == FALSE)
# sum(unique(popDat$FAZ) %in% fazNames == FALSE)

mazNames <- c("WHG", "NHG", "NSKEst", "HStr", "SFj", "WVI", "GStr")
# # Check they're all in; should be zero
# sum(mazNames %in% unique(popDat$MAZ) == FALSE)
# sum(unique(popDat$MAZ) %in% mazNames == FALSE)


#------------------------------------------------------------------------------
# Habitat variables
#------------------------------------------------------------------------------
habPressures <- popDat[, habNames]
# habPressures[which(habPressures < 0, arr.ind = TRUE)] <- 0