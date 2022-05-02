library(maptools)
library(rgdal)
library(sp)
library(mapdata)
library(ggplot2)
library(PNWColors)
#------------------------------------------------------------------------------
# Load spatial data
#------------------------------------------------------------------------------
col2 <- pnw_palette("Bay", n = 2)
# Load mapping data
library(PBSmapping)
gshhg <- "~/Google Drive/Mapping/gshhg-bin-2.3.7/"
xlim <- c(-134, -118)
ylim <- c(48.17, 57.5)
# five resolutions: crude(c), low(l), intermediate(i), high(h), and full(f).
res <- "i"
land <- importGSHHS(paste0(gshhg,"gshhs_", res, ".b"), xlim = xlim + 360, ylim = ylim, maxLevel = 2, useWest = TRUE)
rivers <- importGSHHS(paste0(gshhg,"wdb_rivers_", res, ".b"), xlim = xlim  + 360, ylim = ylim, useWest = TRUE)
borders <- importGSHHS(paste0(gshhg,"wdb_borders_", res, ".b"), xlim = xlim  + 360, ylim = ylim, useWest = TRUE, maxLevel = 1)

#------------------------------------------------------------------------------
# Import faz and convert to LL
#------------------------------------------------------------------------------

# faz <- importShapefile("~/Google Drive/DFO/DFO_popHab/data/ignore/PSF/faz/FreshwaterAdaptiveZones.shp")
# faz_sp <- readOGR(dsn = "~/Google Drive/DFO/DFO_popHab/data/ignore/PSF/faz/", layer = "FreshwaterAdaptiveZones")
# 
# crslatlong <- CRSargs(CRS("+init=epsg:4326"))
# fazLL_sp <- spTransform(faz_sp, CRS(crslatlong))  
# 
# # Convert to Polyset
# faz$FAZ_Name <- NA
# faz$FAZ_Acrony <- NA
# faz$FAZ_Code <- NA
# 
# for(i in 1:33){
# 	for(j in unique(faz$SID[which(faz$PID == i)])){
# 		LLcoords.ij <- fazLL_sp@polygons[[i]]@Polygons[[j]]@coords
# 		if(dim(LLcoords.ij)[1] == length(which(faz$PID == i & faz$SID == j))){
# 			faz[which(faz$PID == i & faz$SID == j), c('X', 'Y')] <- LLcoords.ij
# 			faz[which(faz$PID == i & faz$SID == j), c("FAZ_Name", "FAZ_Acrony", "FAZ_Code")] <- fazLL_sp@data[i, c("FAZ_Name", "FAZ_Acrony", "FAZ_Code")]
# 		} else {
# 			stop("Stop: number of points don't match.")
# 		}
# }}
# 
# write.csv(faz, file = "~/Google Drive/DFO/DFO_popHab/data/ignore/PSF/fazLL.csv")
# 
# faz <- as.PolySet(read.csv("~/Google Drive/DFO/DFO_popHab/data/ignore/PSF/fazLL.csv"), projection = "LL")
# 
# 
# dim(faz_l)
# dim(faz)
# 
# faz$FAZ_Acrony[faz$FAZ_Acrony == "QCI"] <- "HG"
# unique(faz$FAZ_Name)
# faz$FAZ_Name[faz$FAZ_Name == "Queen Charlottes"] <- "Haida Gwaii"
# 
# unique(popDat$FAZ) %in% sort(unique(faz$FAZ_Acrony))
# 
# faz$inDat <- faz$FAZ_Acrony %in% unique(popDat$FAZ) 
# 
# pidIn <- unique(faz$PID[faz$inDat == TRUE])
# 
# faz_l <- thinPolys(faz, tol = 1) #thin faz borders to 1 km
# 

# faz_l <- faz_l[faz_l$PID %in% pidIn,]
# faz_l$FAZ_Name <- faz$FAZ_Name[match(faz_l$PID, faz$PID)] #c("FAZ_Name", "FAZ_Acrony", "FAZ_Code")
# faz_l$FAZ_Acrony <- faz$FAZ_Acrony[match(faz_l$PID, faz$PID)]
# faz_l$FAZ_Code <- faz$FAZ_Code[match(faz_l$PID, faz$PID)]
# 
# write.csv(faz_l, file = "~/Google Drive/DFO/DFO_popHab/data/ignore/PSF/fazLL_thinned_inDat.csv")

faz <- as.PolySet(read.csv("data/ignore/PSF/fazLL_thinned_inDat.csv"), projection = "LL")

#------------------------------------------------------------------------------
# Import MAZ and convert to LL
#------------------------------------------------------------------------------

# maz <- importShapefile("~/Google Drive/DFO/DFO_popHab/data/ignore/PSF/MarineAdaptiveZones/MarineAdaptiveZones.shp")
# maz_sp <- readOGR(dsn = "~/Google Drive/DFO/DFO_popHab/data/ignore/PSF/MarineAdaptiveZones/", layer = "MarineAdaptiveZones")
# 
# crslatlong <- CRSargs(CRS("+init=epsg:4326"))
# mazLL_sp <- spTransform(maz_sp, CRS(crslatlong))
# 
# # Convert to Polyset
# maz$MAZ_name <- NA
# maz$MAZ_acrony <- NA
# maz$MAZ_code <- NA
# 
# inclMAZ <- sort(c(12, 6, 7, 5, 10, 13, 4))
# for(i in c(6,7,10,12,13)){
# 	for(j in unique(maz$SID[which(maz$PID == i)])){
# 		LLcoords.ij <- mazLL_sp@polygons[[i]]@Polygons[[j]]@coords
# 		if(dim(LLcoords.ij)[1] == length(which(maz$PID == i & maz$SID == j))){
# 			maz[which(maz$PID == i & maz$SID == j), c('X', 'Y')] <- LLcoords.ij
# 			maz[which(maz$PID == i & maz$SID == j), c("MAZ_name", "MAZ_acrony", "MAZ_code")] <- mazLL_sp@data[i, c("MAZ_name", "MAZ_acrony", "MAZ_code")]
# 		} else {
# 			stop("Stop: number of points don't match.")
# 		}
# 	}
# 	print(paste0("finished ", i))}
# 
# write.csv(maz, file = "~/Google Drive/DFO/DFO_popHab/data/ignore/PSF/mazLL.csv")
# 
# #------------------------------------------------------------------------------
# # maz <- as.PolySet(read.csv("~/Google Drive/DFO/DFO_popHab/data/ignore/PSF/mazLL.csv"), projection = "LL")
# maz <- as.PolySet(maz, projection = "LL")
# 
# # Remove MAZs not in our analysis
# unique(popDat$MAZ) %in% sort(unique(maz$MAZ_acrony))
# maz$inDat <- maz$MAZ_acrony %in% unique(popDat$MAZ)
# pidIn <- unique(maz$PID[maz$inDat == TRUE])
# maz <- maz[maz$PID %in% pidIn,]
# 
# # Change Queen Charlotte Islands acrony 
# maz$MAZ_acrony[maz$MAZ_acrony == "NQCI"] <- "NHG"
# maz$MAZ_acrony[maz$MAZ_acrony == "WQCI"] <- "WHG"
# unique(maz$MAZ_name)
# 
# # thin Polygon to 1 km resolution
# maz_l <- thinPolys(maz, tol = 1) #thin maz borders to 1 km
# 
# # Select only MAZs in our analysis and change to capital Acrony, Name etc. to be
# # consistent with FAZ poly naming
# # maz_l <- maz_l[maz_l$PID %in% pidIn,]
# maz_l$MAZ_Name <- maz$MAZ_name[match(maz_l$PID, maz$PID)] 
# maz_l$MAZ_Acrony <- maz$MAZ_acrony[match(maz_l$PID, maz$PID)]
# maz_l$MAZ_Code <- maz$MAZ_code[match(maz_l$PID, maz$PID)]
# 
# write.csv(maz_l, file = "~/Google Drive/DFO/DFO_popHab/data/ignore/PSF/mazLL_thinned_inDat.csv")

maz <- as.PolySet(read.csv("data/ignore/PSF/mazLL_thinned_inDat.csv"), projection = "LL")

# take out maz points that directly overlap shoreline data
head(maz)
head(land)
# No direct matches to perhaps round?
ind <- which(((round(maz$X, 2) %in% round(land$X, 2)) & (round(maz$Y, 3) %in% round(land$Y, 3)))==FALSE)

mazInterior <- maz[ind, ]

#------------------------------------------------------------------------------
# Import population locations
#------------------------------------------------------------------------------

popDat <- read.csv("data/popDat.csv")
popDat <- popDat[popDat$forestDD == 0, ]

cu <- read.csv("data/conservation_unit_system_sites.csv")

cu <- cu[cu$POP_ID %in% popDat$POP_ID, ]
length(unique(cu$POP_ID))
length(unique(popDat$POP_ID))
# Almost all of them in there

popDat$X_LONGT <- cu$X_LONGT[match(popDat$POP_ID, cu$POP_ID)]
popDat$Y_LAT <- cu$Y_LAT[match(popDat$POP_ID, cu$POP_ID)]

popDat[is.na(popDat$X_LONGT), ]

range(popDat$X_LONGT, na.rm = TRUE)
range(popDat$Y_LAT, na.rm = TRUE)

uniqueLoc <- sort(unique(paste(popDat$Y_LAT, popDat$X_LONGT)), decreasing = TRUE)
indStream <- numeric(length(uniqueLoc))
for(i in 1:length(uniqueLoc)){
	indStream[i] <- which(paste(popDat$Y_LAT, popDat$X_LONGT) == uniqueLoc[i])[1]
}

indStream <- indStream[2:length(indStream)] # Remove first with blank LL

streamDat <- as.EventData(data.frame(
	EID = c(1:length(indStream)),
	X = popDat$X_LONGT[indStream],
	Y = popDat$Y_LAT[indStream],
	POPULATION = popDat$POPULATION[indStream],
	FAZ = popDat$FAZ[indStream],
	MAZ = popDat$MAZ[indStream]
))


#------------------------------------------------------------------------------
# Plot maps of faz and maz borders
#------------------------------------------------------------------------------
polyCol <- pnw_palette("Bay", n = 2)

fazCentroid <- calcCentroid(faz, rollup = 1)
fazCentroid$FAZ_Acrony <- faz$FAZ_Acrony[match(fazCentroid$PID, faz$PID)]
fazCentroid[fazCentroid$FAZ_Acrony == "EVI",c("X","Y")] <- c(-126.7291, 50.22792)
fazCentroid[fazCentroid$FAZ_Acrony == "HecLow",c("X","Y")] <- c(-126.7291, 50.22792)


mazCentroid <- calcCentroid(maz, rollup = 1)
mazCentroid$MAZ_Acrony <- maz$MAZ_Acrony[match(mazCentroid$PID, maz$PID)]
mazCentroid[mazCentroid$MAZ_Acrony == "WHG",c("X","Y")] <- c(-129.2388, 50.41241)

fazCol <- pnw_palette("Bay", n = 22)[c(12,10,6,4,15,17,5,20,11,8, 13,1,7,21,18,22,3,16,14,19,2 ,9)]
mazCol <- pnw_palette("Bay", n = 8)[c(1,8,3,7,5,4,6,2)]


# quartz(width = 3.54, height = 6, pointsize = 10)
# png(file = "figures/faz_maz.png", width = 1400, height = 750, res= 150)

quartz(width = 4, height = 4, pointsize = 10)

# par(mfrow = c(2,1), mar = c(2,3,2,1), oma= c(4, 0, 0, 0))

# MAZ map
plotMap(land, col = "white", bg = "#0f85a050", las = 1, border = "#0f85a0", lwd = 0.6, xlab = "", ylab = "", ylim = c(48.2, 57.2))

addLines(rivers, col = col2[1], lwd = 0.5)
# addLines(borders, lwd = 1.5)

#addPolys(maz, col = paste0(mazCol[match(mazCentroid$MAZ_Acrony, mazNames)], "80"), lwd = 0.5, border = mazCol[match(mazCentroid$MAZ_Acrony, mazNames)]) # Colourful FAZs
addPolys(maz, col = "#00000000", border = 1, lwd = 1)
addPoints(streamDat, pch = 19, cex = 0.5, col = "#dd412460")

addLines(mazInterior, col = 3)

addPolys(maz, col = "#00000000", border = 1, lwd = 0.5)

text(mazCentroid$X, mazCentroid$Y, mazCentroid$MAZ_Acrony, font = 2, cex = 0.8, col = 1) #col = fazCol[match(fazNames, fazCentroid$FAZ_Acrony)],
mtext(side = 3, line = 1, "a) Marine Adaptive Zones (MAZs)", adj = 0)

# FAZ map
plotMap(land, ylim = c(48.2, 57.2),	col = "white", bg = grey(0.9), las = 1, border = grey(0.8), lwd = 0.6, xlab = "", ylab = "")
addLines(rivers, col = grey(0.8), lwd = 0.5)
addLines(borders, lwd = 1.5)
mtext(side = 1, outer= TRUE, expression(paste(degree, "Longitude")), line = 1, cex = 0.8)
#addPolys(faz, col = paste0(fazCol[match(fazCentroid$FAZ_Acrony, fazNames)], "80"), lwd = 0.5, border = fazCol[match(fazCentroid$FAZ_Acrony, fazNames)]) # Colourful FAZs
addPolys(faz, col = "#00000000", border = grey(0.4), lwd = 0.5)

text(fazCentroid$X, fazCentroid$Y, fazCentroid$FAZ_Acrony, font = 2, cex = 0.8)#, col = 

mtext(side = 3, line = 1, "b) Freshwater Adaptive Zones (FAZs)", adj = 0)
mtext(side = 2, outer = TRUE, expression(paste(degree, "Latitude")))


# dev.off()

# # Colourful full faz names (old fig)
# png(file = "figures/fazNames.png", width = 900, height = 750, res= 150)
# plot(1,1,"n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", ylim = c(0.5, nFAZ+0.5))
# for(i in 1:nFAZ){
# 	text(1, nFAZ - i + 1, sort(unique(paste(faz$FAZ_Acrony, faz$FAZ_Name, sep = ": ")))[match(fazNames, sort(fazCentroid$FAZ_Acrony))][i], adj = 0, col = fazCol[i])
# }
# dev.off()


#----------


# MAZ map
plotMap(land, col = "white", bg = "#0f85a050", las = 1, border = "#0f85a0", lwd = 0.6, xlab = "", ylab = "", ylim = c(48.2, 57.2))
addLines(rivers, col = col2[1], lwd = 0.5)
addPolys(maz, col = "#00000000", border = 1, lwd = 1)
addPoints(streamDat, pch = 19, cex = 0.5, col = "#dd412460")

# FAZ map
plotMap(land, col = "white", bg = "#0f85a050", las = 1, border = "#0f85a0", lwd = 0.6, xlab = "", ylab = "", ylim = c(48.2, 57.2))
addLines(rivers, col = col2[1], lwd = 0.5)
addPolys(faz, col = "#00000000", border = 1, lwd = 1)
addPoints(streamDat, pch = 19, cex = 0.5, col = "#dd412460")

legend("topright", bg = "white", pch = 19, pt.cex = 0.5, col = "#dd412460", legend = "Included population")
