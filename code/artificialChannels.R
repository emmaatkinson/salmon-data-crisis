###############################################################################
# Major artifical channels listed by DFO SEP
# https://www.pac.dfo-mpo.gc.ca/sep-pmvs/projects-projets/index-eng.html#spawning
###############################################################################

#------------------------------------------------------------------------------
# Fulton River (1965, 1971)
# Exclusively sockeye
#------------------------------------------------------------------------------

pop <- unique(nusedsCU$POPULATION[grep("Fulton", nusedsCU$POPULATION)])
# There are Chinook, Coho, Pink, and Sockeye populations potentially affected

spPops <- list(Chinook = c(1:2), Coho = c(3:4), Pink = c(5:6), Sockeye = c(7:9))

for(i in 1:length(spPops)){
	p <- unique(nusedsCU$POP_ID[nusedsCU$POPULATION %in% pop[spPops[[i]]]])
	pp(p[1], xlim = c(1950, 2019), cex = 1.5, ylim = c(0, max(nusedsCU$MAX_ESTIMATE[nusedsCU$POP_ID %in% p])))
	
	pp(p[2], pch = 1, cex = 1.5, add = TRUE)
	if(length(p) == 2) legend("topleft", pch = c(19,1), paste(pop[spPops[[i]]], p, sep = "-"))
	
	if(length(p) == 3){
		pp(p[3], pch = 2, cex = 1.5, add = TRUE)
		legend("topleft", pch = c(19,1,2), pop[spPops[[i]]])
	}
	
	abline(v = c(1965, 1971), lty = 2)
}

#------------------------------------------------------------------------------
# Gates Creek Spawning Channel (1968)
#------------------------------------------------------------------------------
pop <- unique(nusedsCU$POPULATION[grep("Gates", nusedsCU$POPULATION)])
id <-  unique(nusedsCU$POP_ID[nusedsCU$POPULATION %in% pop])

G <- data.frame(year = c(1968:2019), channel = NA, creek = NA)
for(i in 1:dim(G)[1]){
	ind1 <- which(nusedsCU$POP_ID == 45102 & nusedsCU$ANALYSIS_YR == G$year[i])
	if(length(ind1) == 1)	G[i, 'channel'] <- nusedsCU$MAX_ESTIMATE[ind1]
	
	ind <- which(nusedsCU$POP_ID == 47211 & nusedsCU$ANALYSIS_YR == G$year[i])
	if(length(ind) == 1)	G[i, 'creek'] <- nusedsCU$MAX_ESTIMATE[ind]
}

# Mainly sockeye spawning channel; look at that population
p <- unique(nusedsCU$POP_ID[nusedsCU$POPULATION %in% pop[4:5]])
pp(p[2], xlim = c(1950, 2019), cex = 1.5, ylim = c(0, max(nusedsCU$MAX_ESTIMATE[nusedsCU$POP_ID %in% p])))
pp(p[1], pch = 1, cex = 1.5, add = TRUE)
legend("topleft", pch = c(1,19), paste(pop[4:5], p, sep = "-"))

par(new = TRUE)
plot(G$year, G$channel/G$creek, "l", col = 2, lwd = 2, xlim = c(1950, 2019), ylim = c(0,10), yaxt = "n", xaxt = "n")
axis(side = 4)
abline(h = 0.3, lty = 2, col = 2)

# Other species
par(mfrow = c(1,3))
id <-  unique(nusedsCU$POP_ID[which(nusedsCU$POPULATION %in% pop & nusedsCU$SPECIES == "Pink")])
pp(id)
abline(v = 1968)

#------------------------------------------------------------------------------
# Horsefly Channel (1989)
#------------------------------------------------------------------------------
pop <- unique(nusedsCU$POPULATION[grep("Horsefly", nusedsCU$POPULATION)])
# Focus on sockeye
pop <- pop[4:8]
id <-  unique(nusedsCU$POP_ID[nusedsCU$POPULATION %in% pop])
	
lims <- cbind(X = extendrange(cu$X_LONGT[cu$POP_ID %in%id], f = 0.5), 
							Y = extendrange(cu$Y_LAT[cu$POP_ID %in%id], f= 0.2))
waterClipped <- subset(waterHigh, waterHigh$X > lims[1,"X"] & waterHigh$X < lims[2,"X"] & waterHigh$Y > 52 & waterHigh$Y < 53)
river50clipped <- subset(river50, river50$X > lims[1,"X"] & river50$X < lims[2,"X"] & river50$Y > 52 & river50$Y < 53) # Is this more detailed?
	
plotMap(land, xlim = lims[, "X"], ylim = lims[, "Y"],	col = grey(0.8), bg = grey(0.8), las = 1, border = grey(0.8))
addPolys(river50clipped, col = "aliceblue", border = grey(0.6), colHoles = grey(0.8))

points(cu$X_LONGT[cu$POP_ID %in%id], cu$Y_LAT[cu$POP_ID %in%id], pch = 19)
points(cu$X_LONGT[cu$POP_ID %in%id[2]], cu$Y_LAT[cu$POP_ID %in%id[2]], pch = 2, col = 2, lwd = 2)
text(cu$X_LONGT[cu$POP_ID %in%id], cu$Y_LAT[cu$POP_ID %in%id], cu$SYSTEM_SITE[cu$POP_ID %in%id], cex = 0.8, pos = c(4,4,1,3,3), xpd = NA)

par(mfrow = c(3,2), mar = c(0,0,0,0), oma = c(5,5,2,1))
# plot(nusedsCU$ANALYSIS_YR[nusedsCU$POP_ID %in% id], nusedsCU$MAX_ESTIMATE[nusedsCU$POP_ID %in% id], "n", xlab = "", ylab = "MAX_ESTIMATE")
for(i in 1:5){
	plot(nusedsCU$ANALYSIS_YR[nusedsCU$POP_ID %in% id], nusedsCU$MAX_ESTIMATE[nusedsCU$POP_ID %in% id], "n", xlab = "", ylab = "MAX_ESTIMATE", xaxt = "n", yaxt = "n")
	mtext(side = 3, line = -2, unique(cu$SYSTEM_SITE[cu$POP_ID %in% id[i]]), cex = 0.8)
	points(nusedsCU$ANALYSIS_YR[nusedsCU$POP_ID %in% id[i]], nusedsCU$MAX_ESTIMATE[nusedsCU$POP_ID %in% id[i]], pch = c(19,2,19,19,19)[i], col = c(1,2,1,1,1)[i])
	if(i >3) axis(side = 1) else axis(side = 1, labels = FALSE)
	if(i %in% c(1,3,5)) axis(side = 2) else axis(side = 2, labels = FALSE)
}

#------------------------------------------------------------------------------
# Nadina River Spawning Channel (1973) 
# Exclusively sockeye
#------------------------------------------------------------------------------
pop <- unique(nusedsCU$POPULATION[which(grepl("Nadina", nusedsCU$POPULATION) & nusedsCU$SPECIES == "Sockeye")])
id <-  unique(nusedsCU$POP_ID[nusedsCU$POPULATION %in% pop])

pp(id[1], pch = 19, xlim = c(1950, 2019))
pp(id[2], pch = 1, cex = 1.5, add = TRUE)
abline(v = 1973, lty = 2)

#------------------------------------------------------------------------------
# Pinkut Creek Spawning Channel (1968)
# Exclusively sockeye
#------------------------------------------------------------------------------
pop <- unique(nusedsCU$POPULATION[which(grepl("Pinkut", nusedsCU$POPULATION) & nusedsCU$SPECIES == "Sockeye")])
id <-  unique(nusedsCU$POP_ID[nusedsCU$POPULATION %in% pop])

# Map it - not very useful. ALl close together
lims <- cbind(X = c(-122.2, -121.5), 
							Y = c(49.1, 49.5))
river50clipped <- subset(river50, river50$X > lims[1, "X"] & river50$X < lims[2, "X"] & river50$Y > lims[1,"Y"] & river50$Y < lims[2,"Y"]) # Is this more detailed?

plotMap(land, xlim = lims[, "X"], ylim = lims[, "Y"],	col = grey(0.8), bg = grey(0.8), las = 1, border = grey(0.8))
addPolys(river50clipped, col = "aliceblue", border = grey(0.6), colHoles = grey(0.8))

points(cu$X_LONGT[cu$POP_ID %in%id], cu$Y_LAT[cu$POP_ID %in%id], pch = 19)
points(cu$X_LONGT[cu$POP_ID %in%id[2]], cu$Y_LAT[cu$POP_ID %in%id[2]], pch = 2, col = 2, lwd = 2)
text(cu$X_LONGT[cu$POP_ID %in%id], cu$Y_LAT[cu$POP_ID %in%id], cu$SYSTEM_SITE[cu$POP_ID %in%id], cex = 0.8, pos = c(4,4,1,3,3), xpd = NA)

# Population estimates

pp(id[1], pch = 19, xlim = c(1950, 2019), ylim = c(0, 250000))
pp(id[2], pch = 1, cex = 1.5, add = TRUE)
pp(id[3], pch = 2, cex = 1.5, add = TRUE)
abline(v = 1968, lty = 2)
legend("topleft", pch = c(19,1,2), pt.cex = c(1,1.5, 1.5), pop)

# Classification
cu$GFE_TYPE[cu$POP_ID %in% id]

#------------------------------------------------------------------------------
# Weaver Creek Spawning Channel (1965)
# Mostly sockeye, but in smaller numbers chum and pink
#------------------------------------------------------------------------------
pop <- unique(nusedsCU$POPULATION[grep("Weaver", nusedsCU$POPULATION)]) 
id <-  unique(nusedsCU$POP_ID[nusedsCU$POPULATION %in% pop])

# Map it - not very useful. ALl close together
lims <- cbind(X = c(-125.6, -125.3), 
							Y = c(54.4, 54.5))
river50clipped <- subset(river50, river50$X > -126 & river50$X < -125 & river50$Y > lims[1,"Y"] & river50$Y < lims[2,"Y"]) # Is this more detailed?

plotMap(land, xlim = lims[, "X"], ylim = lims[, "Y"],	col = grey(0.8), bg = grey(0.8), las = 1, border = grey(0.8))
addPolys(river50clipped, col = "aliceblue", border = grey(0.6), colHoles = grey(0.8))

points(cu$X_LONGT[cu$POP_ID %in%id], cu$Y_LAT[cu$POP_ID %in%id], pch = 19)
points(cu$X_LONGT[cu$POP_ID %in%id[6:7]], cu$Y_LAT[cu$POP_ID %in%id[6:7]], pch = 2, col = 2, lwd = 2)

# Population estimates
par(mfrow = c(3,2), mar = rep(0, 4), oma = c(5,5,2,1))
for(i in 1:5){
	pop.i <- unique(nusedsCU$POPULATION[which(grepl("Weaver", nusedsCU$POPULATION) & nusedsCU$SPECIES == c("Chinook", "Coho", "Chum", "Pink", "Sockeye")[i])]) 
	id.i <-  unique(nusedsCU$POP_ID[nusedsCU$POPULATION %in% pop.i])
	
	plot(nusedsCU$ANALYSIS_YR[nusedsCU$POP_ID %in% id], nusedsCU$MAX_ESTIMATE[nusedsCU$POP_ID %in% id], "n", xlab = "", ylab = "MAX_ESTIMATE", xaxt = "n", yaxt = "n", ylim = c(0, 100000))
	mtext(side = 3, line = -2, c("Chinook", "Coho", "Chum", "Pink", "Sockeye")[i], cex = 0.8)
	points(nusedsCU$ANALYSIS_YR[nusedsCU$POP_ID %in% id.i[1]], nusedsCU$MAX_ESTIMATE[nusedsCU$POP_ID %in% id.i[1]], pch = 1)
	if(length(id.i) >1){
		for(j in 2:length(id.i)){
			points(nusedsCU$ANALYSIS_YR[nusedsCU$POP_ID %in% id.i[j]], nusedsCU$MAX_ESTIMATE[nusedsCU$POP_ID %in% id.i[j]], pch = j)
		}}
	if(i >3) axis(side = 1) else axis(side = 1, labels = FALSE)
	if(i %in% c(1,3,5)) axis(side = 2) else axis(side = 2, labels = FALSE)
	abline(v = 1965, lty = 2)
}
legend(2030, 800000, pch = c(1,2))

###############################################################################
# Artificial channels and the river populations they affect
###############################################################################

# Load NuSEDS CU data steps 1 - 2 
# Removing NA years and truncating data


cu <- read.csv("data/conservation_unit_system_sites.csv")

# Then remove the non-artifical channels that have < 10 years of data
# since these will be removed anyway
nYr <- tapply(nusedsCU$SQ_POP_ID, nusedsCU$SQ_POP_ID, length)
hist(nYr)

SQ_POP_ID10 <- names(nYr[which(nYr >= 10)])
POP_ID10 <- unique(nusedsCU$POP_ID[match(SQ_POP_ID10, nusedsCU$SQ_POP_ID)])

rmPop <- unique(cu$POP_ID[(cu$GFE_TYPE != "Artificial Channel") & ((cu$POP_ID %in% POP_ID10) == FALSE)])
nusedsCU <- nusedsCU[ - which(nusedsCU$POP_ID %in% rmPop), ]
popCount()

#-----
cu <- cu[cu$POP_ID %in% unique(nusedsCU$POP_ID), ]

# Dataframe of artificial channels
ac <- data.frame(
	POP_ID = unique(cu$POP_ID[cu$GFE_TYPE == "Artificial Channel"]),
	meanMAX_ESTIMATE = NA,
	nYears = NA,
	minYear = NA,
	maxYear = NA,
	POPULATION = NA
)
n <- dim(ac)[1]

for(i in 1:n){
	y <- nusedsCU$ANALYSIS_YR[nusedsCU$POP_ID == ac$POP_ID[i]]
	ac$meanMAX_ESTIMATE[i] <- mean(nusedsCU$MAX_ESTIMATE[nusedsCU$POP_ID == ac$POP_ID[i]])
	ac$nYears[i] <- length(y)
	ac$minYear[i] <- min(y)
	ac$maxYear[i] <- max(y)
	ac$POPULATION[i] <- nusedsCU$POPULATION[nusedsCU$POP_ID == ac$POP_ID[i]][1]
}

ac1 <- ac[which(ac$meanMAX_ESTIMATE > 1000 & ac$nYears >= 5), ]
channelPOP_ID <- ac1$POP_ID

# Load mapping data
library(PBSmapping)
gshhg <- "~/Google Drive/Mapping/gshhg-bin-2.3.7/"
xlim <- c(-141, -118) + 360
ylim <- c(48, 67)
land <- importGSHHS(paste0(gshhg,"gshhs_i.b"), xlim = xlim, ylim = ylim, maxLevel = 2, useWest = TRUE)
# waterHigh <- importShapefile("data/ignore/FWAtlas_water/NTS_BC_RIV_LAKE_WET_POLYS_125M/250_WATPS_polygon.shp", projection = "LL") #1:250,000
waterHigh <- importShapefile("data/ignore/FWAtlas_water/BCGW_7113060B_1618544330829_17396/WDIC_WATERBODY_POLY_SVW/WSA_WB_PLY_polygon.shp", projection = "LL") #1:50,000

# For each population, map artificial channel (red star) and all other populations in that CU
# (black circles)
# Add ring for 10 km radius around artificial channel, and identify other populations within 
# 10 km (as the crow flies) in the same river system
plotCircle <- function(
	LonDec, #LonDec = longitude in decimal degrees
	LatDec, #LatDec = latitude in decimal degrees of the center of the circle
	Km, #Km = radius of the circle in kilometers
	plot = TRUE # Plot the circle or return coordinates?
) {
	ER <- 6371 #Mean Earth radius in kilometers. C
	AngDeg <- seq(1:360) #angles in degrees 
	Lat1Rad <- LatDec*(pi/180) #Latitude of the center of the circle in radians
	Lon1Rad <- LonDec*(pi/180) #Longitude of the center of the circle in radians
	AngRad <- AngDeg*(pi/180) #angles in radians
	Lat2Rad <-asin(sin(Lat1Rad)*cos(Km/ER)+cos(Lat1Rad)*sin(Km/ER)*cos(AngRad)) #Latitude of each point of the circle rearding to angle in radians
	Lon2Rad <- Lon1Rad+atan2(sin(AngRad)*sin(Km/ER)*cos(Lat1Rad),cos(Km/ER)-sin(Lat1Rad)*sin(Lat2Rad))#Longitude of each point of the circle rearding to angle in radians
	Lat2Deg <- Lat2Rad*(180/pi)#Latitude of each point of the circle rearding to angle in degrees (conversion of radians to degrees deg = rad*(180/pi) )
	Lon2Deg <- Lon2Rad*(180/pi)#Longitude of each point of the circle rearding to angle in degrees (conversion of radians to degrees deg = rad*(180/pi) )
	if(plot == FALSE){
		return(cbind(X = Lon2Deg, Y = Lat2Deg))
	} else {
		polygon(Lon2Deg, Lat2Deg,lty=2)
	}
}


closePopulations <- data.frame(
	channelPOP_ID = NA,
	riverPOP_ID = NA
)

###############################################################################
i <- 9
# for(i in 1:length(channelPOP_ID)){

	lims <- apply(plotCircle(
		LonDec = cu$X_LONGT[cu$POP_ID == channelPOP_ID[i]], 
		LatDec =  cu$Y_LAT[cu$POP_ID == channelPOP_ID[i]],
		Km = 20,
		plot = FALSE), 2, range)
	
	circle <- plotCircle(
		LonDec = cu$X_LONGT[cu$POP_ID == channelPOP_ID[i]], 
		LatDec =  cu$Y_LAT[cu$POP_ID == channelPOP_ID[i]],
		Km = 10,
		plot = FALSE)
	circle <- as.PolySet(x = data.frame(PID = rep(1, 360), POS = c(1:360), X = circle[,"X"], Y = circle[, "Y"]), projection = "LL")
	
	waterClipped <- subset(waterHigh, waterHigh$X > lims[1,"X"] & waterHigh$X < lims[2,"X"] & waterHigh$Y > lims[1,"Y"] & waterHigh$Y < lims[2, "Y"])
	
	cu.i <- cu[cu$FULL_CU_IN == cu$FULL_CU_IN[which(cu$POP_ID == channelPOP_ID[i])], ]
	
	# Which populations in that CU would be susceptible to influence?
	# I.e., have an abundance that is less than channel/0.3?
	# channel/river > 0.3 then remove
	# so only look at pops that have river < channel/0.3
	nuseds.i <- nusedsCU[which(nusedsCU$CU_NAME == cu.i$CU_NAME[1]), ]
	meanMAX_ESTIMATE <- tapply(nuseds.i$MAX_ESTIMATE, nuseds.i$POP_ID, mean)
	channelMean <- meanMAX_ESTIMATE[which(names(meanMAX_ESTIMATE) == channelPOP_ID[i])]
	POP_ID <- names(which((meanMAX_ESTIMATE < (channelMean/0.3)) == TRUE))
	
	# nonChannels <- data.frame(
	# 	EID = c(1:length(which(cu.i$GFE_TYPE != "Artificial Channel"))), 
	# 	X = cu.i$X_LONGT[which(cu.i$GFE_TYPE != "Artificial Channel")],
	# 	Y = cu.i$Y_LAT[which(cu.i$GFE_TYPE != "Artificial Channel")],
	# 	POP_ID = cu.i$POP_ID[which(cu.i$GFE_TYPE != "Artificial Channel")]) 
	ind <- which((cu.i$GFE_TYPE != "Artificial Channel") & (cu.i$POP_ID %in% POP_ID))
	nonChannels <- data.frame(
		EID = c(1:length(ind)),
		X = cu.i$X_LONGT[ind],
		Y = cu.i$Y_LAT[ind],
		POP_ID = cu.i$POP_ID[ind])
	
	close2channel <- findPolys(events = nonChannels, polys = circle)
	
	# if(length(close2channel) > 0){
		
		# Inset map
		plotMap(land, xlim = c(-135, -118), ylim = c(48.17, 58),	col = grey(0.8), bg = "aliceblue", las = 1, border = grey(0.6), lwd = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
		polygon(x = lims[c(1,2,2,1), 'X'], y = lims[c(1,1,2,2), 'Y'], col = NA, lwd = 1.5)
		quartz.save(file = paste0("figures/artificialChannels/", cu$SYSTEM_SITE[cu$POP_ID == channelPOP_ID[i]], channelPOP_ID[i], "_inset.png"), type = "png", device = dev.cur())
		
		# Channel map
		plotMap(land, xlim = lims[, "X"], ylim = lims[, "Y"],	col = grey(0.8), bg = "aliceblue", las = 1, border = grey(0.6))
		mtext(side = 3, nusedsCU$POPULATION[nusedsCU$POP_ID == channelPOP_ID[i]][1], line = 2)
		mtext(side = 3, cu$FULL_CU_IN[cu$POP_ID == channelPOP_ID[i]][1], line = 0.5)
		
		addPolys(waterClipped, col = "aliceblue", border = grey(0.6), colHoles = grey(0.8))
		
		addPoints(nonChannels)
		polygon(circle$X, circle$Y, lty = 2)
		addPoints(nonChannels[nonChannels$EID %in% close2channel$EID, ], pch = 19)
		points(cu.i$X_LONGT[which(cu.i$POP_ID == channelPOP_ID[i])], cu.i$Y_LAT[which(cu.i$POP_ID == channelPOP_ID[i])], pch = 2, col = 2, cex = 1, lwd = 2)
		
		systemPoly <- locatePolys(type = "l", col = 4)
		close2channelSameSystem <- findPolys(events = nonChannels[nonChannels$EID %in% close2channel$EID, ], polys = systemPoly)
		
		#----------------------------------------------------------------------------
		# if(length(close2channelSameSystem) > 0){
			addPoints(nonChannels[nonChannels$EID %in% close2channelSameSystem$EID, ], cex = 2, col = 4)
		
		#------------------
		zz <- nonChannels$POP_ID[nonChannels$EID %in% close2channelSameSystem$EID]
		
		closePopulations <- rbind(closePopulations, cbind(channelPOP_ID = rep(channelPOP_ID[i], length(zz)), riverPOP_ID = zz))
		
		quartz.save(file = paste0("figures/artificialChannels/", cu$SYSTEM_SITE[cu$POP_ID == channelPOP_ID[i]], channelPOP_ID[i], ".png"), type = "png", device = dev.cur())
	
		# Plot populations
		popNames <- c(unique(nusedsCU$POPULATION[nusedsCU$POP_ID == channelPOP_ID[i]]))
		
		pp(channelPOP_ID[i], col = 2, xlim = c(1950, 2019), ylim = c(0, max(nusedsCU$MAX_ESTIMATE[nusedsCU$POP_ID %in% c(channelPOP_ID[i], zz)])))
		for(j in 1:length(zz)){
			pp(zz[j], pch = j + 1,add = TRUE, col = 1)
			popNames <- c(popNames, unique(nusedsCU$POPULATION[nusedsCU$POP_ID == zz[j]]))
		}
		legend("topleft", pch = c(19, 1:length(zz)), col = c(2, rep(1, length(zz))), legend = popNames)
		quartz.save(file = paste0("figures/artificialChannels/", cu$SYSTEM_SITE[cu$POP_ID == channelPOP_ID[i]], channelPOP_ID[i], "_spawners.png"), type = "png", device = dev.cur())
		
		print(c(unique(nusedsCU$POPULATION[nusedsCU$POP_ID == channelPOP_ID[i]])))
		
	# # If there were no populations wihtin a 10 km radius	
	# } else {
	# 	closePopulations <- rbind(closePopulations,
	# 														cbind(channelPOP_ID = channelPOP_ID[i], riverPOP_ID = "none on same river/stream"))
	# }} else {
	# 	closePopulations <- rbind(closePopulations,
	# 														cbind(channelPOP_ID = channelPOP_ID[i], riverPOP_ID = "none in radius"))
	# }
	
	# Periodically save
	write.csv(closePopulations, file = "data/ignore/artificialChannels/closePopulations.csv")
	
# }
