##############################################################################
###############################################################################
#
# Assembling spawner abundance data to visualise trends in monitoring
# coverage and data quality for wild Pacific salmon in BC
# Emma Atkinson <emmamargaretatkinson@gmail.com>
#
# Generating monitoring coverage figures
#
# Data were downloaded from Open Data Canada
# https://open.canada.ca/data/en/dataset/c48669a3-045b-400d-b730-48aafe8c5ee6
#
##############################################################################
##############################################################################

### Input ###
rm(list=ls())
graphics.off()

library(here)
setwd(here())

source("code/functions.R")

#------------------------------------------------------------------------------
# Read in NuSEDS data
#------------------------------------------------------------------------------

### Read in cleaned up data ###
nuseds_raw <- read.csv("data_EA/NuSEDS_2023-June.csv")
nuseds <- read.csv("data_EA/Conservation_Unit_Data_filtered_20230627.csv")
popdat <- read.csv("data_EA/popDat_20230627.csv")

popdat <- popdat[which(is.na(popdat$region)==FALSE),]

# Reading in collated spawner-recruit data by species #
nuseds$GEOGRAPHICAL_EXTNT_OF_ESTIMATE <- as.character(nuseds$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)

# Set up dataframe for summary numbers #
data <- nuseds

all.rivers <- unique(data$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
n.rivers <- length(all.rivers)

pink.rivers <- unique(data[data$SPECIES=="Pink",]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
n.pink <- length(pink.rivers)

podd.rivers <- unique(data[data$SPECIES=="Pink" & data$ANALYSIS_YR%%2==1,]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
n.podd <- length(pink.rivers)

peven.rivers <- unique(data[data$SPECIES=="Pink" & data$ANALYSIS_YR%%2==0,]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
n.peven <- length(pink.rivers)

chum.rivers <- unique(data[data$SPECIES=="Chum",]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
n.chum <- length(chum.rivers)

coho.rivers <- unique(data[data$SPECIES=="Coho",]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
n.coho <- length(coho.rivers)

chinook.rivers <- unique(data[data$SPECIES=="Chinook",]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
n.chinook <- length(chinook.rivers)

sock.rivers <- unique(data[data$SPECIES=="Sockeye",]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
n.sock <- length(sock.rivers)

# Calculating annual monitoring coverage #
ny <- length(1950:2021)

m <- data.frame(year = rep(c(1950:2021), 6),
								species = rep(c("Pink","Chum", "Coho", "Chinook", "Sockeye","all"),each=ny),
								odd_even = c(rep(c("O","E"), length.out=ny), rep(NA, 5*ny)),
								total.systems = c(rep(n.podd, ny), rep(c(n.chum, n.coho, n.chinook, n.sock, n.rivers), each=ny)),
								no.counted = rep(NA, 6*ny),
								prop.counted = rep(NA, 6*ny))

# Population data frame with monitoring coverage (number of systems/proportion of systems counted) #
for (species in unique(m$species)) {
	if(species=="all") {
		for (i in 1950:2021){
			m[m$species==species & m$year==i,]$no.counted = length(unique(data[data$ANALYSIS_YR==i,]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE))
		}
	}
	else { 
		for (i in 1950:2021){
			m[m$species==species & m$year==i,]$no.counted = length(unique(data[data$SPECIES==species & data$ANALYSIS_YR==i,]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)) }
	}
}


m$prop.counted = m$no.counted/m$total.systems

# Write to file #
# dir.table <- "C:/Users/Emma Atkinson/Desktop/Research/SoS/Tables"
# setwd(dir.table)
# write.csv(m, "monitoring_coverage_table_Feb-2019.csv")

# Plotting #
pdf("figures/River-level_monitoring_coverage_coastwide_Jun-2023.pdf", 6,6)
#tiff("monitoring_coverage_Sept-2019.tiff", width=4, height=4, units="in", pointsize=9, res=600)
#windows(9,9)

#showtext_auto()

par(mfrow=c(1,1))#, family="Avenir")#, mgp=c(2.5, 0.8, 0), mar=c(3.5, 3.5, 1, 2))

plot(unique(m$year), m[m$species=="all",]$no.counted, "l", las =1, bty="l", xlab = "Year", ylab="Number of systems monitored by DFO", lwd=1, ylim=c(0, 1500), col=grey(0.8), main="Coastwide")
# axis(1, at=seq(1950,2020,5), lab=seq(1950,2020,5))
# axis(2, at=seq(0,70,10), lab=seq(0,1,.1), las=2)

points(2021, 530, pch = 19, col="#2E8B57")
text(2021, 530, "n=530", font = 2, cex=0.8, xpd=NA, pos=4, col="#2E8B57")

lines(unique(m$year), m[m$species=="Chum",]$no.counted, col="grey20", lty=2, lwd=2)
lines(unique(m$year), m[m$species=="Coho",]$no.counted, col="black", lty=3, lwd=1)
lines(unique(m$year), m[m$species=="Chinook",]$no.counted, col="grey40", lty=4, lwd=1)
lines(unique(m$year), m[m$species=="Sockeye",]$no.counted, col="grey50", lty=5, lwd=1)
lines(unique(m$year), m[m$species=="Pink",]$no.counted, col="black", lwd=2)
lines(unique(m$year), m[m$species=="all",]$no.counted, lwd=4, col="#2E8B57")

legend("topright", lty= c(1,1,2,3,4,5), lwd=c(4,2,1,1,1,1), col=c("#2E8B57", "black","grey20","black","grey40","grey50"), 
			 c("All systems (n=80)", "Pink", "Chum","Coho","Chinook","Sockeye"), bty="n", cex=.7)


# plot(unique(m$year), m[m$species=="all",]$no.counted, ylim=c(0,80), type="l", lwd=2, col="white", ylab="Number of systems counted",xlab="Year", xaxt="n", yaxt="n")
# axis(1, at=seq(1950,2020,5), lab=seq(1950,2020,5))
# axis(2, at=seq(0,100,10), lab=seq(0,100,10), las=2)
# lines(unique(m$year), m[m$species=="Chum",]$no.counted, col="orange")
# lines(unique(m$year), m[m$species=="Coho",]$no.counted, col="green")
# lines(unique(m$year), m[m$species=="Chinook",]$no.counted, col="blue")
# lines(unique(m$year), m[m$species=="Sockeye",]$no.counted, col="red")
# lines(unique(m$year), m[m$species=="Pink",]$no.counted, col="magenta")
# lines(unique(m$year), m[m$species=="all",]$no.counted, lwd=2)
# 
# legend("topright", lty= rep(1,6), lwd=(c(2, rep(1,5))), col=c("black", "magenta","orange","green","blue","red"), c("All systems (80)", "Pink (68)", "Chum (63)","Coho (73)","Chinook (23)","Sockeye (34)"), bty="n", cex=0.7)
#showtext_auto(FALSE)
dev.off()

##################################################################################
### GENERATING REGIONALLY-SPECIFIC FIGURES 
##################################################################################
for (reg in unique(popdat$region)){
    
		p = popdat[popdat$region==reg,]$POP_ID
	
		# Set up dataframe for summary numbers #
		data <- nuseds[nuseds$POP_ID %in% p,]
		
		all.rivers <- unique(data$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
		n.rivers <- length(all.rivers)
		
		pink.rivers <- unique(data[data$SPECIES=="Pink",]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
		n.pink <- length(pink.rivers)
		
		podd.rivers <- unique(data[data$SPECIES=="Pink" & data$ANALYSIS_YR%%2==1,]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
		n.podd <- length(pink.rivers)
		
		peven.rivers <- unique(data[data$SPECIES=="Pink" & data$ANALYSIS_YR%%2==0,]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
		n.peven <- length(pink.rivers)
		
		chum.rivers <- unique(data[data$SPECIES=="Chum",]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
		n.chum <- length(chum.rivers)
		
		coho.rivers <- unique(data[data$SPECIES=="Coho",]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
		n.coho <- length(coho.rivers)
		
		chinook.rivers <- unique(data[data$SPECIES=="Chinook",]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
		n.chinook <- length(chinook.rivers)
		
		sock.rivers <- unique(data[data$SPECIES=="Sockeye",]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)
		n.sock <- length(sock.rivers)
		
		# Calculating annual monitoring coverage #
		ny <- length(1950:2021)
		
		m <- data.frame(year = rep(c(1950:2021), 6),
										species = rep(c("Pink","Chum", "Coho", "Chinook", "Sockeye","all"),each=ny),
										odd_even = c(rep(c("O","E"), length.out=ny), rep(NA, 5*ny)),
										total.systems = c(rep(n.podd, ny), rep(c(n.chum, n.coho, n.chinook, n.sock, n.rivers), each=ny)),
										no.counted = rep(NA, 6*ny),
										prop.counted = rep(NA, 6*ny))
		
		# Population data frame with monitoring coverage (number of systems/proportion of systems counted) #
		for (species in unique(m$species)) {
			if(species=="all") {
				for (i in 1950:2021){
					m[m$species==species & m$year==i,]$no.counted = length(unique(data[data$ANALYSIS_YR==i,]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE))
				}
			}
			else { 
				for (i in 1950:2021){
					m[m$species==species & m$year==i,]$no.counted = length(unique(data[data$SPECIES==species & data$ANALYSIS_YR==i,]$GEOGRAPHICAL_EXTNT_OF_ESTIMATE)) }
			}
		}
		
		
		m$prop.counted = m$no.counted/m$total.systems
		
		# Write to file #
		# dir.table <- "C:/Users/Emma Atkinson/Desktop/Research/SoS/Tables"
		# setwd(dir.table)
		# write.csv(m, "monitoring_coverage_table_Feb-2019.csv")
		
		# Plotting #
		pdf(paste("figures/River-level_monitoring_coverage_Jun-2023_",reg,".pdf",sep=""), 6,6)
		#tiff("monitoring_coverage_Sept-2019.tiff", width=4, height=4, units="in", pointsize=9, res=600)
		#windows(9,9)
		
		#showtext_auto()
		
		par(mfrow=c(1,1))#, family="Avenir")#, mgp=c(2.5, 0.8, 0), mar=c(3.5, 3.5, 1, 2))
		
		plot(unique(m$year), m[m$species=="all",]$no.counted, "l", las =1, bty="l", xlab = "Year", ylab="Number of systems monitored by DFO", lwd=1, ylim=c(0, 1.1*max(m$total.systems)), col=grey(0.8), main=paste(reg))
		# axis(1, at=seq(1950,2020,5), lab=seq(1950,2020,5))
		# axis(2, at=seq(0,70,10), lab=seq(0,1,.1), las=2)
		
		points(2021, tail(m$no.counted, 1), pch = 19, col="#2E8B57")
		text(2021, tail(m$no.counted, 1), paste("n=",tail(m$no.counted, 1),sep=""), font = 2, cex=0.8, xpd=NA, pos=4, col="#2E8B57")
		
		lines(unique(m$year), m[m$species=="Chum",]$no.counted, col="grey20", lty=2, lwd=2)
		lines(unique(m$year), m[m$species=="Coho",]$no.counted, col="black", lty=3, lwd=1)
		lines(unique(m$year), m[m$species=="Chinook",]$no.counted, col="grey40", lty=4, lwd=1)
		lines(unique(m$year), m[m$species=="Sockeye",]$no.counted, col="grey50", lty=5, lwd=1)
		lines(unique(m$year), m[m$species=="Pink",]$no.counted, col="black", lwd=2)
		lines(unique(m$year), m[m$species=="all",]$no.counted, lwd=4, col="#2E8B57")
		
		legend("topright", lty= c(1,1,2,3,4,5), lwd=c(4,2,1,1,1,1), col=c("#2E8B57", "black","grey20","black","grey40","grey50"), 
					 c("All systems (n=80)", "Pink", "Chum","Coho","Chinook","Sockeye"), bty="n", cex=.7)
		
		# plot(unique(m$year), m[m$species=="all",]$no.counted, ylim=c(0,80), type="l", lwd=2, col="white", ylab="Number of systems counted",xlab="Year", xaxt="n", yaxt="n")
		# axis(1, at=seq(1950,2020,5), lab=seq(1950,2020,5))
		# axis(2, at=seq(0,100,10), lab=seq(0,100,10), las=2)
		# lines(unique(m$year), m[m$species=="Chum",]$no.counted, col="orange")
		# lines(unique(m$year), m[m$species=="Coho",]$no.counted, col="green")
		# lines(unique(m$year), m[m$species=="Chinook",]$no.counted, col="blue")
		# lines(unique(m$year), m[m$species=="Sockeye",]$no.counted, col="red")
		# lines(unique(m$year), m[m$species=="Pink",]$no.counted, col="magenta")
		# lines(unique(m$year), m[m$species=="all",]$no.counted, lwd=2)
		# 
		# legend("topright", lty= rep(1,6), lwd=(c(2, rep(1,5))), col=c("black", "magenta","orange","green","blue","red"), c("All systems (80)", "Pink (68)", "Chum (63)","Coho (73)","Chinook (23)","Sockeye (34)"), bty="n", cex=0.7)
		#showtext_auto(FALSE)
		dev.off()
		
}			

##################################################################################
### GENERATING REGIONALLY-SPECIFIC FIGURES AT A CU-LEVEL
##################################################################################
for (reg in unique(popdat$region)){
	
		p = popdat[popdat$region==reg,]$POP_ID
		
		# Set up dataframe for summary numbers #
		data <- nuseds[nuseds$POP_ID %in% p,]
		
		ymax=2021
		ny=length(1950:2021)
		
		all.rivers <- unique(data$CU_NAME)
		n.rivers <- length(all.rivers)
		
		peven.cus <- unique(data[data$SPECIES=="Pink" & data$ANALYSIS_YR%%2==0,]$CU_NAME)
		n.peven <- length(peven.cus)
		
		podd.cus <- unique(data[data$SPECIES=="Pink" & data$ANALYSIS_YR%%2==1,]$CU_NAME)
		n.podd <- length(podd.cus)
		
		chum.cus <- unique(data[data$SPECIES=="Chum",]$CU_NAME)
		n.chum <- length(chum.cus)
		
		coho.cus <- unique(data[data$SPECIES=="Coho",]$CU_NAME)
		n.coho <- length(coho.cus)
		
		chinook.cus <- unique(data[data$SPECIES=="Chinook",]$CU_NAME)
		n.chinook <- length(chinook.cus)
		
		sock.cus <- unique(data[data$SPECIES=="Sockeye",]$CU_NAME)
		n.sock <- length(chum.cus)
		
		cu.frame = data.frame(year=c(1950:ymax),
													species = rep(c("Pink","Chum", "Coho", "Chinook", "Sockeye","all"),each=ny),
													odd_even = c(rep(c("E","O"), length.out=ny), rep(NA, 5*ny)),
													total.cus = c(rep(n.podd, ny), rep(c(n.chum, n.coho, n.chinook, n.sock, n.rivers), each=ny)),
													no.counted = rep(NA, 6*ny),
													prop.counted = rep(NA, 6*ny))
		
		for (k in 1:nrow(cu.frame)){
			
			spp = cu.frame$species[k]
			year = cu.frame$year[k]
			
			if (spp=="all"){ CUs = unique(data$CU_NAME) } else { CUs = unique(data[data$SPECIES==spp,]$CU_NAME) }
				
				temp = data[data$CU_NAME %in% CUs & data$ANALYSIS_YR==year,]
				monitored = c(rep(NA, length(CUs)))
			
			
			for (q in 1:length(unique(CUs))){
				monitored[q] = nrow((temp[temp$CU_NAME==CUs[q],]))
			}
			
			cu.frame[k,]$no.counted = length(monitored[monitored > 0])
			cu.frame[k,]$prop.counted = length(monitored[monitored > 0])/length(monitored)
			
		}
		
		# PLOTTING #
		pdf(paste("figures/CU-level_monitoring_coverage_Jun-2023_",reg,".pdf",sep=""), 6,6)
		par(mfrow=c(1,1))#, family="Avenir")#, mgp=c(2.5, 0.8, 0), mar=c(3.5, 3.5, 1, 2))
		
		plot(unique(cu.frame$year), cu.frame[cu.frame$species=="all",]$no.counted, "l", las =1, bty="l", xlab = "Year", ylab="Number of CUs with at least one pop'n counted", lwd=1, ylim=c(0, 1.1*max(cu.frame$total.cus)), col=grey(0.8), main=paste(reg))
		# axis(1, at=seq(1950,2020,5), lab=seq(1950,2020,5))
		# axis(2, at=seq(0,70,10), lab=seq(0,1,.1), las=2)
		
		points(2021, tail(cu.frame$no.counted, 1), pch = 19, col="#2E8B57")
		text(2021, tail(cu.frame$no.counted, 1), paste("n=",tail(cu.frame$no.counted, 1),sep=""), font = 2, cex=0.8, xpd=NA, pos=4, col="#2E8B57")
		
		lines(unique(cu.frame$year), cu.frame[cu.frame$species=="Chum",]$no.counted, col="grey20", lty=2, lwd=2)
		lines(unique(cu.frame$year), cu.frame[cu.frame$species=="Coho",]$no.counted, col="black", lty=3, lwd=1)
		lines(unique(cu.frame$year), cu.frame[cu.frame$species=="Chinook",]$no.counted, col="grey40", lty=4, lwd=1)
		lines(unique(cu.frame$year), cu.frame[cu.frame$species=="Sockeye",]$no.counted, col="grey50", lty=5, lwd=1)
		lines(unique(cu.frame$year), cu.frame[cu.frame$species=="Pink",]$no.counted, col="black", lwd=2)
		lines(unique(cu.frame$year), cu.frame[cu.frame$species=="all",]$no.counted, lwd=4, col="#2E8B57")
		
		legend("topright", lty= c(1,1,2,3,4,5), lwd=c(4,2,1,1,1,1), col=c("#2E8B57", "black","grey20","black","grey40","grey50"), 
					 c("All systems (n=80)", "Pink", "Chum","Coho","Chinook","Sockeye"), bty="n", cex=.7)
		
		dev.off()
}

################################################################################
### GENERATING WHOLE COAST CU FIGURE
################################################################################
data = nuseds

ymax=2021
ny=length(1950:2021)

all.rivers <- unique(data$CU_NAME)
n.rivers <- length(all.rivers)

peven.cus <- unique(data[data$SPECIES=="Pink" & data$ANALYSIS_YR%%2==0,]$CU_NAME)
n.peven <- length(peven.cus)

podd.cus <- unique(data[data$SPECIES=="Pink" & data$ANALYSIS_YR%%2==1,]$CU_NAME)
n.podd <- length(podd.cus)

chum.cus <- unique(data[data$SPECIES=="Chum",]$CU_NAME)
n.chum <- length(chum.cus)

coho.cus <- unique(data[data$SPECIES=="Coho",]$CU_NAME)
n.coho <- length(coho.cus)

chinook.cus <- unique(data[data$SPECIES=="Chinook",]$CU_NAME)
n.chinook <- length(chinook.cus)

sock.cus <- unique(data[data$SPECIES=="Sockeye",]$CU_NAME)
n.sock <- length(chum.cus)

cu.frame = data.frame(year=c(1950:ymax),
											species = rep(c("Pink","Chum", "Coho", "Chinook", "Sockeye","all"),each=ny),
											odd_even = c(rep(c("E","O"), length.out=ny), rep(NA, 5*ny)),
											total.cus = c(rep(n.podd, ny), rep(c(n.chum, n.coho, n.chinook, n.sock, n.rivers), each=ny)),
											no.counted = rep(NA, 6*ny),
											prop.counted = rep(NA, 6*ny))

for (k in 1:nrow(cu.frame)){
	
	spp = cu.frame$species[k]
	year = cu.frame$year[k]
	
	if (spp=="all"){ CUs = unique(data$CU_NAME) } else { CUs = unique(data[data$SPECIES==spp,]$CU_NAME) }
	
	temp = data[data$CU_NAME %in% CUs & data$ANALYSIS_YR==year,]
	monitored = c(rep(NA, length(CUs)))
	
	
	for (q in 1:length(unique(CUs))){
		monitored[q] = nrow((temp[temp$CU_NAME==CUs[q],]))
	}
	
	cu.frame[k,]$no.counted = length(monitored[monitored > 0])
	cu.frame[k,]$prop.counted = length(monitored[monitored > 0])/length(monitored)
	
}

# PLOTTING #
pdf(paste("figures/CU-level_monitoring_coverage_Jun-2023_whole_coast.pdf",sep=""), 6,6)
par(mfrow=c(1,1))#, family="Avenir")#, mgp=c(2.5, 0.8, 0), mar=c(3.5, 3.5, 1, 2))

plot(unique(cu.frame$year), cu.frame[cu.frame$species=="all",]$no.counted, "l", las =1, bty="l", xlab = "Year", ylab="Number of CUs with at least one pop'n counted", lwd=1, ylim=c(0, 1.1*max(cu.frame$total.cus)), col=grey(0.8), main="Coastwide")
# axis(1, at=seq(1950,2020,5), lab=seq(1950,2020,5))
# axis(2, at=seq(0,70,10), lab=seq(0,1,.1), las=2)

points(2021, tail(cu.frame$no.counted, 1), pch = 19, col="#2E8B57")
text(2021, tail(cu.frame$no.counted, 1), paste("n=",tail(cu.frame$no.counted, 1),sep=""), font = 2, cex=0.8, xpd=NA, pos=4, col="#2E8B57")

lines(unique(cu.frame$year), cu.frame[cu.frame$species=="Chum",]$no.counted, col="grey20", lty=2, lwd=2)
lines(unique(cu.frame$year), cu.frame[cu.frame$species=="Coho",]$no.counted, col="black", lty=3, lwd=1)
lines(unique(cu.frame$year), cu.frame[cu.frame$species=="Chinook",]$no.counted, col="grey40", lty=4, lwd=1)
lines(unique(cu.frame$year), cu.frame[cu.frame$species=="Sockeye",]$no.counted, col="grey50", lty=5, lwd=1)
lines(unique(cu.frame$year), cu.frame[cu.frame$species=="Pink",]$no.counted, col="black", lwd=2)
lines(unique(cu.frame$year), cu.frame[cu.frame$species=="all",]$no.counted, lwd=4, col="#2E8B57")

legend("topright", lty= c(1,1,2,3,4,5), lwd=c(4,2,1,1,1,1), col=c("#2E8B57", "black","grey20","black","grey40","grey50"), 
			 c("All systems (n=80)", "Pink", "Chum","Coho","Chinook","Sockeye"), bty="n", cex=.7)

dev.off()
										 