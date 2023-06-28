###############################################################################
###############################################################################
#
# Testing for broad-scale relationships between freshwater habitat pressure 
# indicators and Pacific salmon population trends
# Stephanie Peacock <speacock@psf.ca>
#
#
# Functions used in exploring and compiling nuSEDS CU data
# 
##############################################################################
###############################################################################


#------------------------------------------------------------------------------
# Function to plot timeseries for POP_ID = p
#------------------------------------------------------------------------------

# Points colour coded by Estimate Classification
estCol <- colorRampPalette(c("forest green", "gold", "red", "black"))(n = 10)

# Function for plot without legend
pp <- function(
	p, # POP_ID
	add = FALSE, # add points to current plot?
	pch = 19, 
	cex = 1,
	col = NA,
	xlim = NA, 
	ylim = NA, 
	dat = nusedsCU, # Can change to nuseds
	incl.legend = FALSE){
	
	if(is.na(xlim)) xlim <- range(dat$ANALYSIS_YR[dat$POP_ID == p], na.rm = TRUE)
	if(is.na(ylim)) ylim <- range(dat$MAX_ESTIMATE[dat$POP_ID == p], na.rm = TRUE)
	
	if(is.na(col)) col <- estCol[as.numeric(dat$ESTIMATE_CLASSIFICATION[dat$POP_ID == p])]
	
	if(add == TRUE){
		points(dat$ANALYSIS_YR[dat$POP_ID == p], dat$MAX_ESTIMATE[dat$POP_ID == p], col = col, pch = pch, cex = cex)
	} else {
		
		if(incl.legend == TRUE){
			par(mfrow = c(1,1), mar = c(4, 4, 2, 22), oma = rep(0, 4))
		} else {
			par(mfrow = c(1,1), mar = c(4, 4, 2, 2), oma = rep(0, 4))
		}
		
		plot(dat$ANALYSIS_YR[dat$POP_ID == p], dat$MAX_ESTIMATE[dat$POP_ID == p], xlab = "", ylab = "MAX_ESTIMATE", main = paste(dat$POPULATION[dat$POP_ID == p][1], " - POP_ID ", p, sep = ""), col = col, pch = pch, xlim = xlim, ylim = ylim, cex = cex)
		if(incl.legend == TRUE) legend(par('usr')[2] + (par('usr')[2] - par('usr')[1]) * 0.1, par('usr')[4], col = estCol, legend = levels(dat$ESTIMATE_CLASSIFICATION), pch = 19, xpd = NA)
		
	}
}

#------------------------------------------------------------------------------
# Function to count populations
#------------------------------------------------------------------------------

popCount <- function(){
	length(unique(nusedsCU$SQ_POP_ID))
}

# popFiltered <- data.frame(Step = 0, Descr = "Raw NuSEDS CU", nPop = popCount())
recordStep <- function(popFiltered = popFiltered, newDescr = NULL){
	
	if(length(newDescr) == 0){
		newDescr <- readline("Enter description for step:")
	}
	
	n <- popCount()
	
	popFiltered <- rbind(popFiltered, data.frame(
		Step = max(popFiltered$Step) + 1,
		Descr = newDescr,
		nRemoved = popFiltered$nPop[dim(popFiltered)[1]] - n,
		nPop = n))
	
	return(popFiltered)
}

#------------------------------------------------------------------------------
# Function to get mode (for stream order centering)
#------------------------------------------------------------------------------

getmode <- function(v) {
	uniqv <- unique(v)
	uniqv[which.max(tabulate(match(v, uniqv)))]
}

