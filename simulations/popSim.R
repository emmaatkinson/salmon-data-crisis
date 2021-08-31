###############################################################################
# This code is for a mini-simulation study to explore the combined effects of:
# 1) Fishing
# 2) Density dependence
# on our ability to detect shifts in productivity (e.g., due to changes in 
# habitat)

library(gsl) # for lambert_W0 function

###############################################################################
# Simulation functions
###############################################################################

b <- 1/100000

simPop <- function(a = 0.5, b = 1/100000, sigProc = sqrt(0.26), rho = 0.25, sigObs = sqrt(0.04), ny = 64, harvestRate = rep(0, 64), seed = NA){
	
	if(is.na(seed) == FALSE) set.seed(seed)
	
	if(length(a) == 1) a <- rep(a, ny)
	if(length(b) == 1) b <- rep(b, ny)
	
	spawners <- rep(NA, ny)
	spawners[1:4] <- runif(4, 20000, 80000)
	recruits <- rep(NA, ny)
	
	nu <- rnorm(ny, 0, sd = sigProc)
	epsilon <- numeric(ny)
	
	for(i in 5:ny){
		epsilon[i] <- rho * epsilon[i-1] + nu[i]
		recruits[i] <- spawners[i - 4] * exp(a[i] - b[i] * spawners[i - 4] + epsilon[i])
		spawners[i] <- (1 - harvestRate[i]) * recruits[i]
	}	
	
	spawnersObs <- spawners * exp(rnorm(ny, mean = 0, sd = sigObs))
	
	return(data.frame(
		year = 1:ny, 
		# period = c(rep("init", 12), rep("eval", 52), rep("post", 12)), 
		spawnersObs,
		spawners, 
		recruits, 
		epsilon, 
		harvestRate, 
		a, 
		b))
}

calcSmsy <- function(a, b) {
	
	Smsy = (1 - lambert_W0(exp(1 - a))) / b
	
	return(as.numeric(Smsy))
}

smoothGen <- function(x, gen = 4){
	xSmoothed <- rep(NA, length(x))
	for(i in 1:length(x)){
		ind <- c(max(1, i - gen + 1):i)
		xSmoothed[i] <- prod(x[ind]) ^ (1/length(ind))
	}
	return(xSmoothed)
}

calcDecl <- function(x, gen = 4, assessmentYear = 64, metric = "PercentDecline"){
	
	xSmoothed <- smoothGen(x, gen)
	
	if(metric == "PercentDecline"){
		xLogged <- log(xSmoothed)
	
		mod <- lm(xLogged[13:assessmentYear] ~ c(13:assessmentYear))
		
		xPredict <- exp(predict(mod))
		percDecl <- (tail(xPredict, 1) - xPredict[1]) / xPredict[1]
		
		return(round(percDecl, 4))
		
	} else if(metric == "AnnualChange"){
	
		mod <- lm(xSmoothed[13:assessmentYear] ~ c(13:assessmentYear))
		return(round(coefficients(mod)[2], 4))
}
}

#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------
plotSim <- function(sim1, change = "a"){
	par(mfrow= c(3,1), mar = c(0,5,0,5), oma = c(4, 0, 2, 0))

	#------------------------------------------------------------------------------
	# Harvest rate and productivity
	plot(1:ny, sim1$harvestRate, "o", col = grey(0.6), las = 1, ylab = "", xaxt="n", lwd = 1.5)
	axis(side = 1, labels = FALSE)
	mtext(side = 2, line=3.5, "Harvest rate", col = grey(0.6))
	
	abline(v = c(12.5, 64.5), lty = 2, lwd = 2)
	
	par(new = TRUE)
	
	if(change == "a") {
		plot(1:ny, sim1$a, "o", pch = 19, cex = 0.5, col = 4, bty = "n", xlab = "", ylab = "", xaxt="n", yaxt="n")
		axis(side = 4, las = 1, col = 4)
		abline(h = 0, lty = 3, col = 4)
		mtext(side = 4, "Productivity (a)", col= 4, line = 3)
	}
	if(change == "b"){
		plot(1:ny, sim1$b, "o", pch = 19, cex = 0.5, col = 3, bty = "n", xlab = "", ylab = "", xaxt="n", yaxt="n")
		axis(side = 4, las = 1, col = 3)
		abline(h = 0, lty = 3, col = 3)
		mtext(side = 4, "Density dependence (b)", col= 3, line = 3)
	}
	
	#------------------------------------------------------------------------------
	# Spawners and recruits
	
	plot(1:ny, sim1$spawnersObs, "o", cex = 0.8, col = grey(0.6), xlab = "", ylab = "", las = 1, ylim = range(sim1$spawners, sim1$recruits, na.rm = TRUE),  pch = 21, bg ="white", xaxt="n")
	abline(v = c(12.5, 64.5), lty = 2, lwd = 2)
	lines(1:ny, sim1$spawners,"o", cex=0.8, pch = 19)
	axis(side = 1, labels = FALSE)
	lines(1:(ny-4), sim1$recruits[5:ny], "o", cex = 0.8, col = 2, pch = 19)
	legend("topleft", pch = c(19, 21, 19), pt.bg = c(NA, "white", NA), pt.cex = 0.8, lty = 1, col = c(1,grey(0.6), 2), c("Spawners", "Obs. spawners", "Recruits"))
	lines(1:ny, sim1$a/b, lty = 3)
	mtext(side = 2, "Abundance", line = 3.5)
	
	#------------------------------------------------------------------------------
	# Smoothed abundance
	ySmoothed <- smoothGen(sim1$spawnersObs)
	yLogged <- log(ySmoothed)
	mod <- lm(yLogged[13:64] ~ c(13:64))
	
	plot(1:ny, sim1$spawnersObs, "o", pch = 21, col= grey(0.6), cex = 0.8, las = 1, ylab = "")
	abline(v = c(12.5, 64.5), lty = 2, lwd = 2)
	lines(c(1:ny), ySmoothed, lwd = 1.5)
	u <- par('usr')
	points(c(13, 64), exp(c(predict(mod)[1], tail(predict(mod), 1))), cex = 2, pch = 19)
	segments(x0 = rep(u[1], 2), x1 = c(13, 64), y0 = exp(c(predict(mod)[1], tail(predict(mod), 1))), y1= exp(c(predict(mod)[1], tail(predict(mod), 1))))
	
	text(1, exp(predict(mod)[1]), round(exp(predict(mod)[1])), pos = c(1,3)[as.numeric(predict(mod)[1]>tail(predict(mod),1))+1])
	text(1, exp(tail(predict(mod), 1)), round(exp(tail(predict(mod),1))), pos = c(3,1)[as.numeric(predict(mod)[1]>tail(predict(mod),1))+1])
	
	arrows(x0=1, x1=1, y0 = exp(predict(mod)[1]), y1 = exp(tail(predict(mod), 1)), length = 0.08, code = 3)
	text(1, exp(mean(predict(mod))), round((exp(tail(predict(mod), 1))-exp(predict(mod)[1]))/exp(predict(mod)[1]), 3), pos = 4)
	lines(c(13:64), exp(predict(mod)), lwd = 2, lty = 2)
	
	mtext(side = 2, "Smoothed abundance", line = 3.5)
	
	axis(side = 4, at = exp(seq(round(min(log(sim1$spawnersObs))), round(max(log(sim1$spawnersObs))), 1)), labels = seq(round(min(log(sim1$spawnersObs))), round(max(log(sim1$spawnersObs))), 1), las = 1)
	mtext(side = 4, "Log(smoothed abundance)", line = 3.5)
	
	mtext(side=1, line = 2.5, "Brood year")
}

###############################################################################
# Metric explanation
###############################################################################
ny <- 64

sim1 <- simPop(a = rep(0.3, ny), harvestRate = rep(0.2, ny), seed = 4589)

#------------------------------------------------------------------------------
# Metric #1 : Metrc 4 from d'Eon Eggerson

ySmoothed <- smoothGen(sim1$spawnersObs)
yLogged <- log(ySmoothed)
mod <- lm(yLogged[13:64] ~ c(13:64))


par(mfrow= c(3,1), mar = c(0,5,0,5), oma = c(4, 0, 2, 0))

plot(1:ny, sim1$spawnersObs, "o", pch = 21, col= grey(0.6), cex = 0.8, las = 1, ylab = "", xaxt="n")
abline(v = c(12.5, 64.5), lty = 2, lwd = 2)
mtext(side = 2, "Smoothed abundance", line = 3.5)
lines(c(1:ny), ySmoothed, lwd = 1.5)

plot(1:ny, yLogged, "l", col = 3, lwd = 1.5, las = 1, ylab = "", xaxt="n", yaxt="n")
abline(v = c(12.5, 64.5), lty = 2, lwd = 2)
u <- par('usr')
axis(side = 4, col = 3, las = 1)
mtext(side = 4, "Log(smoothed abundance)", line = 3.5)
lines(c(13:64), predict(mod), col = 3, lwd = 2, lty = 2)
points(c(13, 64), c(predict(mod)[1], tail(predict(mod), 1)), cex = 2, pch = 21, bg = 3)

plot(1:ny, sim1$spawnersObs, "n", las = 1, ylab = "", xaxt="n")
abline(v = c(12.5, 64.5), lty = 2, lwd = 2)
lines(c(1:ny), ySmoothed, lwd = 1.5)
mtext(side = 2, "Smoothed abundance", line = 3.5)
points(c(13, 64), exp(c(predict(mod)[1], tail(predict(mod), 1))), cex = 2, pch = 19)
segments(x0 = rep(u[1], 2), x1 = c(13, 64), y0 = exp(c(predict(mod)[1], tail(predict(mod), 1))), y1= exp(c(predict(mod)[1], tail(predict(mod), 1))))
text(1, exp(predict(mod)[1]), round(exp(predict(mod)[1])), pos = 3)
text(1, exp(tail(predict(mod), 1)), round(exp(tail(predict(mod),1))), pos = 1)
arrows(x0=1, x1=1, y0 = exp(predict(mod)[1]), y1 = exp(tail(predict(mod), 1)), length = 0.08, code = 3)
text(1, exp(mean(predict(mod))), round((exp(tail(predict(mod), 1)) - exp(predict(mod)[1]))/exp(predict(mod)[1]), 3), pos = 4)
# lines(c(13:64), exp(predict(mod)), lwd = 2, lty = 2)

#------------------------------------------------------------------------------
# Metric #2: Annual rate of change for standardized timeseries
# yStandard <- ySmoothed - mean(ySmoothed, na.rm = TRUE)
mod2 <- lm(ySmoothed[13:64] ~ c(13:64))

par(mfrow= c(2,1), mar = c(0,5,0,5), oma = c(4, 0, 2, 0))

plot(1:ny, sim1$spawnersObs, "o", pch = 21, col= grey(0.6), cex = 0.8, las = 1, ylab = "", xaxt="n")
abline(v = c(12.5, 64.5), lty = 2, lwd = 2)
mtext(side = 2, "Smoothed abundance", line = 3.5)
lines(c(1:ny), ySmoothed, lwd = 1.5)

# Standardized
plot(1:ny, ySmoothed, "l", lwd = 1.5,  yaxt = "n", ylim = range(ySmoothed[13:64]), ylab = "")
abline(v = c(12.5, 64.5), lty = 2, lwd = 2)
u <- par('usr')
axis(side = 4, col = 4, las = 1)
lines(13:64, predict(mod2), lty = 2, lwd = 2, col = 4)
text(45, 15000, substitute(paste(Delta, "y/", Delta, "x = ", D), list(D = round(coefficients(mod2)[2], 2))), cex =1.5, col = 4)


#------------------------------------------------------------------------------
# Compare over 1000 simulations
popMetric <- matrix(NA, nrow = 1000, ncol = 2)

for(i in 1:1000){
	sim.i <- simPop(a = rep(rnorm(1, 0.8, 0.5), ny), harvestRate = rep(0.2, ny))
	
	ySmoothed <- smoothGen(sim.i$spawnersObs)
	yLogged <- log(ySmoothed)
	mod <- lm(yLogged[13:64] ~ c(13:64))
	mod2 <- lm(ySmoothed[13:64] ~ c(13:64))
	
	popMetric[i, 1] <- (exp(tail(predict(mod), 1)) - exp(predict(mod)[1]))/exp(predict(mod)[1])
	popMetric[i, 2] <- 	coefficients(mod2)[2]																																				 
}

par(mfrow= c(1,1), mar = c(4,4,2,1), oma = rep(0, 4))
plot(popMetric[, 1], popMetric[, 2], xlab = "Percent decline", ylab = "Annual rate of change")
abline(h = 0)
abline(v = 0)

cor.test(popMetric[, 1], popMetric[, 2])
# Couple of examples
sim.j <- simPop(a = rep(0.3, ny), harvestRate = rep(0.2, ny))
ySmoothed <- smoothGen(sim.j$spawnersObs)
yLogged <- log(ySmoothed)
mod <- lm(yLogged[13:64] ~ c(13:64))
mod2 <- lm(ySmoothed[13:64] ~ c(13:64))

par(mfrow= c(1,1), mar = c(0,5,0,5), oma = c(4, 0, 2, 0))

plot(1:ny, sim.j$spawnersObs, "o", pch = 21, col= grey(0.6), cex = 0.8, las = 1, ylab = "", ylim = range(ySmoothed[13:64], exp(c(predict(mod)[1], tail(predict(mod), 1)))))
abline(v = c(12.5, 64.5), lty = 2, lwd = 2)
mtext(side = 2, "Smoothed abundance", line = 3.5)
lines(c(1:ny), ySmoothed, lwd = 1.5)

lines(13:64, predict(mod2), lty = 2, lwd = 2, col = 4)
text(45, coefficients(mod2)[1] + coefficients(mod2)[2]*45, substitute(paste(Delta, "y/", Delta, "x = ", D), list(D = round(coefficients(mod2)[2], 2))), cex =1.5, col = 4, pos = 4, xpd = NA)

points(c(13, 64), exp(c(predict(mod)[1], tail(predict(mod), 1))), cex = 2, pch = 21, bg = 3)

par(new = TRUE)
plot(1:ny, yLogged, "l", col = 3, lwd = 1.5, ylim = range(predict(mod)), yaxt="n", xaxt = "n", bty = "n", ylab = "")
lines(c(13:64), predict(mod), col = 3, lwd = 2, lty = 2)
u <- par('usr')
segments(x0 = 13, x1 = u[2]+2, y0 = predict(mod)[1], y1 = predict(mod)[1], lty = 3, xpd =NA)
segments(x0 = 64, x1 = u[2]+2, y0 = tail(predict(mod),1), y1 = tail(predict(mod),1), lty = 3, xpd = NA)
text(u[2]+3, u[3] + (u[4]-u[3])/2, round((exp(tail(predict(mod), 1)) - exp(predict(mod)[1]))/exp(predict(mod)[1]), 3), xpd = NA, col = 3, cex = 1.5)
###############################################################################
# Base simulation to illustrate metrics
#    - constant productivity
#    - constant fishing pressure
###############################################################################

#------------------------------------------------------------------------------
# Example plot
#------------------------------------------------------------------------------

sim1 <- simPop(a = rep(1, ny), harvestRate = rep(0.2, ny), seed = 987)
plotSim(sim1)

#------------------------------------------------------------------------------
# Over 1000 simulations
#------------------------------------------------------------------------------
change1 <- numeric(1000)
for(i in 1:1000){
	sim.i <- simPop(a = rep(1, ny), harvestRate = rep(0.2, ny))
	change1[i] <- calcDecl(x = sim.i$spawnersObs)
}

quartz(width = 3.5, height = 4, pointsize = 10)
par(mfrow = c(1,1), mar = c(4,4,2,1))
hist(change1, breaks = seq(-2, 5, 0.1), xlim = c(-1.2, 2.2), las = 1, yaxs = "i", main = "", xlab = "Percent change")
abline(v = 0, col = 2, lwd = 1.5)
# abline(v = median(change1))

###############################################################################
#    - linear decline in productivity
#    - constant fishing pressure
###############################################################################

#------------------------------------------------------------------------------
# Example plot
#------------------------------------------------------------------------------

sim2 <- simPop(a = c(rep(1, 12), seq(1, -0.2, length.out = ny - 12)), harvestRate = rep(0.2, ny), seed = 987)
plotSim(sim2)

#------------------------------------------------------------------------------
# Over 1000 simulations
#------------------------------------------------------------------------------
change2 <- numeric(1000)
for(i in 1:1000){
	sim.i <- simPop(a = c(rep(1, 12), seq(1, -0.2, length.out = ny - 12)), harvestRate = rep(0.2, ny))
	change2[i] <- calcDecl(x = sim.i$spawnersObs)
}

quartz(width = 3.5, height = 4, pointsize = 10)
par(mfrow = c(1,1), mar = c(4,4,2,1))
hist(change2, breaks = seq(-2, 5, 0.1), xlim = c(-1.2, 2.2), las = 1, yaxs = "i", main = "", xlab = "Percent change")
abline(v = 0, col = 2, lwd = 1.5)
# abline(v = median(change1))

#------------------------------------------------------------------------------
# Offset for change in prod
#------------------------------------------------------------------------------
start.year <- c(13:52)
change2a <- matrix(nrow = 1000, ncol = length(start.year))
for(j in 1:length(start.year)){
	for(i in 1:1000){
		sim.i <- simPop(a = c(rep(1, start.year[j] - 1), seq(1, -0.2, length.out = (ny - start.year[j] + 1))), 
										harvestRate = rep(0.2, ny))
		change2a[i, j] <- calcDecl(x = sim.i$spawnersObs)
	}
}


FP2a <- numeric(length(start.year))
for(j in 1:length(start.year)){
	FP2a[j] <- round(length(which(change2a[, j] > 0))/1000, 3)
}

###############################################################################
#    - linear decline in productivity
#    - linear decline in fishing pressure
###############################################################################

#------------------------------------------------------------------------------
# Example plot
#------------------------------------------------------------------------------

sim3 <- simPop(
	a = c(rep(1, 12), seq(1, -0.2, length.out = ny - 12)), 
	harvestRate = c(rep(0.6, 12), seq(0.6, 0.1, length.out = ny - 12)), 
	seed = 987)
plotSim(sim3)

#------------------------------------------------------------------------------
# Over 1000 simulations
#------------------------------------------------------------------------------
change3 <- numeric(1000)
for(i in 1:1000){
	sim.i <- simPop(a = c(rep(1, 12), seq(1, -0.2, length.out = ny - 12)), 
									harvestRate = c(rep(0.6, 12), seq(0.6, 0.1, length.out = ny - 12)))
	change3[i] <- calcDecl(x = sim.i$spawnersObs)
}

quartz(width = 3.5, height = 4, pointsize = 10)
par(mfrow = c(1,1), mar = c(4,4,2,1))
hist(change3, breaks = seq(-2, 5, 0.1), xlim = c(-1.2, 2.2), las = 1, yaxs = "i", main = "", xlab = "Percent change")
abline(v = 0, col = 2, lwd = 1.5)
# abline(v = median(change1))
text(1, 150, paste("FP = ", round(length(which(change3 > 0))/length(change3), 3)))

#------------------------------------------------------------------------------
# Offset for change in prod
#------------------------------------------------------------------------------
start.year <- c(13:52)
change3a <- matrix(nrow = 1000, ncol = length(start.year))
for(j in 1:length(start.year)){
	for(i in 1:1000){
		sim.i <- simPop(a = c(rep(1, start.year[j] - 1), seq(1, -0.2, length.out = (ny - start.year[j] + 1))), 
										harvestRate = c(rep(0.6, 12), seq(0.6, 0.1, length.out = ny - 12)))
		change3a[i, j] <- calcDecl(x = sim.i$spawnersObs)
	}
}


FP3a <- numeric(length(start.year))
for(j in 1:length(start.year)){
	FP3a[j] <- round(length(which(change3a[, j] > 0))/1000, 3)
}

par(mfrow = c(1,1), mar = c(4,4,2,1), oma = rep(0, 4))
plot(start.year, FP3a, "o", las = 1, ylab = "Proportion of trials with positive trend", xlab = "Year decline starts", ylim = c(0, 1), yaxs = "i", pch = 19)
points(start.year, FP2a, "o", pch = 19, xpd = NA, col = 2)
legend("topleft", pch = 19, lwd = 1, col = c(1,2), c("Decline in harvest", "Constant harvest"))


#---------------------
# Plot single simulation to illustrate
j <- 25
sim3a <- simPop(a = c(rep(1, start.year[j] - 1), seq(1, -0.2, length.out = (ny - start.year[j] + 1))), 
									 harvestRate = c(rep(0.6, 12), seq(0.6, 0.1, length.out = ny - 12)), seed = 987)
plotSim(sim3a)

par(mfrow = c(1,1), mar = c(4,4,2,1))
hist(change3a[, j], breaks = seq(-2, max(change3a[, j])+0.1, 0.1), xlim = c(-1.2, 2.2), las = 1, yaxs = "i", main = "", xlab = "Percent change")
abline(v = 0, col = 2, lwd = 1.5)
# abline(v = median(change1))
text(1.5, 60, paste("FP = ", FP3a[j]))

#------------------------------------------------------------------------------
# Offset for change in harvest
#------------------------------------------------------------------------------
start.year <- c(13:52)
change3b <- matrix(nrow = 1000, ncol = length(start.year))
for(j in 1:length(start.year)){
	for(i in 1:1000){
		sim.i <- simPop(a = c(rep(1, 12), seq(1, -0.2, length.out = ny - 12)), 
										harvestRate = c(rep(0.6, start.year[j] - 1), seq(0.6, 0.1, length.out = (ny - start.year[j] + 1))))
		change3b[i, j] <- calcDecl(x = sim.i$spawnersObs)
	}
}


FP3b <- numeric(length(start.year))
for(j in 1:length(start.year)){
	FP3b[j] <- round(length(which(change3b[, j] > 0))/1000, 3)
}

par(mfrow = c(1,1), mar = c(4,4,2,1), oma = rep(0, 4))
plot(start.year, FP3b, "o", las = 1, ylab = "Proportion of trials with positive trend", xlab = "Year decline starts", ylim = c(0, 1), yaxs = "i", pch = 19)


#---------------------
# Plot single simulation to illustrate
j <- 25
sim3b <- simPop(a = c(rep(1, 12), seq(1, -0.2, length.out = ny - 12)), 
								harvestRate = c(rep(0.6, start.year[j] - 1), seq(0.6, 0.1, length.out = (ny - start.year[j] + 1))), seed = 987)
plotSim(sim3b)

par(mfrow = c(1,1), mar = c(4,4,2,1))
hist(change3a[, j], breaks = seq(-2, max(change3a[, j])+0.1, 0.1), xlim = c(-1.2, 2.2), las = 1, yaxs = "i", main = "", xlab = "Percent change")
abline(v = 0, col = 2, lwd = 1.5)
# abline(v = median(change1))
text(1.5, 60, paste("FP = ", FP3a[j]))


###############################################################################
#    - step down in productivity
#    - linear decline in fishing pressure
###############################################################################

#------------------------------------------------------------------------------
# Example plot
#------------------------------------------------------------------------------

sim4 <- simPop(
	a = c(rep(1, 32), rep(-0.2, 32)), 
	harvestRate = c(rep(0.6, 12), seq(0.6, 0.1, length.out = ny - 12)), 
	seed = 987)
plotSim(sim4)

#------------------------------------------------------------------------------
# Over 1000 simulations
#------------------------------------------------------------------------------
change4 <- numeric(1000)
for(i in 1:1000){
	sim.i <- simPop(a = c(rep(1, 32), rep(-0.2, 32)), 
									harvestRate = c(rep(0.6, 12), seq(0.6, 0.1, length.out = ny - 12)))
	change4[i] <- calcDecl(x = sim.i$spawnersObs)
}

quartz(width = 3.5, height = 4, pointsize = 10)
par(mfrow = c(1,1), mar = c(4,4,2,1))
hist(change4, breaks = seq(-2, 5, 0.1), xlim = c(-1.2, 2.2), las = 1, yaxs = "i", main = "", xlab = "Percent change")
abline(v = 0, col = 2, lwd = 1.5)
# abline(v = median(change1))
text(1, 150, paste("FP = ", round(length(which(change4 > 0))/1000, 3)))

#------------------------------------------------------------------------------
# Effect of step size
#------------------------------------------------------------------------------
end.a <- seq(-0.2, 2, 0.1)
change4a <- matrix(nrow = 1000, ncol = length(end.a))
for(j in 1:length(end.a)){
	for(i in 1:1000){
		sim.i <- simPop(a = c(rep(1, 32), rep(end.a[j], 32)), 
									harvestRate = c(rep(0.6, 12), seq(0.6, 0.1, length.out = ny - 12)))
		change4a[i, j] <- calcDecl(x = sim.i$spawnersObs)
	}
}

FP4a <- numeric(length(end.a))
for(j in 1:length(end.a)){
	FP4a[j] <- round(length(which(change4a[, j] > 0))/1000, 3)
}
par(mfrow = c(1,1), mar = c(4,4,2,1), oma = rep(0, 4))
plot(end.a, FP4a, "o", las = 1, ylab = "Proportion of trials with positive trend", xlab = "Final productivity (a)", ylim = c(0, 1), yaxs = "i", pch = 19, xpd = NA)
u <- par('usr')
abline(v = 1, col = 2)
abline(h = 0.5, lty=2, col = 2)
text(1, 0.3, "Initial productivity\n a = 1.0", col = 2, pos = 4)

###############################################################################
#    Bring all scenarios together and calculate trend detection
###############################################################################

#------------------------------------------------------------------------------
# Productivity & harvest scenarios
#------------------------------------------------------------------------------
prod.all <- matrix(nrow = ny, ncol = 8)
prod.all[, 1] <- rep(1.0, ny)

prod.all[, 2] <- c(rep(1, 12), seq(1, 0, length.out = ny - 12))
prod.all[, 3] <- c(rep(1, 12*2), seq(1, 0, length.out = ny - 12*2))
prod.all[, 4] <- c(rep(1, 12*3), seq(1, 0, length.out = ny - 12*3))

prod.all[, 5] <- c(rep(1, 12), rep(0, ny - 12))
prod.all[, 6] <- c(rep(1, 12*2), rep(0, ny - 12*2))
prod.all[, 7] <- c(rep(1, 12*3), rep(0, ny - 12*3))
prod.all[, 8] <- c(rep(1, 12*4), rep(0, ny - 12*4))

harv.all <- matrix(nrow = ny, ncol = 3)
harv.all[, 1] <- rep(0.2, ny)
harv.all[, 2] <- c(rep(0.6, 12), seq(0.6, 0.1, length.out = ny - 12))
harv.all[, 3] <- c(rep(0.6, 12*2), rep(0.1, ny - 12*2))

posChange_prod <- list(PercentDecline = matrix(NA, nrow = dim(b.all)[2], ncol = dim(harv.all)[2]), AnnualChange = matrix(NA, nrow = dim(b.all)[2], ncol = dim(harv.all)[2]))

for(i in 1:dim(prod.all)[2]){
	for(j in 1:dim(harv.all)[2]){
		
		change.ij <- matrix(NA, nrow = 1000, ncol = 2)
		for(n in 1:1000){
			sim.ijn <- simPop(a = prod.all[, i], harvestRate = harv.all[, j])
			change.ij[n, 1] <- calcDecl(x = sim.ijn$spawnersObs, metric = "PercentDecline")
			change.ij[n, 2] <- calcDecl(x = sim.ijn$spawnersObs, metric = "AnnualChange")
		}
		
		posChange_prod$PercentDecline[i, j] <- length(which(change.ij[, 1] >= 0))/length(change.ij[, 1])
		posChange_prod$AnnualChange[i, j] <- length(which(change.ij[, 2] >= 0))/length(change.ij[, 2])
	}
}

#------------------------------------------------------------------------------
# Density dependence & harvest scenarios
#------------------------------------------------------------------------------

b.all <- matrix(nrow = ny, ncol = 8)
b.all[, 1] <- rep(10^-5, ny)

b.all[, 2] <- c(rep(10^-5, 12), 1/seq(100000, 10000, length.out = ny - 12))
b.all[, 3] <- c(rep(10^-5, 12*2), 1/seq(100000, 10000, length.out = ny - 12*2))
b.all[, 4] <- c(rep(10^-5, 12*3), 1/seq(100000, 10000, length.out = ny - 12*3))

b.all[, 5] <- c(rep(10^-5, 12), rep(10^-4, ny - 12))
b.all[, 6] <- c(rep(10^-5, 12*2), rep(10^-4, ny - 12*2))
b.all[, 7] <- c(rep(10^-5, 12*3), rep(10^-4, ny - 12*3))
b.all[, 8] <- c(rep(10^-5, 12*4), rep(10^-4, ny - 12*4))

posChange_densDep <- list(PercentDecline = matrix(NA, nrow = dim(b.all)[2], ncol = dim(harv.all)[2]), AnnualChange = matrix(NA, nrow = dim(b.all)[2], ncol = dim(harv.all)[2]))

for(i in 1:dim(b.all)[2]){
	for(j in 1:dim(harv.all)[2]){
		
		change.ij <- matrix(NA, nrow = 1000, ncol = 2)
		for(n in 1:1000){
			sim.ijn <- simPop(a = 1, b = b.all[, i], harvestRate = harv.all[, j])
			change.ij[n, 1] <- calcDecl(x = sim.ijn$spawnersObs, metric = "PercentDecline")
			change.ij[n, 2] <- calcDecl(x = sim.ijn$spawnersObs, metric = "AnnualChange")
		}
		
		posChange_densDep$PercentDecline[i, j] <- length(which(change.ij[, 1] >= 0))/length(change.ij[, 1])
		posChange_densDep$AnnualChange[i, j] <- length(which(change.ij[, 2] >= 0))/length(change.ij[, 2])
	}
}

#------------------------------------------------------------------------------
# What is going on in scenario 5, with a sudden step down in productivity?

simb5 <- simPop(
	a = 1,
	b = b.all[, 5], 
	harvestRate = rep(0.2, ny))#, 
	# seed = 987)
plotSim(simb5, change = "b")
