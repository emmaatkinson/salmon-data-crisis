x <- seq(-8, 8, 0.01)

par(mfrow = c(1,1), mar = c(4,4,2,6), oma = rep(0, 4))

plot(x, dnorm(x, mean = 0, sd = 1), "l", lwd = 2, xlim = c(-3, 4.5), ylim = c(0, 0.5), bty = "l", las =1, ylab = "Marginal posterior density", xlab = "Parameter value")
abline(v = 0, lwd = 1, lty = 3)

points(0, 0.5, lwd = 2, xpd = NA)
for(j in 1:3){
	segments(
		x0 = qnorm((1 - thresh[j])/2, 0, 1), 
		x1 = qnorm(1-(1 - thresh[j])/2, 0, 1),
		y0 = 0.5, y1 = 0.5, xpd = NA, lwd = 1 + c(j - 1)*2, lend = 1)
}

# Thresholds
thresh <- c(0.95, 0.80, 0.65)

for(i in 1:3){
	y <- qnorm(1-(1 - thresh[i])/2, 0, 1)
	# lines(x, dnorm(x, mean = qnorm((1 - thresh[i])/2, 0, 1), sd = 1), col = c(3,4,2)[i], lwd = 2)
	lines(x, dnorm(x, mean = y, sd = 1), col = c(3,4,2)[i], lwd = 2)
	
	points(y, 0.42+c(0, 0.025, 0.05)[i], col = c(3,4,2)[i], lwd = 2)
	
	for(j in 1:3){
		segments(
			x0 = qnorm((1 - thresh[j])/2, y, 1), 
			x1 = qnorm(1-(1 - thresh[j])/2, y, 1),
			y0 = 0.42+c(0, 0.025, 0.05)[i], 
			y1 = 0.42+c(0, 0.025, 0.05)[i], 
			col = c(3,4,2)[i], xpd = NA, lwd = 1 + c(j - 1)*2, lend = 1)
	}
	}
	
}

legend(4.5, 0.33, fill = c(3,4,2,1), title = "Evidence", legend = c("Strong", "Moderate", "Weak", "No"), border = NA, xpd = NA)
legend(4.5, 0.5, lwd = c(1,3,5), title = "Credible interval", legend = c("95%", "80%", "65%"), xpd = NA)


#----------------
# Same effect size, different sd
#----------------

x <- seq(-8, 8, 0.01)

par(mfrow = c(1,1), mar = c(4,4,2,6), oma = rep(0, 4))

plot(x, dnorm(x, mean = 1.5, sd = 2), "n", xlim = c(-3, 5), ylim = c(0, 0.8), bty = "l", las =1, ylab = "Marginal posterior density", xlab = "Parameter value")
abline(v = 0, lwd = 1, lty = 3)
abline(v = 1.5, col = grey(0.8), lwd = 3)

for(i in 1:4){
	
	lines(x, dnorm(x, mean = 1.5, sd = c(2,1.5, 1, 0.5)[i]), col = c(1,3,4,2)[i], lwd = 2)

	for(j in 1:3){
		segments(
			x0 = qnorm((1 - thresh[j])/2, 1.5, c(2, 1.5, 1, 0.5)[i]), 
			x1 = qnorm(1-(1 - thresh[j])/2, 1.5, c(2, 1.5, 1, 0.5)[i]),
			y0 = 0.8+c(0, 0.025, 0.05, 0.075)[i], 
			y1 = 0.8+c(0, 0.025, 0.05, 0.075)[i], 
			col = c(1,3,4,2)[i], xpd = NA, lwd = 1 + c(j - 1)*2, lend = 1)
	}
}
legend(5, 0.4, fill = c(2,4,3,1), title = "Evidence", legend = c("Strong", "Moderate", "Weak", "No"), border = NA, xpd = NA)
legend(5, 0.7, lwd = c(1,3,5), title = "Credible interval", legend = c("95%", "80%", "65%"), xpd = NA)

#----------------
# Compare
#----------------
plot(x, dnorm(x, mean = 1.5, sd = 2), "n", xlim = c(-3, 7), ylim = c(0, 0.8), bty = "l", las =1, ylab = "Marginal posterior density", xlab = "Parameter value")
abline(v = 0, lwd = 1, lty = 3)

lines(x, dnorm(x, mean = 1.5, sd = 0.7655), col = 2, lwd = 2)
lines(x, dnorm(x, mean = qnorm(1-(1 - 0.95)/2, 0, 2), sd = 2), col = 2, lty = 2)
	
	
	for(j in 1:3){
		segments(
			x0 = qnorm((1 - thresh[j])/2, 1.5, c(2, 1.5, 1, 0.5)[i]), 
			x1 = qnorm(1-(1 - thresh[j])/2, 1.5, c(2, 1.5, 1, 0.5)[i]),
			y0 = 0.8+c(0, 0.025, 0.05, 0.075)[i], 
			y1 = 0.8+c(0, 0.025, 0.05, 0.075)[i], 
			col = c(1,3,4,2)[i], xpd = NA, lwd = 1 + c(j - 1)*2, lend = 1)
	}
}