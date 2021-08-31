# ###############################################################################
# # What does the logistic look like
# ###############################################################################
# plot(c(0:100), -0.1 - 0.2 / (1 + exp(-0.8 * (c(0:100) - 20))), "l", lwd = 3, ylim = c(-0.3, 0), xaxt = "n", yaxt = "n", bty = "n", xaxs = "i", ylab = "", xlab = "")
# lines(c(0:100), -0.1 - 0.2 / (1 + exp(-0.2 * (c(0:100) - 20))), col = grey(0.8), lwd = 3)
# u <- par('usr')
# arrows(x0 = u[1], x1 = u[2], y0 = 0, y1 = 0)
# segments(x0 = u[1], x1 = u[1], y0 = u[3], y1 = u[4], xpd = NA)
# abline(v = 20, lty = 2)
# abline(h= -0.3, lty = 3)
# abline(h= -0.1, lty = 3)
# axis(side = 2, at = c(0, -0.1), labels = c(0, expression(beta[0])), las = 1)
# text(26, -0.32, expression(beta[2] == 0.8), xpd = NA)
# text(32, -0.26, expression(beta[2] == 0.2), col = grey(0.6), xpd = NA)
# arrows(x0 = 50, x1 = 50, y0 = -0.3, y1 = -0.1, length = 0.08, code = 3)
# text(50, -0.2, pos = 4, expression(beta[1]))
# text(20, 0, pos = 3, expression(beta[3]), xpd = NA)
# mtext(side = 1, "Habitat pressure value", line = 2)
# mtext(side = 2, "Trend in spawner abundance", line = 2)
# 
# k <- 0.3
# m <- 20
# plot(c(0:100), exp(k * (c(0:100) - m)) / (1 + exp(k * (c(0:100) - m))), "l", lwd = 3, bty = "n", xaxs = "i")
# 
# 
# y <- a  + b*x
# 
# p <- exp(y)/(1 + exp(y))
# p <- exp(a  + b*x)/(1 + exp(a  + b*x))

###############################################################################
# Linear threshold model
###############################################################################

model_thresh <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
			}
		
	}
	
	for(i in 1:nSpawn){
		# For habitat indicators that are percents (0 - 100)
		for(j in c(1,2,3,8,9,10)){
			beta2[i, j] ~ dunif(0, 100)
		}
		
		# For habitat indicators that are densities
		for(j in c(4:7)){
			beta2[i, j] ~ dlnorm(-1, pow(1, -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		
		for(j in 1:nHab){
			habEffect[i, j] <- (habPressures[i, j] >= beta2[spawnEco[i], j]) * (beta1[spawnEco[i], j] +  thetaFAZ[faz[i], j]) * (habPressures[i, j] - beta2[spawnEco[i], j])
		}
		
		trendPred[i] <- beta0	+ thetaMAZ[maz[i], rearEco[i]] + sum(habEffect[i, ])
		
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * weight[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * weight[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model

###############################################################################
# Basic model - no stream order
###############################################################################

model_basic_noSO <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
		}
		
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		
		for(j in 1:nHab){
			slopeBeta[i, j] <- beta1[spawnEco[i], j] +  thetaFAZ[faz[i], j]
		}
		
		trendPred[i] <- beta0	+ thetaMAZ[maz[i], rearEco[i]] + inprod(slopeBeta[i, ], t(habPressures[i, ]))
		
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * weight[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * weight[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model

###############################################################################
# Basic model - simpler StreamOrder but species dependent
###############################################################################

model_11 <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
			phi[i, j] ~ dnorm(0, pow(3, -2))
		}
		
		
		
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		
		for(j in 1:nHab){
			slopeBeta[i, j] <- beta1[spawnEco[i], j] + phi[spawnEco[i], j] * streamOrder[i] + thetaFAZ[faz[i], j]
		}
		
		trendPred[i] <- beta0	+ thetaMAZ[maz[i], rearEco[i]] + inprod(slopeBeta[i, ], t(habPressures[i, ]))
		
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * weight[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * weight[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model

###############################################################################
# Basic model - simpler StreamOrder
###############################################################################

model_10 <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
		}
		
		phi[j] ~ dnorm(0, pow(3, -2))

	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		
		for(j in 1:nHab){
			slopeBeta[i, j] <- beta1[spawnEco[i], j] + phi[j] * streamOrder[i] + thetaFAZ[faz[i], j]
		}
		
		trendPred[i] <- beta0	+ thetaMAZ[maz[i], rearEco[i]] + inprod(slopeBeta[i, ], t(habPressures[i, ]))
		
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * weight[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * weight[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model



###############################################################################
# Basic model
###############################################################################

model_basic <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
		}
		
		for(i in 1:2){
			phi[i, j] ~ dnorm(0, pow(3, -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		
		for(j in 1:nHab){
			slopeBeta[i, j] <- beta1[spawnEco[i], j] + phi[1, j] * (streamOrder[i] - phi[2, j])^2 + thetaFAZ[faz[i], j]
		}
		
		trendPred[i] <- beta0	+ thetaMAZ[maz[i], rearEco[i]] + inprod(slopeBeta[i, ], t(habPressures[i, ]))
		
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * weight[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * weight[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model


###############################################################################
# Zero-inflated Model with NO stream order 
###############################################################################

model_logistic3par_noStreamOrder <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(-0.2, pow(1, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))

	# for(j in 1:nHab){
	# 	# Standard deviation among FAZ's (different for each habitat variable)
	# 	sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
	# }
	
	# SLopes: 3 parameters (beta1 = c(max, steepness, midpoint))
	for(i in 1:nSpawn){
		for(j in 1:nHab){
			beta1[1, i, j] ~ dnorm(0, pow(3, -2))
			beta1[2, i, j] <- 0.8#~ dlnorm(-1, pow(1, -2))
			
		}
		# For habitat indicators that are percents (0 - 100)
		for(j in c(1,2,3,8,9,10)){
			beta1[3, i, j] ~ dunif(0, 100)
		}
		
		# For habitat indicators that are densities
		for(j in c(4:7)){
			beta1[3, i, j] ~ dlnorm(-1, pow(1, -2))
		}
	
		}
	
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}

	# for(j in 1:nHab){
	# 	for(k in 1:nFAZ){
	# 		thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
	# 	}
	# }
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		
		for(j in 1:nHab){
		habEffects[i, j] <- (beta1[1, spawnEco[i], j])/(1 + exp(- beta1[2, spawnEco[i], j] * (habPressures[i, j] - beta1[3, spawnEco[i], j]))) #+ thetaFAZ[faz[i], j]
		}
		
		trendPred[i] <- beta0	 + thetaMAZ[maz[i], rearEco[i]] + sum(habEffects[i, ])
		
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * weight[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * weight[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model


###############################################################################
# Zero-inflated Model with NO stream order 
###############################################################################

model_zeroInf4 <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
			gamma1[i,j] ~ dnorm(0, pow(3, -2))
		}
		
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		
		for(j in 1:nHab){
			slopeGamma[i, j] <- gamma1[spawnEco[i], j] + thetaFAZ[faz[i], j]
			slopeBeta[i, j] <- beta1[spawnEco[i], j] 
		}
		
		trendPred[i] <- beta0	+ thetaMAZ[maz[i], rearEco[i]] + inprod(slopeGamma[i, ], t(habZero[i, ])) + inprod(slopeBeta[i, ], t(habPressures[i, ]))
		
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * weight[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * weight[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model

###############################################################################
# Zero-inflated Model with continuous stream order but affecting slopes 
###############################################################################

model_zeroInfStreamOrder2 <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
			gamma1[i,j] ~ dnorm(0, pow(3, -2))
		}
		
		for(i in 1:2){
			phi[i, j] ~ dnorm(0, pow(3, -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		
		for(j in 1:nHab){
			slopeGamma[i, j] <- gamma1[spawnEco[i], j] 
			slopeBeta[i, j] <- beta1[spawnEco[i], j] + phi[1, j] * (streamOrder[i] - phi[2, j])^2 + thetaFAZ[faz[i], j]
		}
		
		trendPred[i] <- beta0	+ thetaMAZ[maz[i], rearEco[i]] + inprod(slopeGamma[i, ], t(habZero[i, ])) + inprod(slopeBeta[i, ], t(habPressures[i, ]))
		
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * weight[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * weight[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model


###############################################################################
# Zero-inflated Model with NO stream order 
###############################################################################

model_zeroInf3 <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
			gamma1[i,j] ~ dnorm(0, pow(3, -2))
		}
		
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		
		for(j in 1:nHab){
			slopeGamma[i, j] <- gamma1[spawnEco[i], j] 
			slopeBeta[i, j] <- beta1[spawnEco[i], j] + thetaFAZ[faz[i], j]
		}
		
		trendPred[i] <- beta0	+ thetaMAZ[maz[i], rearEco[i]] + inprod(slopeGamma[i, ], t(habZero[i, ])) + inprod(slopeBeta[i, ], t(habPressures[i, ]))
		
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * weight[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * weight[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model

###############################################################################
# Zero-inflated Model with continuous stream order but affecting slopes 
###############################################################################

model_zeroInfStreamOrder2 <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
			gamma1[i,j] ~ dnorm(0, pow(3, -2))
		}
		
		for(i in 1:2){
			phi[i, j] ~ dnorm(0, pow(3, -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		
		for(j in 1:nHab){
			slopeGamma[i, j] <- gamma1[spawnEco[i], j] 
			slopeBeta[i, j] <- beta1[spawnEco[i], j] + phi[1, j] * (streamOrder[i] - phi[2, j])^2 + thetaFAZ[faz[i], j]
		}
		
		trendPred[i] <- beta0	+ thetaMAZ[maz[i], rearEco[i]] + inprod(slopeGamma[i, ], t(habZero[i, ])) + inprod(slopeBeta[i, ], t(habPressures[i, ]))
		
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * weight[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * weight[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model

###############################################################################
# Zero-inflated Model with continuous stream order
###############################################################################

model_zeroInfStreamOrder <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
			gamma1[i,j] ~ dnorm(0, pow(3, -2))
		}
		
		for(i in 1:2){
			phi[i, j] ~ dnorm(0, pow(3, -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		
		for(j in 1:nHab){
			slopeGamma[i, j] <- gamma1[spawnEco[i], j] + phi[1, j] * (streamOrder[i] - phi[2, j])^2
			slopeBeta[i, j] <- beta1[spawnEco[i], j] + thetaFAZ[faz[i], j]
		}
		
		trendPred[i] <- beta0	+ thetaMAZ[maz[i], rearEco[i]] + inprod(slopeGamma[i, ], t(habZero[i, ])) + inprod(slopeBeta[i, ], t(habPressures[i, ]))
	
		}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * weight[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * weight[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model

###############################################################################
# Zero-inflated Model with 2-way interaction only, MAZ and new spawn/rear ecotypes
###############################################################################

model_2wayzeroInf <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
			gamma1[i,j] ~ dnorm(0, pow(3, -2))
		}
		
		# for(i in 1:2){
		# 	beta2[i, j] ~ dnorm(0, pow(3, -2))
		# }
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		for(j in 1:nHab){
			slope[i, j] <- beta1[spawnEco[i], j] + thetaFAZ[faz[i], j] #+ beta2[1, j] * streamSmall[i] + beta2[2, j] * streamLarge[i] 
		}
		trendPred[i] <- beta0 + thetaMAZ[maz[i], rearEco[i]] + inprod(gamma1[spawnEco[i], ], t(habZero[i, ])) + inprod(slope[i, ], t(habPressures[i, ]))
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * dataQuality[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * dataQuality[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model

###############################################################################
# Full model linking spawner data
###############################################################################

model_full <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for trend in spawner abundance
	for(i in 1:nPop){
		alpha0[i] ~ dnorm(0, pow(3, -2))
	}
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	# sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		# sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
		}
		
		# for(i in 1:2){
		# 	beta2[i, j] ~ dnorm(0, pow(3, -2))
		# }
	}
	
	# #-----------------------------------------------------------------------------
	# # Random effects
	# for(j in 1:nMAZ){
	# 	for(k in 1:nRear){
	# 		thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
	# 	}}
	# 
	# for(j in 1:nHab){
	# 	for(k in 1:nFAZ){
	# 		thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
	# 	}
	# }
	

	#-----------------------------------------------------------------------------
	# Calculate predicted trend
	#-----------------------------------------------------------------------------
	for(i in 1:nPop){
		for(t in 1:spawner_nt[i]){
			logSp1_pred[i,t] <- alpha0[i] + alpha1[i] * analysis_yr[i, t]
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	#-----------------------------------------------------------------------------
	for(i in 1:nPop){
		for(j in 1:nHab){
			slope[i, j] <- beta1[spawnEco[i], j] #+ thetaFAZ[faz[i], j]+ beta2[1, j] * streamSmall[i] + beta2[2, j] * streamLarge[i] 
		}
		alpha1[i] <- beta0 + inprod(slope[i, ], t(habPressures[i, ]))#+ thetaMAZ[maz[i], rearEco[i]] 
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	#-----------------------------------------------------------------------------
	
	for(i in 1:nPop){
		for(t in 1:spawner_nt[i]){
			spawners[i, t] ~ dnorm(logSp1_pred[i, t], pow(sigma, -2))# * dataQuality[i])
		}
	}
	
}

###############################################################################
# Model with 2-way interaction only, MAZ and new spawn/rear ecotypes
###############################################################################

model_2way2 <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	beta0 ~ dnorm(0, pow(3, -2))
	
	# Standard deviation among coastal regions
	sigmaMAZ ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nSpawn){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
		}
		
	# 	for(i in 1:2){
	# 		beta2[i, j] ~ dnorm(0, pow(3, -2))
	# 	}
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nMAZ){
		for(k in 1:nRear){
			thetaMAZ[j, k] ~ dnorm(0, pow(sigmaMAZ, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		for(j in 1:nHab){
			slope[i, j] <- beta1[spawnEco[i], j]  + thetaFAZ[faz[i], j] #+ beta2[1, j] * streamSmall[i] + beta2[2, j] * streamLarge[i]
		}
		trendPred[i] <- beta0 + thetaMAZ[maz[i], rearEco[i]] + inprod(slope[i, ], t(habPressures[i, ]))
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * dataQuality[i])
		
		# Track log-likelihood for model evaluation
		log_lik0[i] <- logdensity.norm(trend[i], trendPred[i], pow(sigma, -2) * dataQuality[i])
	}
	
	log_lik <- sum(log_lik0)
	
} # end model

###############################################################################
# Model with 2-way interaction only, not 3-way
###############################################################################

model_2way_streamOrdernoRanef <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	for(s in 1:nSpecies){
		beta0[s] ~ dnorm(0, pow(3, -2))
	}
	
	for(j in 1:nHab){
		# Strength of popHab relationship
		for(i in 1:nFWeco){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
		}
		
		# Stream order: base level is middle sized
		# Estimate difference for small and large systems
		beta2[1, j] ~ dnorm(0, pow(3, -2))
		beta2[2, j] ~ dnorm(0, pow(3, -2))
	}
	
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		for(j in 1:nHab){
			slope[i, j] <- beta1[FWeco[i], j] + beta2[1, j] * streamSmall[i] + beta2[2, j] * streamLarge[i]
		}
		trendPred[i] <- beta0[species[i]] + inprod(slope[i, ], t(habPressures[i, ]))
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2))
	}
} # end model

###############################################################################
# Model with 2-way interaction only, not 3-way
###############################################################################

model_2way_noQual <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	for(s in 1:nSpecies){
		beta0[s] ~ dnorm(0, pow(3, -2))
	}
	
	# Standard deviation among coastal regions
	sigmaRegion ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nFWeco){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
		}
		
		for(i in 1:2){
			beta2[i, j] ~ dnorm(0, pow(3, -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nRegion){
		for(k in 1:nSpecies){
			thetaRegion[j, k] ~ dnorm(0, pow(sigmaRegion, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		for(j in 1:nHab){
			slope[i, j] <- beta1[FWeco[i], j] + beta2[1, j] * streamSmall[i] + beta2[2, j] * streamLarge[i] + thetaFAZ[faz[i], j]
		}
		trendPred[i] <- beta0[species[i]] + thetaRegion[region[i], species[i]] + inprod(slope[i, ], t(habPressures[i, ]))
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2))
	}
} # end model

###############################################################################
# Model with 2-way interaction only
# NO stream order interaction
# No data quality
###############################################################################

model_2way_FAZ <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	for(s in 1:nSpecies){
		beta0[s] ~ dnorm(0, pow(3, -2))
	}
	
	# Standard deviation among coastal regions
	sigmaRegion ~ dlnorm(-1, pow(1, -2))
	
	# Slope for habitat predictor
	for(j in 1:nHab){
		
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		for(i in 1:nFWeco){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nRegion){
		for(k in 1:nSpecies){
			thetaRegion[j, k] ~ dnorm(0, pow(sigmaRegion, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		
		for(j in 1:nHab){
			slope[i, j] <- (beta1[FWeco[i], j] + thetaFAZ[faz[i], j])
		}
		
		trendPred[i] <- beta0[species[i]] + thetaRegion[region[i], species[i]] + inprod(slope[i, ], t(habPressures[i, ]))
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2))
	}
} # end model

###############################################################################
# Model with 2-way interaction only
# NO FAZ random effect
# NO stream order interaction
# No data quality
###############################################################################

model_2way_region <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	for(s in 1:nSpecies){
		beta0[s] ~ dnorm(0, pow(3, -2))
	}
	
	# Standard deviation among coastal regions
	sigmaRegion ~ dlnorm(-1, pow(1, -2))
	
	# Slope for habitat predictor
	for(j in 1:nHab){
		for(i in 1:nFWeco){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nRegion){
		for(k in 1:nSpecies){
			thetaRegion[j, k] ~ dnorm(0, pow(sigmaRegion, -2))
		}}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		trendPred[i] <- beta0[species[i]] + thetaRegion[region[i], species[i]] + inprod(beta1[FWeco[i], ], t(habPressures[i, ]))
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2))
	}
} # end model

###############################################################################
# Model with 2-way interaction only, not 3-way
# No random effects
# No stream order interaction
# No data quality weight
###############################################################################

model_2way_basic <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	for(s in 1:nSpecies){
		beta0[s] ~ dnorm(0, pow(3, -2))
	}
	
	# Strength of popHab relationship
	for(j in 1:nHab){
		for(i in 1:nFWeco){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		trendPred[i] <- beta0[species[i]] + inprod(beta1[FWeco[i], ], t(habPressures[i, ]))
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2))
	}
} # end model

###############################################################################
# Model with 2-way interaction only, not 3-way
###############################################################################

model_2way <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	for(s in 1:nSpecies){
		beta0[s] ~ dnorm(0, pow(3, -2))
	}
	
	# Standard deviation among coastal regions
	sigmaRegion ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nFWeco){
			beta1[i, j] ~ dnorm(0, pow(3, -2))
		}
		
		for(i in 1:2){
			beta2[i, j] ~ dnorm(0, pow(3, -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nRegion){
		for(k in 1:nSpecies){
			thetaRegion[j, k] ~ dnorm(0, pow(sigmaRegion, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dnorm(0, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		for(j in 1:nHab){
			slope[i, j] <- beta1[FWeco[i], j] + beta2[1, j] * streamSmall[i] + beta2[2, j] * streamLarge[i] + thetaFAZ[faz[i], j]
		}
		trendPred[i] <- beta0[species[i]] + thetaRegion[region[i], species[i]] + inprod(slope[i, ], t(habPressures[i, ]))
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * dataQuality[i])
	}
} # end model
###############################################################################
# Model with 3-way interaction
###############################################################################

model_3way <- function(){
	
	#-----------------------------------------------------------------------------
	# Priors on parameters
	sigma ~ dlnorm(-1, pow(1, -2))
	
	# Intercept for pop-hab relationship
	# Different for 6 marine ecotypes: 
	# ChinookOcean, ChinookStream, Chum, Coho, Pink, Sockeye 
	for(s in 1:nSpecies){
		beta0[s] ~ dnorm(0, pow(3, -2))
	}
	
	# Standard deviation among coastal regions
	sigmaRegion ~ dlnorm(-1, pow(1, -2))
	
	for(j in 1:nHab){
		# Standard deviation among FAZ's (different for each habitat variable)
		sigmaFAZ[j] ~ dlnorm(-1, pow(1, -2))
		
		# Strength of popHab relationship
		# Different for each habitat pressure, FW ecotype (7), and stream order.
		for(i in 1:nFWeco){
			for(k in 1:nStreamOrder){
				beta1[j, i, k] ~ dnorm(0, pow(3, -2))
			}
		}
	}
	
	#-----------------------------------------------------------------------------
	# Random effects
	for(j in 1:nRegion){
		for(k in 1:nSpecies){
			thetaRegion[j, k] ~ dnorm(-1, pow(sigmaRegion, -2))
		}}
	
	for(j in 1:nHab){
		for(k in 1:nFAZ){
			thetaFAZ[k, j] ~ dlnorm(-1, pow(sigmaFAZ[j], -2))
		}
	}
	
	#-----------------------------------------------------------------------------
	# Link to habitat data
	# Intercept
	for(i in 1:nPop){
		for(j in 1:nHab){
			slope[i, j] <- (beta1[j, FWeco[i], streamOrder[i]] + thetaFAZ[faz[i], j])
		}
		trendPred[i] <- beta0[species[i]] + thetaRegion[region[i], species[i]] + slope[i, 1:nHab] %*% t(habPressures[i, 1:nHab])
		# trendPred[i] <- beta0[Meco[i]] + thetaRegion[region[i]] + t(beta1[1:nHab]) %*% habPressures[i, 1:nHab]
	}
	
	#-----------------------------------------------------------------------------
	# Likelihood
	for(i in 1:nPop){
		trend[i] ~ dnorm(trendPred[i], pow(sigma, -2) * dataQuality[i])
	}
} # end model


