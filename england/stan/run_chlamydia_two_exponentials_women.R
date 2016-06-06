###################################
# Run STAN model chlamydia_two_exponentials_women.stan using the R interface
# Based on Supplementary Web Appendix 1 to:
# Malcolm J. Price, A. E. Ades, Daniela De Angelis, Nicky J. Welton, John Macleod, Kate Soldan, Katy Turner, Ian Simms and Paddy J. Horner
# Mixture-of-exponentials models to explain heterogeneity in studies of the duration of Chlamydia trachomatis infection.
# Statistics in Medicine 32:1547–1560 (2013)
###################################

rm(list=ls())

library(rstan)

####################
# data
####################
# duration
# study order
# 1 Johhanisson
# 2 Joyner
# 3 Geisler
# 4 Paavonen	
# 5 Rahm
# 6 Sorensen
# 7 McCormack	
# 8 Morre	
# 9 Molano

chlamydia_dat <- list(
	studnum = 9,
	studnum_bytype = c(4,3,2),
	studobs = c(4,5,1,1,3,1,1,5,4),
	cumobs = c(0,cumsum(c(4,5,1,1,3,1,1,5))),
	Nobs = sum(4,5,1,1,3,1,1,5,4),
	r = c(10,7,6,6, 
		2,7,1,0,3, 
		23,
		3, 
		17,0,0,
		8, 
		3, 
		2,2,4,0,2,
		44,23,7,2),
	n = c(
		23,14,14,8,
		12,28,4,8,6, 
		129, 
		15, 
		93,1,1, 
		13, 
		7, 
		20,5,15,1,13, 
		82,37,14,6), 
	t = c(
		0.038,0.058,0.077,0.125, 
		0.012,0.03,0.049,0.088,0.274,  
		0.045, 
		0.083, 
		0.25,0.5,0.75,
		1,
		1.375, 
		0.083,0.5,0.417,0.917,0.5, 
		1,1,1,1), 
	seind = c(
		1,1,1,1, 
		0,0,0,0,0, 
		0,
		1, 
		1,1,1,
		0, 
		1, 
		0,0,0,0,0,
		0,0,0,0), 
	T = c( 
		999,999,999,999, 
		999,999,999,999,999, 
		999,
		999, 
		999,999,999,
		999, 
		999, 
		0,0,0.083,0.083,0.5, 
		0,1,2,3
		)
	)
	
####################
# initial values
####################

# Initial values - 0
init0 <- list(psi = 7/77, lambda_slow = 0.74, p1 = 0.23)

# Initial values - 1
init1 <- list(psi = 0.9, lambda_slow = 0.7, p1 = 0.2)

# Initial values - 2
init2 <- list(psi = 0.6, lambda_slow = 0.1, p1 = 0.5)

####################
# run the model
####################

fit <- stan(file = 'chlamydia_two_exponentials_women.stan', data = chlamydia_dat, 
            iter = 22000, warmup = 2000, init=list(init0, init1, init2), chains = 3, seed=12345
            )

####################
# plots
####################

op <- extract(fit)
plot(op$lp__, type='l', col=rgb(0,0,0,0.1))
pairs(~lp__+lambda_slow+p1+psi,data=op, col=rgb(0,0,0,0.1), pch=16, cex=0.1)
hist(op$lambda_slow)
hist(op$p1)

mean(op$lambda_slow)
quantile(op$lambda_slow, p=c(0.025, 0.5, 0.975))
mean(op$p1)
quantile(op$p1, p=c(0.025, 0.5, 0.975))

# how do estimates of mean and centiles change with increasing number of samples?
centiles <- matrix(nrow=20000, ncol=3)
mean <- matrix(20000, ncol=1)
for(i in 1:20000){
	centiles[i,] <- quantile(op$lambda_slow[1:i], p=c(0.025, 0.5, 0.975), na.rm=TRUE)
	mean[i] <- mean(op$lambda_slow[1:i], na.rm=TRUE)
}
plot(centiles[,1], type='l', ylim=c(0.6, 1))
lines(centiles[,2])
lines(centiles[,3])
lines(mean, col='red')

abline(h=0.61, lty=2) 
abline(h=0.74, lty=2)
abline(h=0.89, lty=2) 

####################
# write 10000 samples to csv file (choose last 3333 or 3334 samples from each chain)
####################

write.csv(
	data.frame(lambda_slow = op$lambda_slow[c(16668:20000, 36668:40000, 56667:60000)], 
		p_fast = op$p1[c(16668:20000, 36668:40000, 56667:60000)]
		),
	file='chlamydia_two_exponentials_women.csv', row.names=FALSE)

