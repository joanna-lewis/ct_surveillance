###################################
# Run STAN model chlamydia_two_exponentials_men.stan using the R interface
# Model based on:
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
# 1 Handsfield 1976
# 2 Prentice 1976
# 3 Johannisson 1979
# 4 Paavonen 1980
# 5 Joyner 2002
# 6 Geisler 2008
# 7 Stamm 1986
# 8 van den Brule 2002	

chlamydia_dat <- list(
	studnum = 8,
	studnum_bytype = c(6,2,0),
	studobs = c(1,1,4,1,5,1,1,4),
	cumobs = c(0,cumsum(c(1,1,4,1,5,1,1))),
	Nobs = sum(c(1,1,4,1,5,1,1,4)),
	r = c( # number who cleared infection at each time point
		0,
		4, 
		3,13,3,2,
		7,
		3,2,1,0,1,
		5,
		1,0,0,0,
		1
		),
	n = c( # number tested at each time point
		10,
		13,
		17,27,6,2,
		21,
		15, 9, 4, 4, 4,
		14,
		5, 2, 2, 1,
		9
		), 
	t = c(
		0.019,
		0.023, 
		0.019, 0.038, 0.058, 0.077, 
		0.077,
		0.012, 0.030, 0.049, 0.088, 0.190,
		0.045,
		0.019, 0.038, 0.058, 0.077,
		0.5
		), 
	seind = c( # did the study use culture as opposed to NAATs?
		1,
		1,
		1, 1, 1, 1,
		1,
		0,0,0,0,0,
		0,
		1,1,1,1,
		0
		), 
	T = rep(999, times=sum(c(1,1,4,1,5,1,1,4)))
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

fit <- stan(file = 'chlamydia_two_exponentials_men.stan', data = chlamydia_dat, 
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
plot(centiles[,1], type='l', ylim=c(0, 1.5))
lines(centiles[,2])
lines(centiles[,3])
lines(mean, col='red')

abline(h=0.61, lty=2) 
abline(h=0.74, lty=2)
abline(h=0.89, lty=2) 

# posterior predictive check
c2.5 <- apply(op$r_sim,2,quantile,p=0.025)/chlamydia_dat$n
c97.5 <- apply(op$r_sim,2,quantile,p=0.975)/chlamydia_dat$n
study_names <- c("Handsfield 1976", "Prentice 1976", "Johannisson 1979", "Paavonen 1980", "Joyner 2002", "Geisler 2008", "Stamm 1986", "van den Brule 2002")
plot(
	chlamydia_dat$t, 
	chlamydia_dat$r/chlamydia_dat$n, 
	pch=16, 
	col=rep(1:8, times = chlamydia_dat$studobs), 
	log='x',
	xlab = 'Mean follow-up time (years)', 
	ylab= 'Proportion cleared infection'
	)
arrows(
	chlamydia_dat$t, 
	y0 = apply(op$r_sim,2,quantile,p=0.025)/chlamydia_dat$n, 
	y1 = apply(op$r_sim,2,quantile,p=0.975)/chlamydia_dat$n, 
	col = rep(1:8, times = chlamydia_dat$studobs), 
	code=3, angle=90, length = 0.1
	)
legend('topleft', col=1:8, pch=16, legend=study_names)

####################
# write to csv file (choose last 3333 or 3334 samples from each chain)
####################

write.csv(
	data.frame(lambda_slow = op$lambda_slow[c(16668:20000, 36668:40000, 56667:60000)], 
		p_fast = op$p1[c(16668:20000, 36668:40000, 56667:60000)]
		),
	file='chlamydia_two_exponentials_men.csv', row.names=FALSE)

