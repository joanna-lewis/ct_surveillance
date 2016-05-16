###################################
# Run STAN model chlamydia_two_exponentials.stan using the R interface
###################################

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
# 9 Mollano

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

fit <- stan(file = 'chlamydia_two_exponentials.stan', data = chlamydia_dat, 
            iter = 22000, warmup = 2000, init=list(init0), chains = 1, seed=12345,
            sample_file = 'chlamydia_two_exponentials_init0_280416')

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

write.csv(
	data.frame(lambda_slow = op$lambda_slow[10001:20000], 
		p_fast = op$p1[10001:20000]
		),
	file='chlamydia_two_exponentials_init0_130516.csv', row.names=FALSE)

# ####################
# # some tests on behaviour with initial parameters
# ####################

# psi = 0.9
# lambda_slow = 0.7
# p1 = 0.2

# theta <- numeric(length=chlamydia_dat$Nobs)
# lambda <- c(120, lambda_slow)

  # for (i in 1:4) { # should find a way of using input data to state study categories
    # for (j in 1:chlamydia_dat$studobs[i]) {
      # theta[chlamydia_dat$cumobs[i]+j] <-  0;
      # for(k in 0:chlamydia_dat$n[chlamydia_dat$cumobs[i]+j]){
        # theta[chlamydia_dat$cumobs[i]+j] <- theta[chlamydia_dat$cumobs[i]+j] + choose(chlamydia_dat$n[chlamydia_dat$cumobs[i]+j],k)  * p1^k * (1-p1)^(chlamydia_dat$n[chlamydia_dat$cumobs[i]+j]-k) *  ((k / (1.0*chlamydia_dat$n[chlamydia_dat$cumobs[i]+j])) * (1 - exp(-lambda[1] * chlamydia_dat$t[chlamydia_dat$cumobs[i]+j])) + (1 - (k / (1.0*chlamydia_dat$n[chlamydia_dat$cumobs[i]+j]))) * (1 - exp(-lambda[2] * chlamydia_dat$t[chlamydia_dat$cumobs[i]+j])));
        # }
        # theta[chlamydia_dat$cumobs[i]+j] <- theta[chlamydia_dat$cumobs[i]+j] / (psi / (1 - ((1-psi) * (chlamydia_dat$seind[chlamydia_dat$cumobs[i]+j] == 0))));
      # }
    # }
   
  # w1 <- (p1 / lambda[1]) / (p1 / lambda[1] + (1 - p1) / lambda[2]);
  # for (i in 5:7) { # should find a way of using input data to state study categories
    # for (j in 1:chlamydia_dat$studobs[i]) {
      # theta[chlamydia_dat$cumobs[i]+j] <- 0;
      # for(k in 0:chlamydia_dat$n[chlamydia_dat$cumobs[i]+j]){
        # theta[chlamydia_dat$cumobs[i]+j] <- theta[chlamydia_dat$cumobs[i]+j] + choose(chlamydia_dat$n[chlamydia_dat$cumobs[i]+j],k) * w1^k * (1-w1)^(chlamydia_dat$n[chlamydia_dat$cumobs[i]+j]-k) *  ((k / (1.0*chlamydia_dat$n[chlamydia_dat$cumobs[i]+j])) * (1 - exp(-lambda[1] * chlamydia_dat$t[chlamydia_dat$cumobs[i]+j])) + (1 - (k / (1.0*chlamydia_dat$n[chlamydia_dat$cumobs[i]+j]))) * (1 - exp(-lambda[2] * chlamydia_dat$t[chlamydia_dat$cumobs[i]+j])));
        # }
        # theta[chlamydia_dat$cumobs[i]+j] <- theta[chlamydia_dat$cumobs[i]+j] / (psi / (1 - ((1-psi) * (chlamydia_dat$seind[chlamydia_dat$cumobs[i]+j] == 0))));
      # }
    # }
    
  # for (i in 8:9) {
    # for (j in 1:chlamydia_dat$studobs[i]) {
        # theta[chlamydia_dat$cumobs[i]+j] <- 0;
        # for (k in 0:chlamydia_dat$n[chlamydia_dat$cumobs[i]+j]){
          
          # pk <- choose(chlamydia_dat$n[chlamydia_dat$cumobs[i]+j],k) * w1^k * (1-w1)^(chlamydia_dat$n[chlamydia_dat$cumobs[i]+j]-k); # temporary variable, containing probability of this k

          # theta[chlamydia_dat$cumobs[i]+j] <- theta[chlamydia_dat$cumobs[i]+j] + pk * (k / (1.0*chlamydia_dat$n[chlamydia_dat$cumobs[i]+1])) * exp(-lambda[1]*chlamydia_dat$T[chlamydia_dat$cumobs[i]+j]) * (1 - exp(-lambda[1]*chlamydia_dat$t[chlamydia_dat$cumobs[i]+j])) / ( (k / (1.0*chlamydia_dat$n[chlamydia_dat$cumobs[i]+1])) * exp(-lambda[1]*chlamydia_dat$T[chlamydia_dat$cumobs[i]+j]) + (1 - k / (1.0*chlamydia_dat$n[chlamydia_dat$cumobs[i]+1])) * exp(-lambda[2]*chlamydia_dat$T[chlamydia_dat$cumobs[i]+j]) );   # fast category
          # theta[chlamydia_dat$cumobs[i]+j] <- theta[chlamydia_dat$cumobs[i]+j] + pk * (1 - k / (1.0*chlamydia_dat$n[chlamydia_dat$cumobs[i]+1])) * exp(-lambda[2]*chlamydia_dat$T[chlamydia_dat$cumobs[i]+j]) * (1 - exp(-lambda[2]*chlamydia_dat$t[chlamydia_dat$cumobs[i]+j])) / ( (k / (1.0*chlamydia_dat$n[chlamydia_dat$cumobs[i]+1])) * exp(-lambda[1]*chlamydia_dat$T[chlamydia_dat$cumobs[i]+j]) + (1 - k / (1.0*chlamydia_dat$n[chlamydia_dat$cumobs[i]+1])) * exp(-lambda[2]*chlamydia_dat$T[chlamydia_dat$cumobs[i]+j]) ) ;  # slow category
          
          # }
      # theta[chlamydia_dat$cumobs[i]+j] <- theta[chlamydia_dat$cumobs[i]+j] / (psi / (1 - ((1-psi) * (chlamydia_dat$seind[chlamydia_dat$cumobs[i]+j] == 0)))); # apply correction for culture tests
    
      # }
    # }
