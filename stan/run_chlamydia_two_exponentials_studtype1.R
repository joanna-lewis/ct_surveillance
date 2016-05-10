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
	studnum = 4,
	studobs = c(4,5,1,1),
	cumobs = c(0,cumsum(c(4,5,1))),
	Nobs = sum(4,5,1,1),
	r = c(10,7,6,6, 
		2,7,1,0,3, 
		23,
		3),
	n = c(
		23,14,14,8,
		12,28,4,8,6, 
		129, 
		15), 
	t = c(
		0.038,0.058,0.077,0.125, 
		0.012,0.03,0.049,0.088,0.274,  
		0.045, 
		0.083), 
	seind = c(
		1,1,1,1, 
		0,0,0,0,0, 
		0,
		1), 
	T = c( 
		999,999,999,999, 
		999,999,999,999,999, 
		999,
		999)
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

fit <- stan(file = 'chlamydia_two_exponentials_studtype1.stan', data = chlamydia_dat, 
            iter = 100000, warmup = 50000, init=list(init0), chains = 1, seed=12345)

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
centiles <- matrix(nrow=100000, ncol=3)
mean <- matrix(100000, ncol=1)
for(i in 1:100000){
	centiles[i,] <- quantile(op$lambda_slow[1:i], p=c(0.025, 0.5, 0.975), na.rm=TRUE)
	mean[i] <- mean(op$lambda_slow[1:i], na.rm=TRUE)
}
plot(centiles[,1], type='l', ylim=c(0, 10))
lines(centiles[,2])
lines(centiles[,3])
lines(mean, col='red')

abline(h=0.61, lty=2) 
abline(h=0.74, lty=2)
abline(h=0.89, lty=2) 

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
