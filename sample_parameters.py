####################################################################
# This script samples model parameters, in the same way as in the 
# notebook england.ipynb.
####################################################################

###################################
# input prior distributions for various things
###################################

from numpy import *
from scipy.stats import beta
from scipy.stats import gamma
from numpy.random import normal
from scipy.optimize import fsolve

n_sample = 10000

######################
# parameters of beta distributions representing the proportion of the population sexually 
# active, by sex and age group
######################
# men, 16-24
[alpha_m_16_24, beta_m_16_24] = fsolve(
    lambda x: array(beta.interval(0.95, x[0], x[1], loc=0, scale=1))
    - (0.8023836019, 0.843403825),
    [1,1]
    )
# women, 16-24
[alpha_f_16_24, beta_f_16_24] = fsolve(
    lambda x: array(beta.interval(0.95, x[0], x[1], loc=0, scale=1))
    - (0.7998634469, 0.837979601),
    [1,1]
    )

######################
# sexually-active population:
######################

random.seed = 12345

# Population, testing and diagnosis data is from http://www.chlamydiascreening.nhs.uk/ps/data.asp (downloaded 17 April 2015).
p_active_m_16_24 = random.beta(alpha_m_16_24, beta_m_16_24, size=n_sample) # 16-24 yo only
pop_active_m_15_24 = random.binomial(3519015, p_active_m_16_24, size=n_sample)
p_active_f_16_24 = random.beta(alpha_f_16_24, beta_f_16_24, size=n_sample) # 16-24 yo only
pop_active_f_15_24 = random.binomial(3388842, p_active_f_16_24, size=n_sample)

######################
# testing and diagnosis rates, per person per year:
######################
# Population, testing and diagnosis data is from http://www.chlamydiascreening.nhs.uk/ps/data.asp (downloaded 17 April 2015).
test_rate_m_15_24 = random.gamma(566908, 1, size=n_sample)/pop_active_m_15_24
diag_rate_m_15_24 = random.gamma(48387, 1, size=n_sample)/pop_active_m_15_24
test_rate_f_15_24 = random.gamma(1205896, 1, size=n_sample)/pop_active_f_15_24
diag_rate_f_15_24 = random.gamma(88101, 1, size=n_sample)/pop_active_f_15_24


# Proportion of incident infections asymptomatic is not known, so we use estimates of the 
# proportion of _prevalent_ infections asymptomatic. Note that these provide an upper bound
# because treatment seeking in symptomatic cases will deplete the symptomatic pool.
p_asymp_m = random.beta(69 + 1, 78 - 69 + 1, size=n_sample) # McKay et al. Lancet (2003) 88% 
p_asymp_f = random.beta(135 + 1, 26 + 1, size=n_sample) # Kahn et al. Sex Transm Dis (2003) 84% NB numbers taken from text, p656

######################
# test performance
######################
p_true_pos_m = random.beta(32+1, 0+1, size=n_sample) # Horner J. Clin. Microbiol (2005): 32 of 32 infected samples tested +ve
p_false_pos_m = random.beta(2+1, 950+1, size=n_sample) # Horner J. Clin. Microbiol (2005): 2 of 952 uninfected samples tested +ve
p_true_pos_f = random.beta(129+1, 12+1, size=n_sample) # Low Health Technol Assess (2007): 129 of 141 infected samples tested +ve
p_false_pos_f = random.beta(4+1, 2323+1, size=n_sample) # Low Health Technol Assess (2007): 4 of 2327 uninfected samples tested +ve

###################################
# Rate of treatment seeking by symptomatic cases
###################################

# Find beta distributions corresponding to 95% CIs reported in Mercer Sex. Transm. Infect. (2007) (see table above).

from numpy import *
from scipy.stats import beta
from scipy.optimize import fsolve

a = empty(5)
b = empty(5)

# < 1 week
[a[0], b[0]] = fsolve(
    lambda x: array(beta.interval(0.95, x[0], x[1], loc=0, scale=1))
    - (0.144, 0.442),
    [1,1]
    )

# 7-13 days
[a[1], b[1]] = fsolve(
    lambda x: array(beta.interval(0.95, x[0], x[1], loc=0, scale=1))
    - (0.061, 0.302),
    [1,1]
    )

# 14-27 days
[a[2], b[2]] = fsolve(
    lambda x: array(beta.interval(0.95, x[0], x[1], loc=0, scale=1))
    - (0.133, 0.310),
    [1,1]
    )

# 28-41 days
[a[3], b[3]] = fsolve(
    lambda x: array(beta.interval(0.95, x[0], x[1], loc=0, scale=1))
    - (0.085, 0.299),
    [1,1]
    )

# 42 days and over
[a[4], b[4]] = fsolve(
    lambda x: array(beta.interval(0.95, x[0], x[1], loc=0, scale=1))
    - (0.055, 0.564),
    [1,1]
    )

##############
# Metropolis-Hastings to get a sample for rate of treatment
##############

i = 0
att_symp = empty(n_sample+1000) # testing rate per person per year. Allow 1000 extra samples for burn-in
ll = empty(n_sample+1000) # log-likelihood
props = empty([n_sample+1000, 5]) # simulated data, for posterior predictive check
old = 0.04 # starting sample value
new = 0.04 # starting sample value

# simulate probabilities corresponding to data

# proportion expected in each time window
tps = array([0., 7., 14., 28., 42., Inf])
simp_old = exp(-old*tps[:5]) - exp(-old*tps[1:])
simp_new = exp(-new*tps[:5]) - exp(-new*tps[1:])

acc=0.
while i < n_sample+1000: # to do samples for p_test_symp
    
    new = random.normal(old, 0.05) # generate a sample from normal distribution
    
    if new < 0:
        att_symp[i] = old # reject
        ll[i] = -1e10
    else:
        simp_old = exp(-old*tps[:5]) - exp(-old*tps[1:])
        simp_new = exp(-new*tps[:5]) - exp(-new*tps[1:])

        if sum(simp_new > 0) != len(tps) - 1:
            att_symp[i] = old # reject
            ll[i] = -1e10
        else:
            # simulate probabilities corresponding to the data
            log_ratio = sum(beta.logpdf(simp_new, a, b, loc=0, scale=1)) - sum(beta.logpdf(simp_old, a, b, loc=0, scale=1))
    
            if log(random.uniform(0,1)) <  log_ratio:
                att_symp[i] = new # accept
                ll[i] = sum(beta.logpdf(simp_new, a, b, loc=0, scale=1))
                old = new
                acc = acc+1
            else:
                att_symp[i] = old # reject
                ll[i] = sum(beta.logpdf(simp_old, a, b, loc=0, scale=1))
    
    props[i] = simp_old
    i = i+1
    
att_symp = att_symp[1000:] # remove burn-in samples
ll = ll[1000:] # log-likelihood
    
print acc/(n_sample+1000) # print the proportion of samples accepted
print mean(att_symp)*365.25

att_symp = att_symp*365.25 # convert rate from day^-1 to year^-1

###################################
# Spontaneous clearance
###################################

from scipy.stats import binom

##############
# Metropolis-Hastings to get a sample for rate of spontaneous clearance in men
# assuming a constant hazard of recovery
##############

i = 0
sc_m = empty(n_sample) # testing rate per person per year
ll_m = empty(n_sample) # log-likelihood
old = 0.1 # starting sample value

acc=0.
while i < n_sample: # to do samples for p_test_symp
    
    new = random.normal(old, 5) # generate a sample from normal distribution
    
    if new < 0:
        sc_m[i] = old # reject
        ll_m[i] = -1e10
    else:
        simp_old = 1 - exp(-array([5, 11.5, 18.5, 32.5, 78])*old/365.25)
        simp_new = 1 - exp(-array([5, 11.5, 18.5, 32.5, 78])*new/365.25)
                
        if sum(simp_new >0) != 5:
            sc_m[i] = old # reject
            ll_m[i] = -1e10
        else:
            # simulate probabilities 
            log_ratio = sum(binom.logpmf([3,2,1,0,1], [15,9,4,4,4], simp_new)) - sum(binom.logpmf([3,2,1,0,1], [15,9,4,4,4], simp_old))            
            
            if log(random.uniform(0,1)) <  log_ratio: 
                sc_m[i] = new # accept
                ll_m[i] = sum(binom.logpmf([3,2,1,0,1], [15,9,4,4,4], simp_new))
                old = new
                acc = acc+1
            else:
                sc_m[i] = old # reject
                ll_m[i] = sum(binom.logpmf([3,2,1,0,1], [15,9,4,4,4], simp_old))
        
    i = i+1
    
print acc/n_sample

##############
# Metropolis-Hastings to get a sample for rate of spontaneous clearance in women
# assuming a constant hazard of recovery
##############

i = 0
sc_f = empty(n_sample) # testing rate per person per year
old = 0.1 # starting sample value

acc=0.
while i < n_sample: # to do samples for p_test_symp
    
    new = random.normal(old, 5) # generate a sample from normal distribution
    
    if new < 0:
        sc_f[i] = old # reject
    else:
        simp_old = 1 - exp(-array([5, 11.5, 18.5, 32.5, 137.5])*old/365.25)
        simp_new = 1 - exp(-array([5, 11.5, 18.5, 32.5, 137.5])*new/365.25)
                
        if sum(simp_new >0) != 5:
            sc_f[i] = old # reject
        else:
            # simulate probabilities 
            log_ratio = sum(binom.logpmf([2,7,1,0,3], [12,28,4,8,6], simp_new)) - sum(binom.logpmf([2,7,1,0,3], [12,28,4,8,6], simp_old))
            
            if log(random.uniform(0,1)) <  log_ratio: 
                sc_f[i] = new # accept
                old = new
                acc = acc+1
            else:
                sc_f[i] = old # reject
        
    i = i+1
    
print acc/n_sample