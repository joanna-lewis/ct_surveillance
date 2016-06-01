####################################################################
# This script samples model parameters, in the same way as in the 
# notebook england.ipynb.
####################################################################

###################################
# input prior distributions for various things
###################################

import numpy as np
from numpy import *
from scipy.stats import beta
from scipy.stats import gamma
from numpy.random import normal
from scipy.optimize import fsolve

rs = random.RandomState(12345)
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

# Population, testing and diagnosis data is from http://www.chlamydiascreening.nhs.uk/ps/data.asp (downloaded 17 April 2015).
p_active_m_16_24 = rs.beta(alpha_m_16_24, beta_m_16_24, size=n_sample) # 16-24 yo only
pop_active_m_15_24 = rs.binomial(3519015, p_active_m_16_24, size=n_sample)
p_active_f_16_24 = rs.beta(alpha_f_16_24, beta_f_16_24, size=n_sample) # 16-24 yo only
pop_active_f_15_24 = rs.binomial(3388842, p_active_f_16_24, size=n_sample)

######################
# testing and diagnosis rates, per person per year:
######################
# Population, testing and diagnosis data is from http://www.chlamydiascreening.nhs.uk/ps/data.asp (downloaded 17 April 2015).
test_rate_m_15_24 = rs.gamma(566908, 1, size=n_sample)/pop_active_m_15_24
diag_rate_m_15_24 = rs.gamma(48387, 1, size=n_sample)/pop_active_m_15_24
test_rate_f_15_24 = rs.gamma(1205896, 1, size=n_sample)/pop_active_f_15_24
diag_rate_f_15_24 = rs.gamma(88101, 1, size=n_sample)/pop_active_f_15_24

######################
# test performance
######################
p_true_pos_m = rs.beta(32+1, 0+1, size=n_sample) # Horner J. Clin. Microbiol (2005): 32 of 32 infected samples tested +ve
p_false_pos_m = rs.beta(2+1, 950+1, size=n_sample) # Horner J. Clin. Microbiol (2005): 2 of 952 uninfected samples tested +ve
p_true_pos_f = rs.beta(129+1, 12+1, size=n_sample) # Low Health Technol Assess (2007): 129 of 141 infected samples tested +ve
p_false_pos_f = rs.beta(4+1, 2323+1, size=n_sample) # Low Health Technol Assess (2007): 4 of 2327 uninfected samples tested +ve

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
    
    new = rs.normal(old, 0.05) # generate a sample from normal distribution
    
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
    
            if log(rs.uniform(0,1)) <  log_ratio:
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
    
#print acc/(n_sample+1000) # print the proportion of samples accepted
#print mean(att_symp)*365.25

att_symp = att_symp*365.25 # convert rate from day^-1 to year^-1

###################################
# Spontaneous clearance
###################################

import csv
sc_m = empty(n_sample) # clearance rate per person per year
with open('../england/stan/chlamydia_two_exponentials_men.csv', 'rU') as m:
    reader = csv.reader(m)
    i=0
    next(reader) # skip the header row
    for row in reader:
        sc_m[i] = row[0]
        i = i+1
        
sc_f = empty(n_sample) # clearance rate per person per year
with open('../england/stan/chlamydia_two_exponentials_women.csv', 'rU') as f:
    reader = csv.reader(f)
    i=0
    next(reader) # skip the header row
    for row in reader:
        sc_f[i] = row[0]
        i = i+1

###################################
# Infer the proportion of incident infections which are asymptomatic, by calibrating to the Natsal-3 prevalence estimates.
###################################

from scipy.stats import beta

# generate samples for prevalence

[alpha_prev_m, beta_prev_m] = fsolve(
    lambda x: array(beta.interval(0.95, x[0], x[1], loc=0, scale=1))
    - (0.015, 0.034), # Natsal-3 prevalence in men
    [1,1]
    )
prev_m = rs.beta(alpha_prev_m, beta_prev_m, size=n_sample)

[alpha_prev_f, beta_prev_f] = fsolve(
    lambda x: array(beta.interval(0.95, x[0], x[1], loc=0, scale=1))
    - (0.022, 0.043), # Natsal-3 prevalence in women
    [1,1]
    )
prev_f = rs.beta(alpha_prev_f, beta_prev_f, size=n_sample)

# incidence, screening and proportion of incident infections asymptomatic in men

inc_m = np.zeros(n_sample)
scr_m = np.zeros(n_sample)
p_asymp_m = np.zeros(n_sample)

for i in xrange(n_sample):
    def tmpfun(inc, scr, p_asymp):
        [tr, dr] = test_diag_fun(
            array([
                inc,
                scr,
                1-p_asymp, # proportion of incident infections which are symptomatic
                sc_m[i], # rate of self-clear 
                att_symp[i],
                p_true_pos_m[i], 
                p_false_pos_m[i]
                ]))
        prev = dyn_fun(
            inc*p_asymp, 
            sc_m[i] + scr*p_true_pos_m[i], 
            inc*(1-p_asymp), 
            scr*p_true_pos_m[i] + att_symp[i]*p_true_pos_m[i]
        ) 
        return (tr - test_rate_m_15_24[i], 
                dr - diag_rate_m_15_24[i], 
                prev - prev_m[i])

    [inc_m[i], scr_m[i], p_asymp_m[i]] = fsolve(lambda x: tmpfun(x[0], x[1], x[2]), [0.09, 0.25, 0.9] )

# incidence, screening and proportion of incident infections asymptomatic in women

inc_f = np.zeros(n_sample)
scr_f = np.zeros(n_sample)
p_asymp_f = np.zeros(n_sample)

for i in xrange(n_sample):
    def tmpfun(inc, scr, p_asymp):
        [tr, dr] = test_diag_fun(
            array([
                inc,
                scr,
                1-p_asymp, # proportion of incident infections which are symptomatic
                sc_f[i], # rate of self-clear 
                att_symp[i],
                p_true_pos_f[i], 
                p_false_pos_f[i]
                ]))
        prev = dyn_fun(
            inc*p_asymp, 
            sc_f[i] + scr*p_true_pos_f[i], 
            inc*(1-p_asymp), 
            scr*p_true_pos_f[i] + att_symp[i]*p_true_pos_f[i]
        ) 
        return (tr - test_rate_f_15_24[i], 
                dr - diag_rate_f_15_24[i], 
                prev - prev_f[i])

    [inc_f[i], scr_f[i], p_asymp_f[i]] = fsolve(lambda x: tmpfun(x[0], x[1], x[2]), [0.09, 0.25, 0.9] )





