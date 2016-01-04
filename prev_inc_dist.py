####################################################################
# This script is to get distributions for incidence, prevalence and screening rate
# based on data from England
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

# read data from the MSTIC study (Aicken et al. BMC Health Serv. Res.; 2015)
import csv

# men
with open('../../../data/MSTIC/mstic_timedata_m.csv', 'rU') as csvfile:
    timereader = csv.reader(csvfile, delimiter=',')
    times = [int(row[0]) for row in timereader]
times_m = times[0:37] # exclude waiting times longer than 1 year
with open('../../../data/MSTIC/mstic_timedata_m.csv', 'rU') as csvfile:
    timereader = csv.reader(csvfile, delimiter=',')
    prop = [int(row[1]) for row in timereader]
prop_m = prop[0:37]
with open('../../../data/MSTIC/mstic_timedata_m.csv', 'rU') as csvfile:
    timereader = csv.reader(csvfile, delimiter=',')
    cumprop = [int(row[2]) for row in timereader]
cumprop_m = cumprop[0:37]

# women
with open('../../../data/MSTIC/mstic_timedata_f.csv', 'rU') as csvfile:
    timereader = csv.reader(csvfile, delimiter=',')
    times = [int(row[0]) for row in timereader]
times_f = times[0:43]
with open('../../../data/MSTIC/mstic_timedata_f.csv', 'rU') as csvfile:
    timereader = csv.reader(csvfile, delimiter=',')
    prop = [int(row[1]) for row in timereader]
prop_f = prop[0:43]
with open('../../../data/MSTIC/mstic_timedata_f.csv', 'rU') as csvfile:
    timereader = csv.reader(csvfile, delimiter=',')
    cumprop = [int(row[2]) for row in timereader]
cumprop_f = cumprop[0:43]

##############
# Metropolis-Hastings to get a sample for rate of treatment in men
##############

i = 0
att_symp_m = empty(n_sample) # testing rate per person per year
ll_m = empty(n_sample) # log-likelihood
old = 15 # starting sample value
new = 15 # starting sample value
# simulate probabilities corresponding to MSTIC data
simp_new = exp(-array(times_m)*new/365.25) - append(exp(-array(times_m[1:len(times_m)])*new/365.25),0)

acc=0.
while i < n_sample: # to do samples for p_test_symp
    
    new = random.normal(old, 2) # generate a sample from normal distribution
    
    if new < 0:
        att_symp_m[i] = old # reject
        ll_m[i] = -1e10
    else:
        simp_old = exp(-array(times_m)*old/365.25) - append(exp(-array(times_m[1:len(times_m)])*old/365.25), 0)
        simp_new = exp(-array(times_m)*new/365.25) - append(exp(-array(times_m[1:len(times_m)])*new/365.25), 0)

        if sum(simp_new > 0) != len(times_m):
            att_symp_m[i] = old # reject
            ll_m[i] = -1e10
        else:
            # simulate probabilities corresponding to the MSTIC data
            m_old = Multinomial(simp_old)
            m_new = Multinomial(simp_new)
            log_ratio = m_new.log_pmf(prop_m) - m_old.log_pmf(prop_m)
    
            if log(random.uniform(0,1)) <  log_ratio:
                att_symp_m[i] = new # accept
                ll_m[i] = m_new.log_pmf(prop_m)
                old = new
                acc = acc+1
            else:
                att_symp_m[i] = old # reject
                ll_m[i] = m_old.log_pmf(prop_m)
        
    i = i+1
    
print acc/n_sample # print the proportion of samples accepted

##############
# Metropolis-Hastings to get a sample for rate of treatment in women
##############

i = 0
att_symp_f = empty(n_sample) # testing rate per person per year
ll_f = empty(n_sample) # log-likelihood
old = 15 # starting sample value
new = 15 # starting sample value
# simulate probabilities corresponding to MSTIC data
simp_new = exp(-array(times_f)*new/365.25) - append(exp(-array(times_f[1:len(times_f)])*new/365.25),0)

acc=0.
while i < n_sample: # to do samples for p_test_symp
    
    new = random.normal(old, 2) # generate a sample from normal distribution
    
    if new < 0:
        att_symp_f[i] = old # reject
        ll_f[i] = -1e10
    else:
        simp_old = exp(-array(times_f)*old/365.25) - append(exp(-array(times_f[1:len(times_f)])*old/365.25), 0)
        simp_new = exp(-array(times_f)*new/365.25) - append(exp(-array(times_f[1:len(times_f)])*new/365.25), 0)

        if sum(simp_new > 0) != len(times_f):
            att_symp_f[i] = old # reject
            ll_f[i] = -1e10
        else:
            # simulate probabilities corresponding to the MSTIC data
            m_old = Multinomial(simp_old)
            m_new = Multinomial(simp_new)
            log_ratio = m_new.log_pmf(prop_f) - m_old.log_pmf(prop_f)
    
            if log(random.uniform(0,1)) <  log_ratio:
                att_symp_f[i] = new # accept
                ll_f[i] = m_new.log_pmf(prop_f)
                old = new
                acc = acc+1
            else:
                att_symp_f[i] = old # reject
                ll_f[i] = m_old.log_pmf(prop_f)
        
    i = i+1
    
print acc/n_sample # print the proportion of samples accepted

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


###################################
# now use these samples to try and infer prevalence in men and women
###################################

# men first...
prev_m_15_24 = zeros(n_sample)
inc_m_15_24 = zeros(n_sample)
scr_m_15_24 = zeros(n_sample)

for i in xrange(n_sample):
    [inc_m_15_24[i], scr_m_15_24[i]] = fsolve(lambda x: test_diag_fun(concatenate([
                    x, array([
                            1-p_asymp_m[i], # proportion of incident infections which are symptomatic
                            sc_m[i], # rate of self-clear 
                            att_symp_m[i],
                            p_true_pos_m[i], 
                            p_false_pos_m[i]
                        ])])) - array([test_rate_m_15_24[i],diag_rate_m_15_24[i]]), [0.09, 0.25])
    prev_m_15_24[i] = dyn_fun(inc_m_15_24[i]*p_asymp_m[i], sc_m[i] + scr_m_15_24[i]*p_true_pos_m[i], inc_m_15_24[i]*(1-p_asymp_m[i]), sc_m[i] + scr_m_15_24[i]*p_true_pos_m[i] + att_symp_m[i]*p_true_pos_m[i])
    
# ... then women
prev_f_15_24 = zeros(n_sample)
inc_f_15_24 = zeros(n_sample)
scr_f_15_24 = zeros(n_sample)

for i in xrange(n_sample):
    [inc_f_15_24[i], scr_f_15_24[i]] = fsolve(lambda x: test_diag_fun(concatenate([
                    x, array([
                            1-p_asymp_f[i], # proportion of incident infections which are symptomatic
                            sc_f[i], # rate of self-clear 
                            att_symp_f[i],
                            p_true_pos_f[i], 
                            p_false_pos_f[i]
                        ])])) - array([test_rate_f_15_24[i],diag_rate_f_15_24[i]]), [0.03, 0.44])
    prev_f_15_24[i] = dyn_fun(inc_f_15_24[i]*p_asymp_f[i], sc_f[i] + scr_f_15_24[i]*p_true_pos_f[i], inc_f_15_24[i]*(1-p_asymp_f[i]), sc_f[i] + scr_f_15_24[i]*p_true_pos_f[i] + att_symp_f[i]*p_true_pos_f[i])