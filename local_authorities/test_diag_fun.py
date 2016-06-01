# this is the preliminary work so that you have a function that gives testing and diagnosis rates per person per unit time

from numpy import *
from sympy import *

##############################
# function to calculate steady-state A, U and S
##############################

A, U, S = symbols("A U S")
alpha_UA, alpha_AU, alpha_US, alpha_SU  = symbols("alpha_UA alpha_AU alpha_US alpha_SU")

model_dyn = [
    alpha_UA*U - alpha_AU*A,
    alpha_AU*A + alpha_SU*S - (alpha_UA + alpha_US)*U,
    alpha_US*U - alpha_SU*S,
    A + U + S - 1
    ]

sol_dyn = solve(model_dyn, A, U, S)

U_fun = lambdify((alpha_UA, alpha_AU, alpha_US, alpha_SU), sol_dyn[U])
A_fun = lambdify((alpha_UA, alpha_AU, alpha_US, alpha_SU), sol_dyn[A])
S_fun = lambdify((alpha_UA, alpha_AU, alpha_US, alpha_SU), sol_dyn[S])
dyn_fun = lambdify((alpha_UA, alpha_AU, alpha_US, alpha_SU), sol_dyn[A] + sol_dyn[S])

##############################
# now the model for observed testing and diagnosis rates
##############################

tsym, dsym, ssym, test_sym, true_pos, false_pos = symbols("tsym dsym ssym test_sym true_pos false_pos")

model_test_diag = [
    tsym - ( ssym + (1 - A - U)*test_sym ),
    dsym - ( A*ssym*true_pos + U*ssym*false_pos + (1 - A - U)*test_sym*true_pos )
    ]

sol_prev = solve(model_test_diag, ssym, U)
sol_test_diag = solve(model_test_diag, tsym, dsym)
test_fun = lambdify((A, U, ssym, test_sym, true_pos, false_pos), sol_test_diag[tsym])
diag_fun = lambdify((A, U, ssym, test_sym, true_pos, false_pos), sol_test_diag[dsym])

def test_diag_fun(parms):
    # parms = (incidence, screening rate, proportion symptomatic, self-cure rate, testing rate in symptomatics, true positive rate, false positive rate, ...)
    inc = parms[0]
    scr = parms[1]
    p_symp = parms[2]
    self_cure = parms[3]
    test_sym = parms[4]
    true_pos = parms[5]
    false_pos = parms[6]
    
    A = A_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, scr*true_pos + test_sym*true_pos)
    U = U_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, scr*true_pos + test_sym*true_pos)
    return [test_fun(A, U, scr, test_sym, true_pos, false_pos), diag_fun(A, U, scr, test_sym, true_pos, false_pos)]
    
##############################
# ...and if symptomatic and asymptomatic diagnoses are observed separately:
##############################

dssym, dasym = symbols("dssym dasym")

model_test_diag_sa = [
    tsym - ( (A + U)*ssym + (1 - A - U)*test_sym ),
    dssym - (1 - A - U)*test_sym*true_pos,
    dasym - ( A*ssym*true_pos + U*ssym*false_pos )
    ]

sol_test_diag_sa = solve(model_test_diag_sa, tsym, dssym, dasym)
test_sa_fun = lambdify((A, U, ssym, test_sym), sol_test_diag_sa[tsym])
diags_fun = lambdify((A, U, ssym, test_sym, true_pos), sol_test_diag_sa[dssym])
diaga_fun = lambdify((A, U, ssym, true_pos, false_pos), sol_test_diag_sa[dasym])

def test_diag_sym_asym_fun(parms):
    # parms = (incidence, screening rate, self-cure rate, testing rate in symptomatics, true positive rate, false positive rate, ...)
    inc = parms[0]
    scr = parms[1]
    p_symp = parms[2]
    self_cure = parms[3]
    test_sym = parms[4]
    true_pos = parms[5]
    false_pos = parms[6]
    
    A = A_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, self_cure + scr*true_pos + test_sym*true_pos)
    U = U_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, self_cure + scr*true_pos + test_sym*true_pos)
    return [test_sa_fun(A, U, scr, test_sym), diags_fun(A, U, scr, test_sym, true_pos), diaga_fun(A, U, scr, true_pos, false_pos)]
    
##############################
# use information on symptoms in _prevalent_ (as opposed to incident) infections
##############################

def test_diag_prev_symp_fun(parms):
    # parms = (incidence, screening rate, proportion symptomatic, self-cure rate, testing rate in symptomatics, true positive rate, false positive rate, ...)
    inc = parms[0]
    scr = parms[1]
    p_symp = parms[2]
    self_cure = parms[3]
    test_sym = parms[4]
    true_pos = parms[5]
    false_pos = parms[6]
    
    A = A_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, self_cure + scr*true_pos + test_sym*true_pos)
    U = U_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, self_cure + scr*true_pos + test_sym*true_pos)
    return [test_fun(A, U, scr, test_sym, true_pos, false_pos), diag_fun(A, U, scr, test_sym, true_pos, false_pos), (1 - A - U)/(1 - U)]
