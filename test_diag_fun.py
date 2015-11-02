# this is the preliminary work so that you have a function that gives testing and diagnosis rates per person per unit time

from numpy import *
from sympy import *

##############################
# function to calculate steady-state A, S and I
##############################

A, S, I = symbols("A S I")
alpha_SA, alpha_AS, alpha_SI, alpha_IS  = symbols("alpha_SA alpha_AS alpha_SI alpha_IS")

model_dyn = [
    alpha_SA*S - alpha_AS*A,
    alpha_AS*A + alpha_IS*I - (alpha_SA + alpha_SI)*S,
    alpha_SI*S - alpha_IS*I,
    A + S + I - 1
    ]

sol_dyn = solve(model_dyn, A, S, I)

S_fun = lambdify((alpha_SA, alpha_AS, alpha_SI, alpha_IS), sol_dyn[S])
A_fun = lambdify((alpha_SA, alpha_AS, alpha_SI, alpha_IS), sol_dyn[A])
I_fun = lambdify((alpha_SA, alpha_AS, alpha_SI, alpha_IS), sol_dyn[I])
dyn_fun = lambdify((alpha_SA, alpha_AS, alpha_SI, alpha_IS), sol_dyn[A] + sol_dyn[I])

##############################
# now the model for observed testing and diagnosis rates
##############################

tsym, dsym, ssym, test_sym, true_pos, false_pos = symbols("tsym dsym ssym test_sym true_pos false_pos")

model_test_diag = [
    tsym - ( (A + S)*ssym + (1 - A - S)*test_sym ),
    dsym - ( A*ssym*true_pos + S*ssym*false_pos + (1 - A - S)*test_sym*true_pos )
    ]

sol_prev = solve(model_test_diag, ssym, S)
sol_test_diag = solve(model_test_diag, tsym, dsym)
test_fun = lambdify((A, S, ssym, test_sym, true_pos, false_pos), sol_test_diag[tsym])
diag_fun = lambdify((A, S, ssym, test_sym, true_pos, false_pos), sol_test_diag[dsym])

def test_diag_fun(parms):
    # parms = (incidence, screening rate, proportion symptomatic, self-cure rate, testing rate in symptomatics, true positive rate, false positive rate, ...)
    inc = parms[0]
    scr = parms[1]
    p_symp = parms[2]
    self_cure = parms[3]
    test_sym = parms[4]
    true_pos = parms[5]
    false_pos = parms[6]
    
    A = A_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, self_cure + test_sym*true_pos)
    S = S_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, self_cure + test_sym*true_pos)
    return [test_fun(A, S, scr, test_sym, true_pos, false_pos), diag_fun(A, S, scr, test_sym, true_pos, false_pos)]
    
##############################
# ...and if symptomatic and asymptomatic diagnoses are observed separately:
##############################

dssym, dasym = symbols("dssym dasym")

model_test_diag_sa = [
    tsym - ( (A + S)*ssym + (1 - A - S)*test_sym ),
    dssym - (1 - A - S)*test_sym*true_pos,
    dasym - ( A*ssym*true_pos + S*ssym*false_pos )
    ]

sol_test_diag_sa = solve(model_test_diag_sa, tsym, dssym, dasym)
test_sa_fun = lambdify((A, S, ssym, test_sym), sol_test_diag_sa[tsym])
diags_fun = lambdify((A, S, ssym, test_sym, true_pos), sol_test_diag_sa[dssym])
diaga_fun = lambdify((A, S, ssym, true_pos, false_pos), sol_test_diag_sa[dasym])

def test_diag_sym_asym_fun(parms):
    # parms = (incidence, screening rate, self-cure rate, testing rate in symptomatics, true positive rate, false positive rate, ...)
    inc = parms[0]
    scr = parms[1]
    p_symp = parms[2]
    self_cure = parms[3]
    test_sym = parms[4]
    true_pos = parms[5]
    false_pos = parms[6]
    
    A = A_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, self_cure + test_sym*true_pos)
    S = S_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, self_cure + test_sym*true_pos)
    return [test_sa_fun(A, S, scr, test_sym), diags_fun(A, S, scr, test_sym, true_pos), diaga_fun(A, S, scr, true_pos, false_pos)]