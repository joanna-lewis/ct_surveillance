#########################
# define a function to calculate multinomial pmf
# see http://stackoverflow.com/questions/13903922/multinomial-pmf-in-python-scipy-numpy
#########################

import math
from scipy.stats import beta
from numpy.random import normal

class Multinomial(object):
  def __init__(self, params):
    self._params = params

  def pmf(self, counts):
    if not(len(counts)==len(self._params)):
      raise ValueError("Dimensionality of count vector is incorrect")

    prob = 1.
    for i,c in enumerate(counts):
#      print self._params[i]
      prob *= self._params[i]**counts[i]
      print prob

    if prob != 0:
      return prob * math.exp(self._log_multinomial_coeff(counts))
    else:
      return 0

  def log_pmf(self,counts):
    if not(len(counts)==len(self._params)):
      raise ValueError("Dimensionality of count vector is incorrect")

    prob = 0.
    for i,c in enumerate(counts):
      prob += counts[i]*math.log(self._params[i])

    return prob + self._log_multinomial_coeff(counts)

  def _log_multinomial_coeff(self, counts):
    return self._log_factorial(sum(counts)) - sum(self._log_factorial(c)
                                                    for c in counts)

  def _log_factorial(self, num):
    if not round(num)==num and num > 0:
      raise ValueError("Can only compute the factorial of positive ints")
    return sum(math.log(n) for n in range(1,num+1))
