from functools import reduce
import operator as op
import numpy as np
import collections
import scipy
import math

ncr_mem = {}

def score(obs_med, n):
    """
    Calculate SRRS score from median observed rank
    """
    exp_med = (n+1)/2
    eff = (exp_med-obs_med)/(exp_med-1)
    return eff


def ncr(n, r):
    """
    Adapted from https://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python/4941846
    """
    if r > n:
        return 0

    if r == 0:
        return 1

    r = min(r, n-r)

    if (n,r) not in ncr_mem:
        numer = reduce(op.mul, range(n, n-r, -1), 1)
        denom = reduce(op.mul, range(1, r+1), 1)
        ncr_mem[(n,r)] = numer // denom

    return ncr_mem[(n,r)]


def p_med(m,n,obs_med):
    """
    Calculates the probability of observing the median value obs_med
    when choosing m spots out of n

    if m is odd, then the median is simply the (m+1)/2 order statistic. Call this Rc order statistic

    if m is even, then the median is the average of the m/2 and (m/2)+1 order statistics
    call these the Rl and Rr order statistics


    when m is odd the probability of observing median value x is the probability that order statistic
    Rc is equal to x

    when m is even the probability of observing median value x is the sum of all probabilities
    P(Rl+Rr = 2x) for all possible pairs of Rl and Rr

    Note that P(Rl and Rr) != P(Rl)*P(Rr) since they are not independent
    Instead,
        P(Rl and Rr) = P(Rr|Rl)P(Rl) #always true
                     = P_(m-l,n-Rl)(R1 = Rr-Rl)P(Rl)

    The intuition is that P(Rr|Rl) can be seen as a subproblem where we've already chosen Rl
    So now instead of 1 <= Rr <= n, we know Rl+1 <= Rr <= n, and 1 <= Rr-Rl <= n-Rl

    Now we just need to calculate the probability that the first rank equals Rr-Rl within the new bounds
    """
    if m%2 == 1:
        #Odd case
        obs_med = int(obs_med)
        k = (m+1)//2
        p = ncr(obs_med-1,k-1)*ncr(n-obs_med,m-k)

    else:
        #Even case
        l = m//2
        r = l+1

        p = 0

        Rl = math.ceil(obs_med-1)
        Rr = math.floor(obs_med+1)

        #Iterate through possible Rl and Rr values for the given median
        while Rl >= m//2 and Rr <= n-m//2+1:
            #simplify P(Rr|Rl) since k=1

            #p_ord_stat(m,n,k,x):
            #new_n = n-Rl
            #new_m = m-l
            #new_k = 1
            #new_x = Rr-Rl

            #p = ncr(x-1,k-1)*ncr(n-x,m-k)/ncr(n,m)
            p_r_given_l = ncr(n-Rr,m-l-1)/ncr(n-Rl,m-l)
            p += p_r_given_l*ncr(Rl-1,l-1)*ncr(n-Rl,m-l)

            Rl -= 1
            Rr += 1

    return p/ncr(n,m)


def calc_var(m,n):
    """
    Calculates the variance in effect sizes for choosing m spots out of n
    Uses the formula below

    Var(X) = E[(X-u)**2]
    Var(X) = E[X**2] #since u is always 0

    E[X**2] = sum(P(X**2=x**2)x**2)

    symmetry under the null allows to calculate only half the possible meds
    """
    if m%2 == 1:
        #m odd case
        min_med = m//2+1
        max_med = n//2
        possible_meds = np.arange(min_med,max_med+1)

    else:
        #even case
        min_med = (m+1)/2
        max_med = n//2
        possible_meds = np.arange(min_med,max_med+1,0.5)

    var = 2*sum(p_med(m,n,med)*score(med,n)**2 for med in possible_meds)

    return var



