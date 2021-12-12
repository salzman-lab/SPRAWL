import pytest
from SRRS import utils
import numpy as np

@pytest.mark.parametrize(
    'm,n',
    [
        (1,1),
        (10,10),
        (100,100),
        (1000,1000),
        (3000,3000),
    ]
)
def test_m_equals_n_var_is_zero(m,n):
    assert utils.calc_var(m,n) == 0


@pytest.mark.parametrize('n',[14,49,892])
def test_given_n_increasing_m_lowers_vars(n):
    #just choose 15 values of m, even for large n
    ms = sorted(np.random.choice(range(1,n,2), 15))

    vs = [utils.calc_var(m,n) for m in ms]

    #reverse sort since ms increase, but vs should decrease
    assert vs == sorted(vs, reverse=True)


@pytest.mark.parametrize('seed',[17,23])
@pytest.mark.parametrize(
    'm,n',
    [
        (1,10),
        (3,157),
        (19,157),
        (155,157),
        (4,17),
        (16,18),
        (14,178),
    ]
)
def test_empirical_var_under_null(seed,m,n):
    np.random.seed(seed) #<-- just so the test always fails or always succeeds for a given code
    theory_var = utils.calc_var(m,n)

    threshold = 5e-3
    its = 10000
    null_meds = [np.median(np.random.choice(n,m,replace=False)) for _ in range(its)]
    null_scores = [utils.score(med,n) for med in null_meds]
    emp_var = np.var(null_scores)

    assert abs(theory_var-emp_var) < threshold


