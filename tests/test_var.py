import pytest
from SRRS import utils
from SRRS import scoring
from SRRS import cell
import numpy as np
import collections

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
        (2,1015),
        (48,1015),
        (6,557),
    ]
)
def test_empirical_exp_var_under_null(seed,m,n):
    np.random.seed(seed) #<-- just so the test always fails or always succeeds
    theory_var = utils.calc_var(m,n)

    threshold = 5e-3
    its = 30000

    #add 1 to medians to have ranks be 1-indexed
    null_meds = [np.median(np.random.choice(n,m,replace=False))+1 for _ in range(its)]
    null_scores = [utils.score(med,n) for med in null_meds]
    emp_var = np.var(null_scores)
    emp_exp = sum(x*f for x,f in collections.Counter(null_scores).items())/its

    assert abs(theory_var-emp_var) < threshold
    assert abs(emp_exp) < threshold


@pytest.mark.parametrize('seed',[23])
@pytest.mark.parametrize(
    'm,n,med',
    [
        (3,6,3.5),
        (1,10,2),
        (1,10,8),
        (1,10,2.5),
        (2,10,1),
        (2,10,7),
        (2,10,1.5),
        (2,10,5.5),
        (20,1527,57.5),
        (801,1527,98.5),
    ]
)
def test_empirical_p_med_under_null(seed,m,n,med):
    np.random.seed(seed) #<-- just so the test always fails or always succeeds

    threshold = 5e-3
    its = 30000

    null_meds = [np.median(np.random.choice(n,m,replace=False))+1 for _ in range(its)]

    theory_p_med = utils.p_med(m,n,med)
    emp_p_med = sum(m == med for m in null_meds)/its
    assert abs(theory_p_med-emp_p_med) < threshold

    theory_p_two_sided = utils.p_two_sided_med(m,n,med)
    emp_p_le = sum(m < med for m in null_meds)/its
    emp_p_ge = sum(m > med for m in null_meds)/its
    emp_p_two_sided = 2*(min(emp_p_le, emp_p_ge)+emp_p_med)
    assert abs(theory_p_two_sided-emp_p_two_sided) < threshold


@pytest.mark.parametrize(
    'm_n_meds',[
        ((4,50),(1,4,20.5,30.5)),
        ((3,50),(1,4,27,49)),
        ((5,150),(18,99)),
        ((20,1527),(57.5,)),
    ],
)
def test_p_twosided_helper(m_n_meds):
    ps = scoring._calc_p_twosided_helper(m_n_meds)

    np.random.seed(1)
    its = 30000
    threshold = 5e-3
    m,n = m_n_meds[0]
    null_meds = [np.median(np.random.choice(n,m,replace=False))+1 for _ in range(its)]

    for (m,n,med),p in ps.items():
        emp_p_med = sum(m == med for m in null_meds)/its
        emp_p_le = sum(m < med for m in null_meds)/its
        emp_p_ge = sum(m > med for m in null_meds)/its
        emp_p_two_sided = 2*(min(emp_p_le, emp_p_ge)+emp_p_med)

        assert abs(emp_p_two_sided-p) < threshold




@pytest.mark.slow
def test_iter_vars(m1s4):
    cells = m1s4.iter_cells()
    cells = scoring._iter_vars(cells)

    for c in cells:
        assert type(c) == cell.Cell
        assert c.gene_vars.keys() == c.gene_counts.keys()

        for i,(g,v) in enumerate(c.gene_vars.items()):
            m = c.gene_counts[g]
            n = c.n
            assert v == utils.calc_var(m,n)

            if i >= 1:
                break #stop after 2 genes as a spot check for each cell



