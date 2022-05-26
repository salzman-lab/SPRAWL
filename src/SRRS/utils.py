import multiprocessing as mp
import operator as op
import pandas as pd
import numpy as np
import collections
import functools
import tempfile
import random
import scipy
import pysam
import math
import os


p_med_cache = {}
ncr_mem = {}

def score(obs_med, n):
    """
    Calculate SRRS score from median observed rank
    """
    if obs_med == None:
        return None

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
        numer = functools.reduce(op.mul, range(n, n-r, -1), 1)
        denom = functools.reduce(op.mul, range(1, r+1), 1)
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
    if (m,n,obs_med) in p_med_cache:
        return p_med_cache[(m,n,obs_med)]

    if m%2 == 1:
        #Odd case
        if not obs_med == int(obs_med):
            #Odd m cant have fractional median
            p = 0
        else:
            obs_med = int(obs_med)
            k = (m+1)//2
            p = ncr(obs_med-1,k-1)*ncr(n-obs_med,m-k)/ncr(n,m)

    else:
        #Even case
        l = m//2
        r = l+1

        Rl = math.ceil(obs_med-1)
        Rr = math.floor(obs_med+1)

        p = 0

        #Iterate through possible Rl and Rr values for the given median
        while Rl >= m//2 and Rr <= n-m//2+1:
            #simplify P(Rr|Rl) since k=1

            #p_ord_stat(m,n,k,x):
            #new_n = n-Rl
            #new_m = m-l
            #new_k = 1
            #new_x = Rr-Rl

            #p = ncr(x-1,k-1)*ncr(n-x,m-k)/ncr(n,m)
            p_r_given_l = ncr(n-Rr,m-l-1)
            p += p_r_given_l*ncr(Rl-1,l-1)/ncr(n,m)

            Rl -= 1
            Rr += 1

    p_med_cache[(m,n,obs_med)] = p
    return p

def _calc_var_helper(m,n):
    """
    Helper function to calculate the variance
    """
    if m%2 == 1:
        #m odd case
        min_med = m//2+1
        max_med = n//2
        possible_meds = np.arange(min_med,max_med+1)

    else:
        #m even case
        min_med = (m+1)/2
        max_med = n//2
        possible_meds = np.arange(min_med,max_med+1,0.5)

    var = 2*sum(p_med(m,n,med)*score(med,n)**2 for med in possible_meds)

    return var


def p_two_sided_med(m,n,obs_med):
    """
    Calculates the two-sided p-value
    of observing such an extreme median score

    Assume
        there are 7 total spots (n = 7)
        there are 3 gene spots (m = 3)
        the obs median is 2 (obs_med = 2)

    There is only one way for this to occur
    the ranks must be [1,2,3] so the p_two_sided_med = p_med(3,7,2)

    There is symmetry, the following situation has the same probability
        7 total spots (n = 7)
        3 gene spots (m = 3)
        the obs median is 6 (obs_med = 6)

    Now
        7 total spots (n = 7)
        3 gene spots (m = 3)
        the obs median is 3 (obs_med = 3)

    then p_two_sided_med = p_med(3,7,2)+p_med(3,7,3)

    the symmetry is accounted for by flipping larger observed medians to be
    on the smaller side, the same distance from the expected median

    then the final p-value is multiplied by 2
    """
    exp_med = (n+1)/2
    p_eq_med = p_med(m,n,obs_med)

    if obs_med > exp_med:
        #symmetrical distribution around expected median, flip to be on < side
        obs_med = 2*exp_med-obs_med

    if m%2 == 1:
        #m odd case
        min_med = m//2+1
        possible_meds = np.arange(min_med,obs_med)

    else:
        #m even case
        min_med = (m+1)/2
        possible_meds = np.arange(min_med,obs_med,0.5)

    p_twosided = 2*(sum(p_med(m,n,med) for med in possible_meds)+p_eq_med)
    return p_twosided



def calc_var(m,n,approx_evens=False):
    """
    Calculates the variance in effect sizes for choosing m spots out of n
    Uses the formula below

    Var(X) = E[(X-u)**2]
    Var(X) = E[X**2] #since u is always 0

    E[X**2] = sum(P(X**2=x**2)x**2)

    symmetry under the null allows to calculate only half the possible meds

    Includes the option to approximate the variance when m is even
    """
    if m == n:
        #special case
        return 0

    if approx_evens and m%2 == 0:
        #m even case approximation as average of R and L odds
        r_var = _calc_var_helper(m-1,n)
        l_var = _calc_var_helper(m+1,n)
        return (r_var+l_var)/2

    else:
        return _calc_var_helper(m,n)


def _cell_var(cell, var_mem={}):
    """
    Helper function to calculate theoretical gene variance
    utilizes manager.dict() shared memory
    adds gene_vars as a member of the cell object
    """
    #If the gene_vars are already calculated, just return the cell
    if cell.gene_vars:
        return cell

    n = cell.n
    for g,m in cell.gene_counts.items():
        if (m,n) not in var_mem:
            var_mem[(m,n)] = calc_var(m,n)
        cell.gene_vars[g] = var_mem[(m,n)]

    return cell


def _iter_vars(cells, processes=2):
    """
    Calculate the theoretical variance of each gene in each cell
    utilizes multiprocessing
    """
    manager = mp.Manager()
    var_mem = manager.dict()

    with mp.Pool(processes=processes) as p:
        f = functools.partial(_cell_var, var_mem = var_mem)
        for cell in p.imap_unordered(f, cells):
            yield cell

#Taken from https://stackoverflow.com/questions/22229796/choose-at-random-from-combinations
def random_combination(iterable, r):
    """
    Random selection from itertools.combinations(iterable, r)
    """
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n), r))
    return tuple(pool[i] for i in indices)


def random_mean_pairs_dist(spots, num_pairs):
    """
    Helper function to choose 'num_pairs' pairs of gene spots and calculate the mean distance for each gene
    Input is an array of spots that it will choose from
    """
    d = 0
    for _ in range(num_pairs):
        (x1,y1),(x2,y2) = random_combination(spots,2)
        d += math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))

    return d/num_pairs

def random_mean_pairs_angle(spots, centroid, num_pairs):
    """
    Helper function to choose 'num_pairs' pairs of gene spots and calculate the mean angle for each gene
    Input is an array of spots that it will choose from and the cell centroid as a tuple of (x,y)
    """
    cx,cy = centroid
    ang_sum = 0
    for _ in range(num_pairs):
        (x1,y1),(x2,y2) = random_combination(spots,2)
        v1 = (x1-cx,y1-cy)
        v2 = (x2-cx,y2-cy)
        ang_sum += angle_between(v1, v2)

    return ang_sum/num_pairs


#taken directly from https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def _map_bam_tag_chr(chrom, bam_path, mapping, key_tag, val_tag, tmp_dir):
    """
    Helper function to map bam tags
    creates a bam output file in a tmp_dir

    returns the path to the tmp bam file
    """
    chr_out_path = os.path.join(tmp_dir.name, '{}.bam'.format(chrom))

    in_bam = pysam.AlignmentFile(bam_path)
    out_bam = pysam.AlignmentFile(chr_out_path, 'wb', template=in_bam)

    for i,read in enumerate(in_bam.fetch(chrom)):

        try:
            k = read.get_tag(key_tag)
        except KeyError:
            continue

        if k not in mapping:
            continue

        read.set_tag(val_tag, mapping[k])
        out_bam.write(read)

    out_bam.close()
    in_bam.close()

    return chr_out_path


def map_bam_tag(bam_path, out_path, mapping, key_tag='CB', val_tag='XO', processes=1):
    """
    Create a new sorted and indexed bam file with an additional bam tag field
    Potentially useful for adding celltype tag to a bam for downstream processing
    Supports multiprocessing by operating on individual chromosomes and then merging

    returns the out_path

    Arguments:
        bam_path [required]
            path to a position-sorted and indexed bam

        out_path [required]
            path such as my_outs/out.bam of where to store output
            a .bam.bai of the same prefix will also be created

        mapping [required]
            dictionary mapping keys from the key_tag to values in the val_tag
            as a concrete example could map cell-barcode (CB) to cell-type (XO) custom tag

        key_tag, val_tag
            the key and val tag names to use to perform the mapping

        processes
            the number of processes to use to multiplex the operation
    """
    tmp_dir = tempfile.TemporaryDirectory()

    #get the chroms for multiplexing
    with pysam.AlignmentFile(bam_path) as bam:
        chroms = bam.references

    #perform the mapping
    with mp.Pool(processes=processes) as p:
        f = functools.partial(
            _map_bam_tag_chr,
            bam_path = bam_path,
            mapping = mapping,
            key_tag = key_tag,
            val_tag = val_tag,
            tmp_dir = tmp_dir,
        )
        chrom_bam_paths = p.map(f, chroms)

    #merge the individual bam files and create index
    merge_args = ['-@',str(processes),str(out_path)] + chrom_bam_paths
    pysam.merge(*merge_args)

    pysam.index(str(out_path))


    #delete the tmpdir containing all the tmp bams
    tmp_dir.cleanup()

    return out_path


def readzs_proxy_score(bam_path, locus, stratify_tag=None, **kwargs):
    """
    Flexible function to create a ReadZs score proxy directly from a bam
    Can optionally stratify by a tag in the bam
    """
    min_ont_reads = kwargs.get('min_ont_reads',0)

    #handling locus
    try:
        chrom,start,end = locus
        start = int(start)
        end = int(end)
    except:
        raise Exception('Must pass in a tuple of locus info such as locus=("chr1",1,100)')


    #Counting reads per tag in this region
    def get_tag(read):
        try:
            return read.get_tag(stratify_tag)
        except KeyError:
            return None

    strat_func = get_tag if stratify_tag else (lambda read: 'All')

    count_data = {
        'strat':[],
        'pos':[],
    }

    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch(chrom,start,end):
            if r.pos < start or r.pos > end:
                continue

            strat = strat_func(r)
            if strat:
                count_data['strat'].append(strat)
                count_data['pos'].append(r.pos)

    count_df = pd.DataFrame(count_data)
    count_df = count_df.groupby('strat').filter(lambda g: len(g) > min_ont_reads)

    exp_med = (end+start)/2
    span = end-start
    ont_to_score = count_df.groupby('strat')['pos'].median().subtract(exp_med).div(span).to_dict()
    return ont_to_score

