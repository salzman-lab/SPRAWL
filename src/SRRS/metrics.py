import shapely.geometry
from scipy import stats

import pandas as pd
import numpy as np
import collections
import h5py
import os

import operator as op
from functools import reduce
import math
import time

from . import utils

def ncr(n, r):
    """
    Adapted from https://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python/4941846
    """
    if r > n:
        return 0

    if r == 0:
        return 1
    
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom

def p_ord_stat(m,n,k,x):
    """
    Calculates the probability that the k-th ranked item in the sample
    has order statistic x when choosing m spots out of n without replacement
    
    Taken from http://www.math.wm.edu/~leemis/2005informsjoc.pdf
    """
    numer = ncr(x-1,k-1)*ncr(n-x,m-k)
    denom = ncr(n,m)
    if numer == 0:
        return 0
    
    return numer/denom
    
def p_le_med(m,n,obs_med):
    """
    Calculates the probability of observing the median value x
    or a smaller value when choosing m spots out of n
    
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
    #If obs_median is smaller than half, calculate area under left-side of PDF
    #otherwise calculate area under right-side and subtract that from 1
    count_left = obs_med <= n/2
   
    if m%2 == 1:
        #Odd case
        k = (m+1)//2
        min_med = k
        max_med = n-k+1
        
        if count_left:
            xs = np.arange(min_med,obs_med+1,dtype=int)
        else:
            xs = np.arange(obs_med+1,max_med+1,dtype=int)
            
        denom = ncr(n,m)
        p = sum(ncr(x-1,k-1)*ncr(n-x,m-k)/denom for x in xs)
        
    else:
        #Even case
        l = m//2
        r = l+1
        
        min_med = (m+1)/2
        max_med = n-(m+1)/2+1
        
        if count_left:
            xs = np.arange(min_med,obs_med+0.5,0.5)
        else:
            xs = np.arange(obs_med+0.5,max_med+0.5,0.5)
            
        p = 0
        Rl_cache = {}
        
        denom = ncr(n,m)
        
        for x in xs:
        
            Rl = math.ceil(x-1)
            Rr = math.floor(x+1)
            
            while Rl >= m//2 and Rr <= n-m//2+1:
                #simplify P(Rr|Rl) since k=1
                
                #p_ord_stat(m,n,k,x):
                #new_n = n-Rl
                #new_m = m-l
                #new_k = 1
                #new_x = Rr-Rl
                
                #p = ncr(x-1,k-1)*ncr(n-x,m-k)/ncr(n,m)
                p_r_given_l = ncr(n-Rr,m-l-1)/ncr(n-Rl,m-l)
                
                if Rl not in Rl_cache:
                    Rl_cache[Rl] = ncr(Rl-1,l-1)*ncr(n-Rl,m-l)/denom
                    
                p += p_r_given_l*Rl_cache[Rl]

                Rl -= 1
                Rr += 1
                       
    #If we were counting the right-side area, then return 1-p
    if count_left:
        return p
    else:
        return 1-p

def calc_exp_midord(m,n):
    """
    Calculates the expected median/midord
    """
    i = (m+1)//2 #median when m is odd, one "left" of the median when m is even
    numerator = sum(x*ncr(x-1,i-1)*ncr(n-x,m-i) for x in range(1,n+1))
    denominator = ncr(n,m)
    return numerator/denominator


########################
#   Metrics guidelines #
########################
#Each metric needs to take input in the form of an hdf5_path and a cell_id
#The metric function calculates the per-gene score for all genes in the cell
#   For the periphery ranking, this will be based on the minimum distance of each spot to the periph
#   For the radial ranking this will be based on the min angle dist of each spot to the gene radial center

#Return value will be a pandas dataframe where each row is a unique gene. The df is expected to have the following columns:
#     'cell_id' the cell id. will be the same for every row
#     'annotation' The cell-type annotation information. same for all rows
#     'num_spots' The number of spots in the whole cell. same for all rows
#     'gene' The gene of the row
#     'num_gene_spots' The number of spots of this spots gene in the cell
#     'metric' peripheral/polar etc
#     'score' value ranging from -1 to 1
#     'variance' theoretical variance under the null

def _test(hdf5_path, cell_id):
    """
    Test metric
    Returns a test df with all columns as expected except
        metric - test
        score - 0
    """

    with h5py.File(hdf5_path,'r') as f:
        cell = f['cells'][cell_id]
        gene_counts = collections.defaultdict(int)

        for z_ind in cell.attrs['zslices']:
            z_counts = collections.Counter(cell['spot_genes'][z_ind])
            gene_counts.update(z_counts)

        #drop a spot for genes with an even number of spots
        n = sum(gene_counts.values())
        ms = [m if m%2 == 1 else m-1 for m in gene_counts.values()]
        vs = [utils.calc_var(m,n) for m in ms]

        df = pd.DataFrame({
            'cell_id':cell_id,
            'annotation':cell.attrs['annotation'],
            'num_spots':n,
            'gene':gene_counts.keys(),
            'num_gene_spots':gene_counts.values(),
            'metric':'test',
            'score':0,
            'variance':vs,
        })

    return df


def peripheral(hdf5_path, cell_id):
    """
    Peripheral metric
    """
    min_periph_dists = []
    min_spot_genes = []

    for zslice in hdf5_cell.attrs.get('zslices'):

        #Calculate dists of each spot to periphery
        boundary = hdf5_cell['boundaries'][zslice]
        spot_coords = hdf5_cell['spot_coords'][zslice]
        spot_genes = hdf5_cell['spot_genes'][zslice]

        poly = shapely.geometry.Polygon(boundary)
        for p,gene in zip(spot_coords,spot_genes):
            dist = poly.boundary.distance(shapely.geometry.Point(p))
            min_periph_dists.append(dist)
            min_spot_genes.append(gene)

    #Rank the spots
    min_spot_genes = np.array(min_spot_genes)
    spot_ranks = np.array(min_periph_dists).argsort().argsort()+1 #add one so ranks start at 1 rather than 0
    tot_spots = len(spot_ranks)

    #Iterate through unique genes to get per-gene score
    genes = np.unique(min_spot_genes)
    num_genes_per_spot = []
    gene_scores = []
    mn_vars = []
    obs_medians = []
    exp_medians = []
    for gene in genes:
        gene_inds = min_spot_genes == gene
        num_spots = int(sum(gene_inds))

        gene_ranks = spot_ranks[gene_inds]

        #Drop one spot randomly (uniformly) if an even number
        if num_spots%2 == 0:
            num_spots -= 1
            gene_ranks = np.random.choice(gene_ranks,num_spots,replace=False)

        num_genes_per_spot.append(num_spots)
        obs_med = np.median(gene_ranks)
        eff = utils.f_effect(obs_med, tot_spots)
        var = utils.calc_var(num_spots,tot_spots)

        gene_scores.append(eff)
        mn_vars.append(var)
        obs_medians.append(obs_med)
        exp_medians.append((tot_spots+1)/2)

    
    cell_info = pd.DataFrame({
        'gene':[g.decode() for g in genes],
        'gene_score':gene_scores,
        'mn_var':mn_vars,
        'obs_midord':obs_medians,
        'exp_midord':exp_medians,
        'annotation':hdf5_cell.attrs.get('annotation'),
        'subclass':hdf5_cell.attrs.get('subclass'),
        'num_genes':hdf5_cell.attrs.get('num_genes'),
        'num_gene_spots':num_genes_per_spot,
        'num_spots':hdf5_cell.attrs.get('num_spots'),
        'volume':hdf5_cell.attrs.get('volume'),
    })

    return cell_info



