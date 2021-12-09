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


def get_ks_p(obs_ranks):
    """
    Calculate ks-p value for observed ranks
    """
    stat,p = stats.kstest(
        rvs = obs_ranks,
        cdf = 'uniform',
        alternative='greater',
        mode='exact',
    )
    return p


########################
#   Metrics guidelines #
########################
#Each metric needs to take input in the form of an hdf5_cell
#The metric function calculates the per-gene score for all genes in the cell
#   For the periphery ranking, this will be based on the minimum distance of each spot to the periph
#   For the centrality ranking, this will be based on the minimum distance of each spot to the centroid
#   For the radial ranking this will be based on the min angle dist of each spot to the ref line

#The dataframe object where each row is a unique gene. The df must have the following columns:
#     'gene' The gene of the row
#     'annotation' The cell-type annotation information. same for all rows
#     'subclass' The cell-type annotation information. same for all rows
#     'num_genes' The number of unique genes in the cell. same for all rows
#     'num_gene_spots' The number of spots of this spots gene in the cell
#     'num_spots' The number of spots in the whole cell. same for all rows
#     'volume' optional, the volume of the cell, same for all rows

def calculate_ks_periphery_scores(hdf5_cell):
    """
    Returns a list/array of min distances of each spot to the cell boundary
    Operates on an individual cell
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

    #Rank and normalize the spots
    min_spot_genes = np.array(min_spot_genes)
    spot_ranks = np.array(min_periph_dists).argsort().argsort()
    norm_spot_ranks = spot_ranks/len(spot_ranks)

    #Iterate through unique genes to get per-gene score
    genes = np.unique(min_spot_genes)
    num_genes_per_spot = []
    gene_scores = []
    for gene in genes:
        gene_inds = min_spot_genes == gene
        num_genes_per_spot.append(sum(gene_inds))
        
        gene_ranks = norm_spot_ranks[gene_inds]
        ks_p = get_ks_p(gene_ranks)
        gene_scores.append(ks_p)

    
    cell_info = pd.DataFrame({
        'gene':[g.decode() for g in genes],
        'gene_score':gene_scores,
        'annotation':hdf5_cell.attrs.get('annotation'),
        'subclass':hdf5_cell.attrs.get('subclass'),
        'num_genes':hdf5_cell.attrs.get('num_genes'),
        'num_gene_spots':num_genes_per_spot,
        'num_spots':hdf5_cell.attrs.get('num_spots'),
        'volume':hdf5_cell.attrs.get('volume'),
    })

    return cell_info


def calculate_ks_centrality_scores(hdf5_cell):
    """
    Returns a list/array of min distances of each spot to the cell centroid 
    Operates on an individual cell

    Calculates (x-x0)^2+(y-y0)^2 where (x0,y0) is the cell centroid
    Don't have to take the sqrt for the L2 norm because I'm just comparing to spots
    """
    min_centroid_dists = []
    min_spot_genes = []
    for zslice in hdf5_cell.attrs.get('zslices'):

        #Calculate dists of each spot to centroid
        boundary = hdf5_cell['boundaries'][zslice]
        spot_coords = hdf5_cell['spot_coords'][zslice]
        spot_genes = hdf5_cell['spot_genes'][zslice]

        centroid = np.mean(boundary,axis=0)
        spot_norms = spot_coords-centroid
        centroid_dists = np.sum(spot_norms*spot_norms,axis=1)
        
        min_centroid_dists.extend(centroid_dists)
        min_spot_genes.extend(spot_genes)

        
    min_spot_genes = np.array(min_spot_genes)
    spot_ranks = np.array(min_centroid_dists).argsort().argsort()
    norm_spot_ranks = spot_ranks/len(spot_ranks)


    #Iterate through unique genes to get per-gene score
    genes = np.unique(min_spot_genes)
    num_genes_per_spot = []
    gene_scores = []
    for gene in genes:
        gene_inds = min_spot_genes == gene
        num_genes_per_spot.append(sum(gene_inds))

        gene_ranks = norm_spot_ranks[gene_inds]
        ks_p = get_ks_p(gene_ranks)
        gene_scores.append(ks_p)

    cell_info = pd.DataFrame({
        'gene':[g.decode() for g in genes],
        'gene_score':gene_scores,
        'annotation':hdf5_cell.attrs.get('annotation'),
        'subclass':hdf5_cell.attrs.get('subclass'),
        'num_genes':hdf5_cell.attrs.get('num_genes'),
        'num_gene_spots':num_genes_per_spot,
        'num_spots':hdf5_cell.attrs.get('num_spots'),
        'volume':hdf5_cell.attrs.get('volume'),
    })

    return cell_info


def calculate_ks_radial_scores(hdf5_cell):
    """
    Returns a list/array of min angle distances of each spot to the reference line
    Operates on an individual cell
    """

    #Step 0: Calculate the absolute value angle of each spot to the line segment (0,0) --> (1,0) for each zslice
    horiz_angs = []
    spot_genes = []
    for zslice in hdf5_cell.attrs.get('zslices'):

        #Calculate dists of each spot to centroid
        z_boundary = hdf5_cell['boundaries'][zslice]
        z_spot_coords = hdf5_cell['spot_coords'][zslice]
        z_spot_genes = hdf5_cell['spot_genes'][zslice]

        z_centroid = np.mean(z_boundary,axis=0)
        centered_spots = z_spot_coords-z_centroid
        x = centered_spots[:,0]
        y = centered_spots[:,1]
        z_horiz_angs = np.abs(np.arctan2(y,x))
       
        horiz_angs.extend(z_horiz_angs)
        spot_genes.extend(z_spot_genes)

    horiz_angs = np.array(horiz_angs)
    spot_genes = np.array(spot_genes)
    
    #Iterate through each GENE doing the following steps:
    #Step 1: Calculate median angle of all the GENE spots only
    #Step 2: Calculate abs-value of angle of ALL spots to this median
    #Step 3: Rank and normalize ALL spots based on their med_ang_dist
    #Step 4: Retrieve the normalized ranks JUST for spots of GENE
    #Step 5: Run KS-test with just the normalized GENE ranks
    genes = np.unique(spot_genes)
    num_genes_per_spot = []
    gene_scores = []
    for gene in genes:
        gene_inds = spot_genes == gene
        num_genes_per_spot.append(sum(gene_inds))

        med_gene_ang = np.median(horiz_angs[gene_inds])
        med_ang_dists = np.abs(horiz_angs-med_gene_ang)
        spot_ranks = med_ang_dists.argsort().argsort()
        norm_spot_ranks = spot_ranks/len(spot_ranks)
 
        norm_gene_spot_ranks = norm_spot_ranks[gene_inds]
        
        ks_p = get_ks_p(norm_gene_spot_ranks)
        gene_scores.append(ks_p)

    
    cell_info = pd.DataFrame({
        'gene':[g.decode() for g in genes],
        'gene_score':gene_scores,
        'annotation':hdf5_cell.attrs.get('annotation'),
        'subclass':hdf5_cell.attrs.get('subclass'),
        'num_genes':hdf5_cell.attrs.get('num_genes'),
        'num_gene_spots':num_genes_per_spot,
        'num_spots':hdf5_cell.attrs.get('num_spots'),
        'volume':hdf5_cell.attrs.get('volume'),
    })

    return cell_info


def calculate_periphery_dists(hdf5_cell):
    """
    Helper function to return min dists to the periphery

    Not used in nf pipelin run
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

    df = pd.DataFrame({
        'gene':[g.decode() for g in min_spot_genes],
        'dist':min_periph_dists,
        'rank':spot_ranks,
        'annotation':hdf5_cell.attrs.get('annotation'),
    })

    return df


def calculate_rs_periphery_scores(hdf5_cell):
    """
    Returns a list/array of min distances of each spot to the cell boundary
    Operates on an individual cell

    Returns a dataframe where each row is a unique gene in this cell
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
    obs_medians = []
    exp_medians = []
    for gene in genes:
        gene_inds = min_spot_genes == gene
        num_spots = int(sum(gene_inds))
        gene_ranks = spot_ranks[gene_inds]

        """
        #hack, randomly dropping a spot if even number of points
        #doing this because calculating probability is two orders of magnitude
        #faster for odd-m than even-m
        if num_spots%2==0:
            num_spots -= 1
            gene_ranks = np.random.choice(gene_ranks,num_spots,replace=False)
        """

        num_genes_per_spot.append(num_spots)
        obs_median = np.median(gene_ranks)
        rs_p = p_le_med(num_spots,tot_spots,obs_median)

        gene_scores.append(rs_p)
        obs_medians.append(obs_median)
        exp_medians.append(tot_spots/2) #don't even really need to calculate this here, could calc on output df

    
    cell_info = pd.DataFrame({
        'gene':[g.decode() for g in genes],
        'gene_score':gene_scores,
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


def calculate_rs_centrality_scores(hdf5_cell):
    pass


def calculate_rs_radial_scores(hdf5_cell):
    pass


def calculate_linear_median_effect_sizes(hdf5_cell):
    """
    Operates on an individual cell

    Returns a dataframe where each row is a unique gene in this cell
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
        var = utils.calc_var_effect_size(num_spots,tot_spots)

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


def _test(hdf5_path, cell_id):
    """
    Test metric
    """
    f = h5py.File(hdf5_path,'r')
    cell = f['cells'][cell_id]
    name = os.path.basename(cell.name)
    f.close()

    return name

def _slow_test(hdf5_path, cell_id):
    """
    Test metric
    """
    f = h5py.File(hdf5_path,'r')
    cell = f['cells'][cell_id]
    name = os.path.basename(cell.name)
    f.close()
    time.sleep(1)
    return name


########################
#   Rank calculator    #
########################
def calculate_ranks(hdf5_path, metric, cell_ids=[]):
    """
    Iterate through all cells in the hdf5 file, calculating the per-gene scores
    the function will calculate each cell individually and yield a dataframe object for each cell

    The columns and rows of the df are described in the comments above the metric functions
    """
    with h5py.File(hdf5_path,'r') as f:

        #If subset of cell_ids is not passed in, then process them all
        if len(cell_ids) == 0:
            cell_ids = f['cell_ids']

        for i,cell_id in enumerate(cell_ids):

            #Apply the correct metric to this cell
            hdf5_cell = f['cells'][cell_id.decode()]
            metric_func = metric_lookup[metric]
            gene_scores = metric_func(hdf5_cell)

            gene_scores['cell_id'] = cell_id
            gene_scores['metric'] = metric

            yield gene_scores


#Metric lookup to map string to metric function
metric_lookup = {
    'ks_periph': calculate_ks_periphery_scores,
    'ks_central': calculate_ks_centrality_scores,
    'ks_radial': calculate_ks_radial_scores,
    'rs_periph': calculate_rs_periphery_scores,
    'rs_central': calculate_rs_centrality_scores,
    'rs_radial': calculate_rs_radial_scores,
    'lin_eff_periph':calculate_linear_median_effect_sizes,
}


