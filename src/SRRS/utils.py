import pandas as pd
import numpy as np
import itertools
import argparse
import h5py
import os

from scipy import stats
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

import operator as op
from functools import reduce
import math
import collections
import scipy


def calc_romano_tn(var_df):
    """
    Code from JuliaO to calculate romano statistic

    Notes from JuliaO (paraphrased by me):
        This function takes in a dataframe with one row per ontology (or whatever group you're analyzing)
        The dataframe must have columns
        * num_cells_ont (the number of cells in the ontology)
        * ont_median (the median of the scores for that ontology - this could be another summary statistic as well, like the mean)
        * ont_var (the sample variance in the scores).

    Returns the T_n value statistic
    """

    # calculate the inner sum that's subtracted
    num = 0
    denom = 0
    for index, row in var_df.iterrows():
      num += row["num_cells_ont"]*row["ont_median"]/row["ont_var"]
      denom += row["num_cells_ont"]/row["ont_var"]
    const = num/denom
    # calculate the outer sum
    tn = 0
    for index, row in var_df.iterrows():
      tn += (row["num_cells_ont"]/row["ont_var"])*(row["ont_median"] - const)**2
    return tn



def romano_anova_permute(df, val_col, test_by_col='gene', cat_col='annotation', num_perms=1000):
    """
    Calculate ANOVA p-values based on the Romano test statistic

    Inputs:
    - df is a dataframe which has one row per gene/cell combination with the following column names
    - val_col is the name of the column to use as values (periph scores for example)
    - test_by_col is the name of the column that will be used to differentiate separate tests (gene name for example)
    - cat_col is the ANOVA category column (cell annotation for example)

    Performs num_perms permutations per test and returns a dictionary keyed by test_by_col and valued by two sided p-value
    """
    p_perms = {}

    #Loop through each group (usually gene) calculating the obs tn and the perm tns
    for gene,g in df.groupby(test_by_col):
       
        #Group by category to get df to plug into calc_romano_tn
        var_df = g.groupby(cat_col).agg(
            num_cells_ont = (val_col,'size'),
            ont_median = (val_col,'median'),
            ont_var = (val_col,'var'),
        )
        
        obs_tn = calc_romano_tn(var_df)
        
        #Permute the category labels for the cells
        num_lt = 0
        for perm_num in range(num_perms):
            perm_var_df = pd.DataFrame({
                cat_col:g[cat_col].sample(g.shape[0], replace=False).values,
                val_col:g[val_col].values,
            }).groupby(cat_col).agg(
                num_cells_ont = (val_col,'size'),
                ont_median = (val_col,'median'),
                ont_var = (val_col,'var'),
            )

            perm_tn = calc_romano_tn(perm_var_df)
            num_lt += perm_tn <= obs_tn
            
        frac_lt = num_lt/num_perms
        p_perm = 2*min(frac_lt, 1-frac_lt)
        p_perms[gene] = p_perm

    return p_perms



def plot_cell_spots_zslices(hdf5_cell, spot_colors={}, default_spot_color='grey'):

    zslices = hdf5_cell.attrs['zslices']
    fig, axs = plt.subplots(figsize=(8,8),nrows=3,ncols=3,sharex=True,sharey=True)
    axs = axs.flatten()

    global_min_x,global_min_y = None,None
    global_max_x,global_max_y = None,None

    for i,zslice in enumerate(zslices):

        #Add the spots
        colors = []
        for gene in hdf5_cell['spot_genes'][zslice]:
            if gene.decode() in spot_colors:
                colors.append(spot_colors[gene.decode()])
            else:
                colors.append(default_spot_color)

        #Draw the cell outline
        boundary = hdf5_cell['boundaries'][zslice][:]
        min_x,min_y = boundary.min(axis=0)
        max_x,max_y = boundary.max(axis=0)

        if not global_min_x or min_x < global_min_x:
            global_min_x = min_x
        if not global_min_y or min_y < global_min_y:
            global_min_y = min_y
        if not global_max_x or max_x > global_max_x:
            global_max_x = max_x
        if not global_max_y or max_y > global_max_y:
            global_max_y = max_y

        polygon = Polygon(boundary, fill=None)
        axs[i].add_artist(polygon)


        axs[i].scatter(
            x = hdf5_cell['spot_coords'][zslice][:,0],
            y = hdf5_cell['spot_coords'][zslice][:,1],
            alpha = 0.8,
            color = colors,
        )


    for used_i in range(i+1):
        axs[used_i].set_xticks([])
        axs[used_i].set_yticks([])
        axs[used_i].set_xlim(global_min_x,global_max_x)
        axs[used_i].set_ylim(global_min_y,global_max_y)
        
        
    for unused_i in range(i+1,len(axs)):
        axs[unused_i].set_axis_off()

    plt.subplots_adjust(wspace=0, hspace=0)

    plt.suptitle('Celltype {}, gene_colors {}'.format(hdf5_cell.attrs['annotation'], spot_colors), y=0.92)
    plt.show()
    plt.close()

def plot_cell_3D(cell, spot_colors={}, default_spot_color='grey'):

    fig = plt.figure(figsize=(6,6))
    ax = plt.axes(projection='3d')

    for z_ind in cell.attrs['zslices'][:-1]:
        z = int(z_ind)*1.5
        
        b_xs = cell['boundaries'][z_ind][:,0]
        b_ys = cell['boundaries'][z_ind][:,1]
        border_color = 'gray'
            
        ax.plot3D(b_xs, b_ys, z, border_color)
        
        s_genes = cell['spot_genes'][z_ind][:]
        s_xs = cell['spot_coords'][z_ind][:,0]
        s_ys = cell['spot_coords'][z_ind][:,1]

        #Plot all spots in gray and then paint over if I have colors
        ax.scatter3D(s_xs, s_ys, z, color='gray')

        for gene,color in spot_colors.items():
            gene_inds = s_genes == gene.encode()
            ax.scatter3D(s_xs[gene_inds], s_ys[gene_inds], z, color=color)


    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    """
    for angle in range(0, 360):
        ax.view_init(30, angle)
        plt.draw()
        plt.pause(.001)
    """

    plt.show()
    plt.close()


def plot_cell_spots(hdf5_cell, spot_colors={}, default_spot_color='grey'):
    """
    Plotting function that returns fig,ax

    Takes as input a cell in the form of an hdf5 group
    The hdf5 group must have the following datasets:
    * boundary
    * spot_genes (for coloring)
    * spot_coords
    * attr for the cell centroid

    Plots cell boundary, centroid, and spots

    For coloring spots by gene name accepts an argument spot_colors as a dict
    
    spot_colors = {'Acta2':'red', 'Plekhg3':'blue'} for example
    
    default_spot_color is also an argument for all genes not mentioned in spot_colors

    This function just creates the plots, does not display or save them
    """
    
    fig, ax = plt.subplots(figsize=(8,8))
    
    #Draw the cell outline
    polygon = Polygon(hdf5_cell['boundary'], fill=None)
    ax.add_artist(polygon)
    
    #Add the centroid as a blue spot
    if 'centroid' in hdf5_cell:
        centroid = hdf5_cell['centroid']
    else:
        centroid = hdf5_cell.attrs['centroid']

    c = Circle(
        centroid,
        radius=0.1,
        color='b',
    )
    ax.add_artist(c)
    
    #Add the spots
    colors = []
    for gene in hdf5_cell['spot_genes']:
        if gene.decode() in spot_colors:
            colors.append(spot_colors[gene.decode()])
        else:
            colors.append(default_spot_color)
                
    plt.scatter(
        x = hdf5_cell['spot_coords'][:,0],
        y = hdf5_cell['spot_coords'][:,1],
        alpha = 0.8,
        color = colors,
    )

    plt.axis('equal')
    return fig,ax


def grouper(iterable, n, fillvalue=None):
    """
    Groups an iterator into chunks of size N
    Stolen from https://docs.python.org/3/library/itertools.html
    """
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def ks_scoring(df):
    """
    Takes in a multi-cell dataframe and calculates a per-cell/per-gene score

    Uses a one-sided uniform ks test

    Returns a multi-index series keyed by cell/gene valued on score
    """
    def get_ks_p(obs_ranks):
        stat,p = stats.kstest(
            rvs = obs_ranks,
            cdf = 'uniform',
            alternative='greater',
            mode='exact',
        )
        return p

    #Calculate the per-cell per-gene periphery scores
    ks_ps = df.groupby(['cell_id','gene'])['spot_rank'].apply(get_ks_p)

    return ks_ps

def str2bool(v):
    #stolen from https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')



def f_effect(obs_med, n):
    """
    Calculate a linear effect size
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
        p = ncr(obs_med-1,k-1)*ncr(n-obs_med,m-k)/ncr(n,m)

    else:
        #Even case
        l = m//2
        r = l+1
        

        p = 0
        
        denom = ncr(n,m)
        
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
            p += p_r_given_l*ncr(Rl-1,l-1)*ncr(n-Rl,m-l)/denom

            Rl -= 1
            Rr += 1
    
    return p


def calc_var_effect_size(m,n):
    """
    Calculates the variance in effect sizes for choosing m spots out of n
    Uses the formula below
    
    Var(X) = E[(X-u)**2]
    Var(X) = E[X**2] #since u is always 0
    
    E[X**2] = sum(P(X**2=x**2)x**2)
    
    The effect size can be -x and +x, so I keep a dictionary
    mapping from eff**2 to cumulative probability
    
    Can likely (definitely) be simplified because p(X=-x) = P(X=x))
    """
    if m%2 == 1:
        #odd case
        min_med = m//2+1
        max_med = n-m//2
        exp_med = n/2
        possible_meds = np.arange(min_med,max_med)

    else:
        #even case
        min_med = (m+1)/2
        max_med = n-(m+1)/2+1
        possible_meds = np.arange(min_med,max_med+0.5,0.5)

        
    p_eff = collections.defaultdict(float)
    
    for med in possible_meds:
        eff_sq = f_effect(med,n)**2
        p_eff[eff_sq] += p_med(m,n,med)*eff_sq
            
    var = sum(p_eff.values())
    
    return var


def normal_mean_effect_sizes(df, m_col='m', n_col='n', eff_col='eff', var_memo={}):
    """
    m_n_eff_table is a dataframe with columns of 'm', 'n', and 'eff', the observed effect size
    """
    eff_sum = df[eff_col].sum()
    sn2 = 0
    
    for i,r in df.iterrows():
        m,n = int(r[m_col]),int(r[n_col])
        k = '{}_{}'.format(m,n)
        if k not in var_memo:
            var_memo[k] = calc_var_effect_size(m,n)
            
        sn2 += var_memo[k]
        
    normal_mean = eff_sum/math.sqrt(sn2)
    return normal_mean, var_memo

def map_merf_anns_to_seq(df, seq_source='SS2'):
    merf_to_ss2_ann_mappings = {
        'Astro_1':'Astro',
        'Astro_2':'Astro',
        'Astro_3':'Astro',
        'Endo':'Endo',
        'L23_IT_1':'L2/3 IT',
        'L23_IT_2':'L2/3 IT',
        'L23_IT_3':'L2/3 IT',
        'L23_IT_4':'L2/3 IT',
        'L23_IT_5':'L2/3 IT',
        'L45_IT_1':None,
        'L45_IT_2':None,
        'L45_IT_3':None,
        'L45_IT_4':None,
        'L45_IT_5':None,
        'L45_IT_SSp_1':None,
        'L45_IT_SSp_2':None,
        'L56_NP_1':'L5/6 NP',
        'L56_NP_2':'L5/6 NP',
        'L5_ET_1':'L5 ET',
        'L5_ET_2':'L5 ET',
        'L5_ET_3':'L5 ET',
        'L5_ET_4':'L5 ET',
        'L5_ET_5':'L5 ET',
        'L5_IT_1':'L5 IT',
        'L5_IT_2':'L5 IT',
        'L5_IT_3':'L5 IT',
        'L5_IT_4':'L5 IT',
        'L6_CT_1':'L6 CT',
        'L6_CT_2':'L6 CT',
        'L6_CT_3':'L6 CT',
        'L6_CT_4':'L6 CT',
        'L6_CT_5':'L6 CT',
        'L6_CT_6':'L6 CT',
        'L6_CT_7':'L6 CT',
        'L6_CT_8':'L6 CT',
        'L6_CT_9':'L6 CT',
        'L6_IT_1':'L6 IT',
        'L6_IT_2':'L6 IT',
        'L6_IT_3':'L6 IT',
        'L6_IT_Car3':'L6 IT Car3',
        'L6b_1':'L6b',
        'L6b_2':'L6b',
        'L6b_3':'L6b',
        'Lamp5_1':'Lamp5',
        'Lamp5_2':'Lamp5',
        'Lamp5_3':'Lamp5',
        'Lamp5_4':'Lamp5',
        'Lamp5_5':'Lamp5',
        'Lamp5_6':'Lamp5',
        'Lamp5_7':'Lamp5',
        'Lamp5_8':'Lamp5',
        'Lamp5_9':'Lamp5',
        'Micro_1':None,
        'Micro_2':None,
        'OPC':None,
        'Oligo_1':None,
        'Oligo_2':None,
        'Oligo_3':None,
        'PVM':None,
        'Peri':None,
        'Pvalb_1':'Pvalb',
        'Pvalb_10':'Pvalb',
        'Pvalb_11':'Pvalb',
        'Pvalb_12':'Pvalb',
        'Pvalb_2':'Pvalb',
        'Pvalb_3':'Pvalb',
        'Pvalb_4':'Pvalb',
        'Pvalb_5':'Pvalb',
        'Pvalb_6':'Pvalb',
        'Pvalb_7':'Pvalb',
        'Pvalb_8':'Pvalb',
        'Pvalb_9':'Pvalb',
        'SMC':'SMC',
        'Sncg_1':'Sncg',
        'Sncg_2':'Sncg',
        'Sst_1':'Sst',
        'Sst_2':'Sst',
        'Sst_3':'Sst',
        'Sst_4':'Sst',
        'Sst_5':'Sst',
        'Sst_6':'Sst',
        'Sst_7':'Sst',
        'Sst_8':'Sst',
        'Sst_Chodl':'Sst',
        'VLMC':'VLMC',
        'Vip_1':'Vip',
        'Vip_10':'Vip',
        'Vip_2':'Vip',
        'Vip_3':'Vip',
        'Vip_4':'Vip',
        'Vip_5':'Vip',
        'Vip_6':'Vip',
        'Vip_7':'Vip',
        'Vip_8':'Vip',
        'Vip_9':'Vip',
        'striatum_1':None,
        'striatum_2':None,
        'unannotated':None,
        'ventricle_1':None,
        'ventricle_2':None,
    }

    df['new_annotation'] = df['annotation'].map(merf_to_ss2_ann_mappings)
    df = df[df['new_annotation'].notnull()] #drop the annotations which don't have mappings
    df['annotation'] = df['new_annotation']
    df = df.drop(columns=['new_annotation'])
    

    return df



