import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns

import sys

def plot_cell_3D(cell, gene_colors={}, color_by_rank=False, default_spot_color='grey'):
    """
    3D plot of a cell
    optionally color spots by gene or rank

    shows plot and also returns figure handle
    """
    #If gene colors given and coloring by rank, just color by rank
    gene_colors = {} if color_by_rank else gene_colors

    fig = plt.figure(figsize=(6,6))
    ax = plt.axes(projection='3d')

    #plot boundaries and buildup all points
    xs = []
    ys = []
    zs = []
    genes = []
    ranks = []

    for z_ind in cell.zslices[:-1]:
        z = int(z_ind)*1.5

        b_xs = cell.boundaries[z_ind][:,0]
        b_ys = cell.boundaries[z_ind][:,1]
        border_color = 'gray'

        ax.plot3D(b_xs, b_ys, z, border_color)

        s_xs = cell.spot_coords[z_ind][:,0]
        s_ys = cell.spot_coords[z_ind][:,1]
        s_zs = [z]*len(s_xs)
        s_genes = cell.spot_genes[z_ind]
        s_ranks = cell.spot_ranks[z_ind]

        xs.extend(s_xs)
        ys.extend(s_ys)
        zs.extend(s_zs)
        genes.extend(s_genes)
        ranks.extend(s_ranks)

    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)
    genes = np.array(genes)

    #Plot all spots in gray and then paint over if I have colors
    ax.scatter3D(xs, ys, zs, color='gray')
    show_legend = False
    for gene,color in gene_colors.items():
        gene_inds = genes == gene

        if sum(gene_inds):
            ax.scatter3D(xs[gene_inds], ys[gene_inds], zs[gene_inds], color=color, label=gene)
            show_legend = True

    #If trying to color by rank, but cell is not yet ranked, give a warning
    if color_by_rank:
        if not cell.ranked:
            sys.stderr.write('Warning: Trying to color by rank, but cell spots not yet ranked\n')
        else:
            cmap = cm.get_cmap('viridis_r', len(ranks))
            norm = mpl.colors.Normalize(vmin=1, vmax=len(ranks))
            ax.scatter3D(xs, ys, zs, color=cmap(ranks))
            fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),label='Rank')


    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    plt.title('Celltype {}'.format(cell.annotation))
    if show_legend:
        plt.legend()

    plt.show()
    plt.close()

    return fig



def plot_cell_zslices(cell, gene_colors={}, color_by_rank=False, default_spot_color='grey'):
    """
    plot of a cell separated by z-slice
    optionally color spots by gene

    shows plot and also returns figure handle
    """
    #If gene colors given and coloring by rank, just color by rank
    gene_colors = {} if color_by_rank else gene_colors

    zslices = cell.zslices
    fig, axs = plt.subplots(figsize=(8,8),nrows=3,ncols=3,sharex=True,sharey=True)
    axs = axs.flatten()

    global_min_x,global_min_y = None,None
    global_max_x,global_max_y = None,None

    handles = []
    labels = []
    cmap = cm.get_cmap('viridis_r', cell.n)
    norm = mpl.colors.Normalize(vmin=1, vmax=cell.n)

    for i,zslice in enumerate(zslices):
        #default coloring
        colors = [default_spot_color]*cell.n_per_z[zslice]

        #color spots by gene
        if gene_colors:
            colors = []
            for gene in cell.spot_genes[zslice]:
                if gene in gene_colors:
                    colors.append(gene_colors[gene])
                    if gene not in labels:
                        labels.append(gene)
                        handle = mpatches.Circle((0.5, 0.5), radius = 0.25, facecolor=gene_colors[gene], edgecolor="none")
                        handles.append(handle)
                else:
                    colors.append(default_spot_color)

        #color by rank
        elif color_by_rank and cell.ranked:
            colors = cmap(cell.spot_ranks[zslice])


        #Draw the cell outline
        boundary = cell.boundaries[zslice][:]
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

        polygon = mpatches.Polygon(boundary, fill=None)
        axs[i].add_artist(polygon)

        axs[i].scatter(
            x = cell.spot_coords[zslice][:,0],
            y = cell.spot_coords[zslice][:,1],
            alpha = 0.8,
            color = colors,
        )


    #If trying to color by rank, but cell is not yet ranked, give a warning
    if color_by_rank:
        if not cell.ranked:
            sys.stderr.write('Warning: Trying to color by rank, but cell spots not yet ranked\n')
        else:
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
            fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),label='Rank', cax=cbar_ax)

    for used_i in range(i+1):
        axs[used_i].set_xticks([])
        axs[used_i].set_yticks([])
        axs[used_i].set_xlim(global_min_x,global_max_x)
        axs[used_i].set_ylim(global_min_y,global_max_y)


    for unused_i in range(i+1,len(axs)):
        axs[unused_i].set_axis_off()

    plt.subplots_adjust(wspace=0, hspace=0)

    plt.suptitle('Celltype {}'.format(cell.annotation), y=0.92)
    if handles and labels:
        fig.legend(handles,labels)
    plt.show()
    plt.close()

    return fig

