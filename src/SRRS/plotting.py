from . import utils

import collections
import pandas as pd
import numpy as np
import pysam

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, PatchCollection
from matplotlib.lines import Line2D
from matplotlib import cm
import seaborn as sns

import sys

def plot_cell_3D(cell, gene_colors={}, color_by_rank=False, default_spot_color='grey', rainbow=False, fig=None, ax=None):
    """
    3D plot of a cell
    optionally color spots by gene or rank

    returns figure handle and ax as tuple
    """
    #If gene colors given and coloring by rank, just color by rank
    gene_colors = {} if color_by_rank else gene_colors

    if not ax or not fig:
        fig = plt.figure(figsize=(6,6))
        ax = plt.axes(projection='3d')

    #plot boundaries and buildup all points
    xs = []
    ys = []
    zs = []
    genes = []
    ranks = []

    for z_ind in cell.zslices:
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
    show_legend = False

    if rainbow:
        #Give each gene a different color
        num_genes = len(np.unique(genes))
        uniq_colors = sns.color_palette('viridis',num_genes)
        gene_colors = {g:c for g,c in zip(sorted(np.unique(genes)),uniq_colors)}
       
        for gene,color in gene_colors.items():
            gene_inds = genes == gene

            if sum(gene_inds):
                ax.scatter3D(xs[gene_inds], ys[gene_inds], zs[gene_inds], color=color, label=gene, s=35)
                show_legend = True

       

    elif color_by_rank:
        #Color all spots by rank
        #If trying to color by rank, but cell is not yet ranked, give a warning
        if not cell.ranked:
            sys.stderr.write('Warning: Trying to color by rank, but cell spots not yet ranked\n')
        else:
            cmap = cm.get_cmap('viridis_r', len(ranks))
            norm = mpl.colors.Normalize(vmin=1, vmax=len(ranks))
            ax.scatter3D(xs, ys, zs, color=cmap(ranks), s=36)
            fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),label='Rank')

    else:
        #Plot spots which have associated colors
        colored_inds = np.array([False]*len(genes))
        for gene,color in gene_colors.items():
            gene_inds = genes == gene
            colored_inds |= gene_inds

            if sum(gene_inds):
                ax.scatter3D(xs[gene_inds], ys[gene_inds], zs[gene_inds], color=color, label=gene, s=100)
                show_legend = True

        #Plot spots without colors in gray (lower zorder to be behind the colored spots)
        ax.scatter3D(
            xs[~colored_inds],
            ys[~colored_inds],
            zs[~colored_inds],
            color='gray', s=10, zorder=-1,
        )

    #ax.view_init(elev=20, azim=0) #changing view angle, can be done on returned ax
    ax.set_axis_off()

    plt.title('Celltype {}'.format(cell.annotation))
    if show_legend:
        plt.legend()


    return fig,ax



def plot_cell_zslices(cell, gene_colors={}, color_by_rank=False, default_spot_color='grey', rainbow=False):
    """
    plot of a cell separated by z-slice
    optionally color spots by gene

    returns figure handle and axs as tuple
    """
    #If gene colors given and coloring by rank, just color by rank
    gene_colors = {} if color_by_rank else gene_colors

    #Make as square a layout as I can
    zslices = cell.zslices
    nslices = len(zslices)
    nrows = int(np.sqrt(nslices))
    ncols = nrows
    while nrows*ncols < nslices:
        ncols += 1

    fig, axs = plt.subplots(figsize=(5*ncols,5*nrows),nrows=nrows,ncols=ncols,sharex=True,sharey=True)
    try:
        axs = axs.flatten()
    except AttributeError:
        #if the grid is 1x1 then an axis is directly returned
        axs = [axs]


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

        #give each gene a different color
        elif rainbow:
            rainbow_colors = sns.color_palette("hls", len(cell.genes))
            rainbow_map = {gene:color for gene,color in zip(cell.genes,rainbow_colors)}

            colors = []
            for gene in cell.spot_genes[zslice]:
                colors.append(rainbow_map[gene])
                if gene not in labels:
                    labels.append(gene)
                    handle = mpatches.Circle((0.5, 0.5), radius = 0.25, facecolor=rainbow_map[gene], edgecolor="none")
                    handles.append(handle)



        #otherwise make spots all default
        else:
            colors = [default_spot_color]*len(cell.spot_coords[zslice])

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

        #Determine which indices are colored/uncolored
        colors = np.array(colors)
        uncolored_inds = np.array([c == default_spot_color for c in colors])
        colored_inds = ~uncolored_inds

        #Plot uncolored spots smaller
        axs[i].scatter(
            x = cell.spot_coords[zslice][uncolored_inds,0],
            y = cell.spot_coords[zslice][uncolored_inds,1],
            alpha = 0.8,
            color = colors[uncolored_inds],
            s = 36,
        )

        #Plot colored spots larger
        axs[i].scatter(
            x = cell.spot_coords[zslice][colored_inds,0],
            y = cell.spot_coords[zslice][colored_inds,1],
            alpha = 0.8,
            color = colors[colored_inds],
            s = 72,
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

    return fig,axs




def plot_tissue_level(cells, color_by_score_gene=None, color_by_ontology=False, region=(False,False,False,False)):
    """
    Plot cells across tissues for a zoomed out view of patterns
    Color by score or by cell-type

    Draws outlines from just a single zslice, the first for each cell, which shouldn't matter in a zoomed out plot

    optionally takes input of boundaries which is a 4-tuple of (xmin,ymin,xmax,ymax) and only plots cell fully in
    this rectangular region

    if not caring about a certain boundary such as min_x, set to False
    """
    r_minx,r_miny,r_maxx,r_maxy = region

    boundaries = []
    cmap_vals = []

    fig, ax = plt.subplots(figsize=(6,6))

    for cell in cells:
        z = cell.zslices[0] #just plot the "top" z-slice
        bounds = cell.boundaries[z]
        xs,ys = bounds[:,0],bounds[:,1]

        #check if a boundary is set (not False) and the cell exceeds it
        outside_region = any([
            r_minx and (min(xs) < r_minx),
            r_miny and (min(ys) < r_miny),
            r_maxx and (max(xs) > r_maxx),
            r_maxy and (max(ys) > r_maxy),
        ])

        if outside_region:
            continue

        #determine coloring
        if color_by_ontology:
            cmap_vals.append(cell.annotation)
        elif color_by_score_gene:
            med_rank = cell.gene_med_ranks.get(color_by_score_gene)
            cmap_val = utils.score(med_rank,cell.n) if med_rank else np.nan
            cmap_vals.append(cmap_val)

        boundaries.append(bounds)


    # Make the collection and add it to the plot
    # do this differently depending on how the coloring is done
    if color_by_ontology:
        uniq_onts = sorted(list(set(cmap_vals)))
        colors = sns.color_palette("hls", len(uniq_onts))
        ont_to_color = {ont:color for ont,color in zip(uniq_onts,colors)}
        colors = [ont_to_color[ont] for ont in cmap_vals]
        coll = PolyCollection(boundaries, facecolors=colors, edgecolors='none')

        markers = [
                Line2D([0], [0], marker='o', color='w', markerfacecolor=ont_to_color[ont], label=ont)
                for ont in uniq_onts
        ]
        ax.legend(handles=markers)

    elif color_by_score_gene:
        cmap = cm.get_cmap('coolwarm')
        cmap.set_bad(color='black')
        coll = PolyCollection(boundaries, array=cmap_vals, cmap=cmap, edgecolors='none')
        fig.colorbar(coll, ax=ax)


    else:
        coll = PolyCollection(boundaries, facecolor='none', edgecolor='black')

    ax.add_collection(coll)
    ax.autoscale_view()
    plt.axis('equal')
    plt.axis('off')

    return fig,ax


def read_buildup_plot(bam_path, locus, ann_df, stratify_tag=None, **kwargs):
    """
    Creates a rough read buildup plot with gene annotations
    useful for creating many plots but not figure quality

    returns a matplotlib axis object

    Arguments:
        bam_path [required]
            path to a position sorted bam with bam.bai index
            read counts can be grouped by custom bam flag such as cell-type

        locus (chrom, start, end) [required]
            these three arguments determine which genomic window to visualize
            example locus = ('chr10', 86344341, 86350020)

        ann_df [required]
            pd.DataFrame object created from gtf/gff file that has the following columns
            chrom/kind/start/end/label/group
            examples of labels are gene names, and groups are transcript_ids

        stratify_tag [optional]
            a bam tag to use to stratify the reads
            for example using 'CB' would create a separate line plot for each cell
            more realistic is using 'XO' which stratifies by cell-type if that tag is present

        **kwargs
            ws: window_size, default 10
    """
    #handling locus
    try:
        chrom,start,end = locus
        start = int(start)
        end = int(end)
    except:
        raise Exception('Must pass in a tuple of locus info locus=("chr1",1,100)')


    #handling kwargs
    ws = kwargs.get('ws',10)
    min_ont_counts = kwargs.get('min_ont_counts',0)

    #Counting reads per tag in this region
    def get_tag(read):
        try:
            return read.get_tag(stratify_tag)
        except KeyError:
            return None

    strat_func = get_tag if stratify_tag else (lambda read: 'All')

    pos_counts = collections.defaultdict(lambda: {p:0 for p in range(start//ws*ws,end//ws*ws+ws,ws)})
    tot_counts = 0

    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch(chrom,start,end):
            if r.pos < start or r.pos > end:
                continue

            strat = strat_func(r)
            if strat:
                pos_counts[strat][r.pos//ws*ws] += 1
                tot_counts += 1

    count_df = (
        pd.DataFrame(pos_counts)
            .rename_axis('Window')
            .reset_index()
            .melt(
                id_vars='Window',
                var_name='Ontology',
                value_name='Count',
            )
    )

    count_df['CDF'] = (
        count_df.groupby('Ontology')['Count']
            .apply(lambda v: v.cumsum()/v.sum())
    )
   
    count_df = count_df.groupby('Ontology').filter(lambda g: g['Count'].sum() >= min_ont_counts)
    onts = count_df.groupby('Ontology')['CDF'].apply(lambda v: sum(v<0.5)).sort_values().index
    ont_colors = {o:c for o,c in zip(onts,sns.color_palette("hls", len(onts)))}


    #Subsetting the gene annotation info to this region
    genes = []
    if len(ann_df) > 0:
        ann_df = ann_df[
            ann_df['chrom'].eq(chrom) &
            (ann_df['start'].between(start,end) | ann_df['end'].between(start,end))
        ]

        genes = ann_df['label'].unique()
        gene_colors = {g:c for g,c in zip(genes,sns.color_palette("hls", len(genes)))}

    #Creating subplots depending on the number of ontologies
    #first row is the gene annotation
    #rows 2 -> n+1 are the n different ontologies
    num_rows = 3+len(onts)
    fig,axs = plt.subplots(
        figsize=(5,num_rows),
        nrows=num_rows, ncols=1,
        sharex=True, sharey=False,
    )

    ann_ax = axs[0]
    ann_ax.axis('off')
    count_ax = axs[1]
    sum_ax = axs[2]

    #Drawing the gene annotations exons/UTRs
    for i,((gene,transcript_id),transcript_df) in enumerate(ann_df.groupby(['label','group'])):
        transcript_df = ann_df[ann_df['label'].eq(gene) & ann_df['group'].eq(transcript_id)]

        y = 2*(i+1)
        #Draw a line for the intron from the total start to total end
        min_x_transcript = transcript_df['start'].min()
        max_x_transcript = transcript_df['end'].max()

        intron_line = mpatches.Rectangle(
            (min_x_transcript, y - 0.1),
            max_x_transcript - min_x_transcript, 0.2,
            linewidth=1,edgecolor='k',facecolor='k',
        )
        ann_ax.add_patch(intron_line)

        for _,feature in transcript_df.iterrows():

            if 'UTR' in feature['kind']:
                height = 0.5

            elif feature['kind'] == 'exon':
                height = 1

            else:
                continue
        
            ymax = y+1
            x = feature['start']
            width = feature['end']-feature['start']

            rect = mpatches.Rectangle(
                (x, y - height/2),
                width, height,
                linewidth=0.5,edgecolor='k',facecolor=gene_colors[feature['label']],
            )
            ann_ax.add_patch(rect)

        ann_ax.set_ylim(-1,y+1)

        #Add legend on top of plot
        handles = [Line2D([0],[0], color=gene_colors[g], lw=4) for g in genes]
        ann_ax.legend(
            handles, genes,
            title='Gene names',
            loc='upper center',
            bbox_to_anchor=(-0.2, 0.5), ncol=len(genes),
        )

    #Plot individual density plots
    for i,ont in enumerate(onts):
        plot_df = count_df[count_df['Ontology'].eq(ont)]
        sns.lineplot(
            x = 'Window',
            y = 'Count',
            color = ont_colors[ont],
            data = plot_df,
            ax = axs[i+3],
        )
        axs[i+3].set_ylabel(ont)
       

    #Plot combined density and cumsum plots
    sns.lineplot(
        x = 'Window',
        y = 'Count',
        hue = 'Ontology',
        hue_order = onts,
        palette = ont_colors,
        data = count_df, 
        ax = count_ax,
    )

    sns.lineplot(
        x = 'Window',
        y = 'CDF',
        hue = 'Ontology',
        hue_order = onts,
        palette = ont_colors,
        data = count_df,
        legend = False,
        ax = sum_ax,
    )
    count_ax.legend(loc='center left', title='Cell-type', bbox_to_anchor=(1, 0.5))
    count_ax.set_xlim(start,end)
    sum_ax.set_xticks([])
    
    return fig


