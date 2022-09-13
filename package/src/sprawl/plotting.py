from . import utils

import collections
import pandas as pd
import scipy as scp
import numpy as np
import pysam

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, PatchCollection
from matplotlib.lines import Line2D
from matplotlib import cm
import matplotlib.backends.backend_pdf
import seaborn as sns

import sys

def plot_cell_3D(cell, gene_colors={}, default_spot_color='grey', rainbow=False, fig=None, ax=None):
    """
    3D plot of a cell
    optionally color spots by gene

    returns figure handle and ax as tuple
    """
    if not ax or not fig:
        fig = plt.figure(figsize=(6,6))
        ax = plt.axes(projection='3d')

    #plot boundaries and buildup all points
    xs = []
    ys = []
    zs = []
    genes = []

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

        xs.extend(s_xs)
        ys.extend(s_ys)
        zs.extend(s_zs)
        genes.extend(s_genes)

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
        plt.legend(ncol=3)


    return fig,ax



def plot_cell_zslices(cell, gene_colors={}, default_spot_color='grey', rainbow=False):
    """
    plot of a cell separated by z-slice
    optionally color spots by gene

    returns figure handle and axs as tuple
    """
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


def read_buildup_plot(bam_path, locus, ann_df, spatial_df=None, stratify_tag=None, **kwargs):
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

    count_df = utils.bam_read_positions(bam_path, locus, stratify_tag=stratify_tag, **kwargs)

    #handling kwargs
    ws = kwargs.get('ws',10)
    strand = kwargs.get('strand',None)

    bins = np.arange(start,end+ws,ws)
    labels = bins[:-1]


    count_df['pos'] = pd.cut(
        count_df['pos'],
        bins=bins,
        labels=labels,
    )
    count_df = count_df.groupby(['strat','pos']).size().reset_index(name='Count')
    count_df = count_df.rename(columns={'strat':'Ontology','pos':'Window'})

    count_df['CDF'] = (
        count_df.groupby('Ontology')['Count']
            .apply(lambda v: v.cumsum()/v.sum())
    )

    #if there is spatial data provided then subset to just shared ontologies
    if spatial_df is not None:
        shared_onts = set(count_df['Ontology']).intersection(spatial_df['annotation'])
        if not shared_onts:
            #sys.stderr.write('Warning: No shared ontologies between spatial_df and count_df, returning None\n')
            return None

        spatial_df = spatial_df[spatial_df['annotation'].isin(shared_onts)]
        count_df = count_df[count_df['Ontology'].isin(shared_onts)]
  
    onts = count_df.groupby('Ontology')['CDF'].apply(lambda v: sum(v<0.5)).sort_values().index
    ont_colors = {o:c for o,c in zip(onts,sns.cubehelix_palette(n_colors=len(onts),start=.5, rot=-.75))}


    #Subsetting the gene annotation info to this region
    genes = []
    if len(ann_df) > 0:
        ann_df = ann_df[
            ann_df['chrom'].eq(chrom) &
            (ann_df['start'].between(start,end) | ann_df['end'].between(start,end))
        ]

        if strand:
            ann_df = ann_df[ann_df['strand'].eq(strand)]

        genes = ann_df['label'].unique()
        gene_colors = {g:c for g,c in zip(genes,sns.color_palette("hls", len(genes)))}

    #Creating subplots depending on the number of ontologies
    #and whether or not there is are spatial data to plot alongside
    #first row is the gene annotation
    #rows 2 -> n+1 are the n different ontologies
    ncols = 1 if spatial_df is None else 2

    num_rows = 3+len(onts)
    fig,axs = plt.subplots(
        figsize=(5*ncols,num_rows),
        nrows=num_rows, ncols=ncols,
        sharex=False, sharey=False,
    )

    ann_ax = axs[0] if spatial_df is None else axs[0][1]
    ann_ax.axis('off')
    count_ax = axs[1] if spatial_df is None else axs[1][1]
    sum_ax = axs[2] if spatial_df is None else axs[2][1]

    #Drawing spatial_df if present
    if spatial_df is not None:
        #merging the top left axes into a larger axis
        gs = axs[0, 0].get_gridspec()
        for ax in axs[[0,1], 0]:
            ax.remove()

        corr_ax = fig.add_subplot(gs[:2, 0])
        axs[2][0].axis('off')

        mean_spatial = spatial_df.groupby('annotation')['score'].mean()

        exp_val = (start+end)/2
        half_span = (end-start)/2
        count_df['Window'] = count_df['Window'].astype(int)

        #Calculating median ReadZs proxy score
        mean_readzs = count_df.groupby('Ontology').apply(
            lambda g: ((g['Window'].multiply(g['Count']).sum()/g['Count'].sum())-exp_val)/half_span
        )

        r,p = scp.stats.pearsonr(mean_spatial, mean_readzs)
        corr_ax.set_title('Pearson r={:.2f} p={:.2f}'.format(r,p))

        corr_df = pd.DataFrame({
            'Spatial':mean_spatial,
            'Genomic':mean_readzs,
        }).reset_index()
        corr_df = corr_df.rename(columns={'index':'ont'})
        
        sns.regplot(
            x = 'Genomic',
            y = 'Spatial',
            color = 'grey',
            ci = None,
            scatter_kws = {'s':0},
            line_kws = {'linestyle':'dashed'},
            data = corr_df,
            ax = corr_ax,
        )

        sns.scatterplot(
            x = 'Genomic',
            y = 'Spatial',
            hue = 'ont',
            palette = ont_colors,
            data = corr_df,
            legend = False,
            ax = corr_ax,
        )

        #Create each spatial plot for the different onts
        for i,ont in enumerate(onts):
            plot_ax = axs[i+3][0]
            sns.boxplot(
                x = 'score',
                y = 'annotation',
                color = ont_colors[ont],
                data = spatial_df[spatial_df['annotation'].eq(ont)],
                ax = plot_ax,
            )
            plot_ax.axvline(0, linestyle='dashed', color='grey')
            plot_ax.set_xlim(-1,1)
            plot_ax.set_xlabel('')
            plot_ax.set_ylabel(ont)
            plot_ax.set_xticks([])
            plot_ax.set_yticks([])

        #Add an xlabel to the last subplot
        plot_ax.set_xticks([-1,0,1])
        plot_ax.set_xlabel('Spatial score')

    #Drawing the gene annotations exons/UTR
    #the following if/else sorts the transcripts to look pretty
    if not strand or strand == '+':
        ann_df['extreme'] = ann_df.groupby(['label','group'])['end'].transform('max')
    else:
        ann_df['extreme'] = ann_df.groupby(['label','group'])['start'].transform('min')

    ann_df = ann_df.sort_values('extreme')

    for i,((gene,transcript_id),transcript_df) in enumerate(ann_df.groupby(['label','group'],sort=False)):
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
    ann_ax.set_xlim(start,end)
    ann_ax.set_xticks([])

    #Add legend on top of plot
    handles = [Line2D([0],[0], color=gene_colors[g], lw=4) for g in genes]
    ann_ax.legend(
        handles, genes,
        title='Gene names' if not strand else f'Gene names for {strand} strand',
        loc='upper center',
        bbox_to_anchor=(0.5, 2.0), 
        ncol=len(genes),
    )

    #Plot individual density plots
    for i,ont in enumerate(onts):
        plot_ax = axs[i+3] if spatial_df is None else axs[i+3][1]
        plot_df = count_df[count_df['Ontology'].eq(ont)]

        sns.lineplot(
            x = 'Window',
            y = 'Count',
            color = ont_colors[ont],
            data = plot_df,
            ax = plot_ax,
        )
        plot_ax.set_ylabel(ont)
        plot_ax.set_xlim(start,end)
        plot_ax.set_xticks([])
        plot_ax.set_xlabel('')

    #Add an xlabel to the last subplot
    plot_ax.set_xlabel(chrom)
    plot_ax.set_xticks([start,end])
  

    #Plot combined density and cumsum plots
    sns.lineplot(
        x = 'Window',
        y = 'Count',
        hue = 'Ontology',
        hue_order = onts,
        palette = ont_colors,
        data = count_df, 
        legend = False,
        ax = count_ax,
    )
    count_ax.set_xlim(start,end)
    count_ax.set_xticks([])
    count_ax.set_xlabel('')

    sns.lineplot(
        x = 'Window',
        y = 'CDF',
        hue = 'Ontology',
        hue_order = onts,
        palette = ont_colors,
        data = count_df,
        ax = sum_ax,
    )
    sum_ax.axhline(0.5, linestyle='dashed', color='grey')
    sum_ax.legend(loc='center left', title='Cell-type', bbox_to_anchor=(1, 0.5))
    sum_ax.set_xlim(start,end)
    sum_ax.set_xticks([])
    sum_ax.set_yticks([0, 0.5, 1])
    sum_ax.set_xlabel('')

    
    return fig


def make_pdf(out_name):
    """
    Silly convenience function because I always forget how to make a pdf handle to save figures to
    (responsibility of the caller to close the handle)
    """
    pdf = matplotlib.backends.backend_pdf.PdfPages(out_name)
    return pdf


