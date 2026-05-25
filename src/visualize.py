"""

Module for visualizing QC, HVGs, PCA, UMAP, Leiden clustering 

"""

from matplotlib.pylab import save
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_qc_metrics(adata, save = False, figdir='../results/'): 
    """
    Plot QC metrics for the AnnData object.

    Parameters:
    adata (AnnData): The input AnnData object.
    save (bool): Whether to save the plot as a PNG file.
    figdir (str): Directory to save the plot if save is True.

    """
    sc.pl.violin(
        adata,
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        jitter=0.4,
        multi_panel=True,
        save = '_qc_metrics.png' if save else None
    )

def plot_hvgs(adata, save = False, figdir='../results/'):
    """
    Plot highly variable genes (HVGs) for the AnnData object.

    Parameters:
     adata (AnnData): The input AnnData object.
    save (bool): Whether to save the plot as a PNG file.
    figdir (str): Directory to save the plot if save is True.

    """
    sc.pl.highly_variable_genes(
        adata,
        save = '_hvgs.png' if save else None
    )

def plot_pca(adata, save = False, figdir='../results/'):
    """
    Plot PCA results for the AnnData object.

    Parameters:
    adata (AnnData): The input AnnData object.
    save (bool): Whether to save the plot as a PNG file.
    figdir (str): Directory to save the plot if save is True.

    """
    sc.pl.pca(
        adata,
        log = True,
        save = '_pca.png' if save else None
    )

def plot_umap(adata, color = 'leiden', save = False, figdir='../results/'):
    """
    Plot UMAP results for the AnnData object.

    Parameters:
    adata (AnnData): The input AnnData object.
    color (str): The column in adata.obs to use for coloring the UMAP plot.
    save (bool): Whether to save the plot as a PNG file.
    figdir (str): Directory to save the plot if save is True.

    """
    sc.pl.umap(
        adata,
        color = color,
        legend_loc = 'on data',
        save = '_umap.png' if save else None
    )

def plot_dotplot(adata, marker_genes, groupby='leiden', save=False, figdir='../results/'):
    """
    Plot a dotplot for the specified marker genes.

    Parameters:
    adata (AnnData): The input AnnData object.
    marker_genes (list): A list of marker genes to include in the dotplot.
    save (bool): Whether to save the plot as a PNG file.
    figdir (str): Directory to save the plot if save is True.

    """
    sc.pl.dotplot(
        adata,
        var_names = marker_genes,
        groupby = groupby,
        save = '_dotplot.png' if save else None
    )

def plot_marker_genes(adata, n_genes = 10, save = False, figdir='../results/'):
    """
    Plot the top marker genes for each cluster.

    Parameters:
    adata (AnnData): The input AnnData object.
    n_genes (int): The number of top marker genes to plot for each cluster.
    save (bool): Whether to save the plot as a PNG file.
    figdir (str): Directory to save the plot if save is True.

    """
    sc.pl.rank_genes_groups(
        adata,
        n_genes = n_genes,
        sharey = False,
        save = '_marker_genes.png' if save else None
    )
