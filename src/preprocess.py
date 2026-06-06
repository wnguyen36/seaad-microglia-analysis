"""
Module for preprocessing single-cell RNA-seq data using Scanpy. This module includes functions for quality control, normalization, and checkpointing of AnnData objects.

"""
import scanpy as sc
import os

# Quality control on the AnnData object
def run_qc(adata, min_genes=200, max_genes=2500, max_mt=5, verbose = True):
    """
    Run quality control on the AnnData object.

    Parameters:
    adata (AnnData): The input AnnData object.
    min_genes (int): Minimum number of genes expressed in a cell to be retained.
    max_genes (int): Maximum number of genes expressed in a cell to be retained. (removes doublets)
    max_mito (float): Maximum percentage of mitochondrial genes expressed in a cell to be retained. (removes dead/dying cells)

    Returns:
    AnnData: The filtered AnnData object after quality control.
    """
    if verbose:
        print(adata)
        print(f"Running QC: {adata.n_obs} cells, {adata.n_vars} genes")
        
    # Flagging MT
    adata.var['mt'] = adata.var_names.str.startswith('MT-')

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Apply filters
    adata = adata[adata.obs.n_genes_by_counts > min_genes]
    adata = adata[adata.obs.n_genes_by_counts < max_genes]
    adata = adata[adata.obs.pct_counts_mt < max_mt]

    # Filter genes that are expressed in less than 3 cells
    sc.pp.filter_genes(adata, min_cells=3)

    if verbose:
        print(f"After QC: {adata.n_obs} cells, {adata.n_vars} genes")
        print(f"Removed: {adata.n_obs} cells")
    
    return adata


# Normalize the data
def run_normalization(adata, verbose = True): 
    """
    Normalize, log transform, select variable genes, regress, scale data
    
    """

    if verbose:
        print("Running normalization...")

    # Normalize total counts per cell
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Logarithmize the data
    sc.pp.log1p(adata)

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    if verbose: 
        print(f"Identified {adata.var.highly_variable.sum()} highly variable genes.")

    # Regress out effects of total counts and percentage of mitochondrial genes
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

    # Scale the data to unit variance and zero mean
    sc.pp.scale(adata, max_value=10)

    if verbose:
        print("Normalization complete.")
    
    return adata

def save_checkpoint(adata, path, verbose = True): 
    """
    Save the AnnData object as an h5ad file.

    Parameters:
    adata (AnnData): The AnnData object to be saved.
    path (str): The file path where the AnnData object will be saved.
    verbose (bool): Whether to print a message after saving.

    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    adata.write_h5ad(path)
    
    if verbose:
        print(f"Checkpoint saved at: {path}")
        
def load_checkpoint(path, verbose = True): 
    """
    Load an AnnData object from an h5ad file.

    Parameters:
    path (str): The file path from which to load the AnnData object.
    verbose (bool): Whether to print a message after loading.

    """
    adata = sc.read_h5ad(path)
    
    if verbose:
        print(f"Checkpoint loaded from: {path}")
    
    return adata