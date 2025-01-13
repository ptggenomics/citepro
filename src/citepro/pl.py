import scanpy as sc
import matplotlib.pyplot as plt
from anndata import AnnData
import numpy as np

from .utils import find_encoding_rna, find_feature


def multiomics_feature_plot(adata: AnnData, protein:str, basis:str= 'X_umap_rna', projection:str='2d'):
    
    rna_name = find_encoding_rna(adata, protein)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    protein_name = protein if protein in adata.var.index else None
    if protein_name: 
        prot_level = find_feature(adata, protein_name)
        vmin = max(2, np.percentile(prot_level, 0.5))
        vmax = min(6, np.percentile(prot_level, 99.5))
        sc.pl.embedding(adata, basis= basis, color = protein, projection=projection, ax=axes[0], show=False, cmap='magma_r',vmin = vmin, vmax = vmax)
    else:
        axes[0].text(0.35, 0.5, f'No such protein {protein}', c='#7777ff')
        axes[0].set_axis_off()

    if rna_name:   
        sc.pl.embedding(adata, basis= basis, color = rna_name, projection=projection, ax=axes[1], show=False, cmap='viridis_r')
    else:
        axes[1].text(0.25, 0.5, f'No encoding RNA found for {protein}', c='#7777ff')
        axes[1].set_axis_off()
    fig.tight_layout()
    fig.show()