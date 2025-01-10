import scanpy as sc
import muon as mu
from mudata import MuData

import numpy as np

from typing import Literal, Optional
from ..io import read_10x_filter
from ..pp import arcsinh_transform, calc_qc_prot_var, calc_cell_qc, gen_adata_celltypist

from tqdm import trange
import logging
logger = logging.getLogger()
"""create_mu module contains code to convert cellranger count matrix into MuData format
"""

def _cluster_umap_by_modal(mudat: MuData, modal = Literal['prot', 'rna'], add_3d:bool = True):
    logger.info(f"{modal} - Run PCA")
    sc.pp.pca(mudat[modal])
    logger.info(f"{modal} - Generate neighbor graph")
    sc.pp.neighbors(mudat[modal])
    logger.info(f"{modal} - Generate UMAP")
    sc.tl.umap(mudat[modal])
    if add_3d:
        logger.info(f"{modal} - Generate 3D UMAP")
        tmp_adata_3d = mudat[modal].copy()
        sc.tl.umap(tmp_adata_3d, n_components=3)
        mudat[modal].obsm['X_umap_3d'] = tmp_adata_3d.obsm['X_umap'].copy()
        del tmp_adata_3d

    for i in trange(5,11, desc=f"leiden cluster for {modal}:"):
        res = i/10
        tmp_adat = mudat[modal].copy()
        sc.tl.leiden(tmp_adat,resolution =res)
        clust_name = 'clust_leiden_{}_{}'.format(modal, str(tmp_adat.uns['leiden']['params']['resolution']))
        mudat[modal].uns[clust_name] = tmp_adat.uns['leiden']
        mudat[modal].obs[clust_name] = tmp_adat.obs['leiden']
        del tmp_adat
    


def create_mudata(path_count: str, 
         path_map_rna:Optional[str] = None, 
         samp_id: Optional[str] = None, 
         prot_norm: Literal['clr','asinh', 'none'] = 'asinh',
         allow_file:Optional[str]=None, 
         block_file:Optional[str]=None,
         celltypist_model:Optional[str]=None, 
         add_3d_umap = True) -> MuData:
    """
    Create MuData object from 10x Genomics Cellranger count matrix, with ability to accomodate barcode allow/block list.

    Parameters
    ----------
    path_count : str
        path to the sample_feature_bc_matrix (h5 or mtx format, .h5 is recommended).
    path_map_rna : Optional[str], optional
        path to the feature_reference.csv file used in the cellranger multi. Must contain the map_rna column. [default: None]
    samp_id : Optional[str], optional
        sample id for this count matrix, this is very useful if multiple count matrices will be combined in downstream analysis. [default: None]
    prot_norm : str, optional
        normalization/transformation method for protein. Options are 'clr' (Centered Log Ratio), 'asinh' (Jittered asinh transformation), or 'none'. [default: 'asinh']
    allow_file : Optional[str], optional
        path to barcode allow list csv file. barcodes in this list will be kept. Can be generated with Loupe Browser. [default: None]
    block_file : Optional[str], optional
        path to barcode block list csv file. barcodes in this list will be removed. Can be generated with Loupe Browser. [default: None]
    celltypist_model : Optional[str], optional
        Celltypist models to be used, 'Immune_All_Low.pkl' is a common choice is the sample type is PBMC. If 'None' is specified, celltypist prediction will not be performed. [default: None]
    add_3d_umap : bool, optional
        Whether to generate 3D UMAP projection. [default: True]

    Returns
    -------
    MuData
        pre-processed multi-omic dataset including cell- and feature-oriented qc matrices, normalized counts, celltype prediction (if specified), and 2D/3D umap projections. 
    """    

    ## Read 
    mudat = read_10x_filter(path_count=path_count, path_map_rna=path_map_rna,
                            allow_file=allow_file, block_file=block_file)
    
    ## Step 1 calculate all the stastics with raw integer counts
    logger.info('Claculating protein UMI count descriptive metadata')
    calc_qc_prot_var(mudat['prot'])

    logger.info('Claculating UMI count descriptive metadata for each cells')
    calc_cell_qc(mudat)
    
    ## Step 2 run celltypist if needed
    if celltypist_model:
        logger.info(f'Predicting celltype using Celltypist with model {celltypist_model}')
        gen_adata_celltypist(mudat, inplace = True)
    

    ## Step 3 add 'sample_id' into data if supplied
    if samp_id:
        logger.info(f'Attaching Sample ID {samp_id}')
        mudat.obs['sample_id'] = samp_id
        for md_key in mudat.mod.keys():
            mudat[md_key].obs['sample_id'] = mudat.obs['sample_id']
    else:
        logger.info('Sample ID not specified')
    
    ## Step 4 normalize the protein and rna data

    match prot_norm:
        case 'clr':
            #print('Running clr')
            logger.info('Protein - Performing CLR.')
            mu.prot.pp.clr(mudat['prot'], inplace = True)
            mudat['prot'].X = np.array(mudat['prot'].X.todense())
            mudat['prot'].layers['clr'] = mudat['prot'].X.copy()
        case 'asinh':
            logger.info('Protein - Performing asinh transformation with jittering.')
            arcsinh_transform(mudat['prot'], densify=True ,to_layers=False)
            mudat['prot'].layers['arcsinh'] = mudat['prot'].X.copy()
        case _:
            logger.info('No normalization for protein.')

    rna_target_sum = 1e4
    logger.info(f'RNA - Performing log1p normalizing for rna with target total count per cell of {rna_target_sum}.')
    sc.pp.normalize_total(mudat['rna'],target_sum=1e4)
    sc.pp.log1p(mudat['rna'])
    mudat['rna'].layers['norm_log1p'] = mudat['rna'].X.copy()
    
    
    ## Step 5 cluster and umap

    for mdkey in mudat.mod.keys():
        _cluster_umap_by_modal(mudat, mdkey, add_3d = add_3d_umap)
    
    ## WNN integration
    logger.info('Integration - Generate integrated protein+rna umap using WNN')
    logger.info("WNN - Generate neighbor graph")
    mu.pp.neighbors(mudat)
    logger.info("WNN - Generate UMAP")
    adat_wnn = mu.tl.umap(mudat, copy = True)
    mudat.obsm['X_umap_wnn'] = adat_wnn.obsm['X_umap'].copy()
    del adat_wnn
    
    if add_3d_umap:
        tmp_mdata_3d = mudat.copy()
        logger.info("WNN - Generate 3DUMAP")
        mu.tl.umap(tmp_mdata_3d, n_components=3)
        mudat.obsm['X_umap_wnn_3d'] = tmp_mdata_3d.obsm['X_umap'].copy()
        del tmp_mdata_3d

    mudat.update()
    return mudat