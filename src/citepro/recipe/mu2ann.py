import muon as mu
from muon import MuData
from anndata import AnnData
import anndata as ad
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix, isspmatrix, hstack


def check_mu_req(mdat: mu.MuData):
    modalities = ['prot','rna']
    #obs_checks = ['sample_id']
    obs_checks = []

    ##  ====== checkup for combined mudata ======

    # if not 'sample_desc' in mdat.uns:
    #     raise KeyError(f'No sample description found')

    # if not isinstance(mdat.uns['sample_desc'],pd.core.frame.DataFrame):
    #     raise TypeError(f'Sample description is not a pandas dataframe')
    
    for obscol in obs_checks:
        if obscol not in mdat.obs_keys():
            raise KeyError(f'No {obscol} found for combined mudata')
    
    if len([1 for k in mdat.obs_keys() if k.startswith('celltype_') ]) < 1:
            raise KeyError('No celltype found')
  
    ## ====== checkup common for both side ======
    for md_key in modalities:
        
        if md_key not in mdat.mod:
            raise KeyError(f'No {md_key} modal found in input mudata')

        if mdat[md_key].layers is None:
            raise KeyError(f'No layers in {md_key} modal')
        
        if 'X_umap' not in mdat[md_key].obsm_keys():
            raise KeyError(f'No umap found for {md_key} modal')
        
        #for obscol in obs_checks:
        #    if obscol not in mdat[md_key].obs_keys():
        #        raise KeyError(f'No {obscol} found for {md_key} modal')
            
        
    ## ====== Protein side checkup =========
    if 'arcsinh' not in mdat['prot'].layers :
        raise KeyError('Layer arcsinh not found on protein side')

    
    ## ====== RNA side checkup ======
    if 'norm_log1p' not in mdat['rna'].layers:
        raise KeyError('Layer norm_log1p not found on rna side')


def _create_csr_matrix(X_mat) -> csr_matrix:
    try: 
        import cupy
        import cupyx
        if isinstance(X_mat, cupy.ndarray):
            return csr_matrix(cupy.asnumpy(X_mat), dtype= np.float32)
        elif cupyx.scipy.sparse.issparse(X_mat):
            return X_mat.get()

    except ImportError:
            return csr_matrix(X_mat, dtype= np.float32)

def collapse_X(mdat: mu.MuData):
    """Generate new sparse matrix from the transformed/normalized data"""
    #if isspmatrix(mdat['prot'].layers['arcsinh']):
    #    spm_prot = mdat['prot'].layers['arcsinh']
    #else:
    #    spm_prot = csr_matrix(mdat['prot'].layers['arcsinh'],dtype= np.float32)
    
    #if isspmatrix(mdat['rna'].layers['norm_log1p']):
    #    spm_rna = mdat['rna'].layers['norm_log1p']
    #else:
    #    spm_rna = csr_matrix(mdat['rna'].layers['norm_log1p'])

    if isspmatrix(mdat['prot'].X):
        spm_prot = mdat['prot'].X
    else:
        spm_prot = _create_csr_matrix(mdat['prot'].X)
    
    if isspmatrix(mdat['rna'].X):
        spm_rna = mdat['rna'].X
    else:
        spm_rna = _create_csr_matrix(mdat['rna'].X)

    return hstack([spm_prot,spm_rna], format='csr')


def collapse_obsm(mudat: mu.MuData):
    new_obsm = {'X_umap_wnn': mudat.obsm['X_umap_wnn'],
               'X_umap_rna': mudat['rna'].obsm['X_umap'],
               'X_umap_prot': mudat['prot'].obsm['X_umap']}
    
    if 'X_umap_wnn_3d' in mudat.obsm_keys():
        new_obsm['X_umap_wnn_3d'] = mudat.obsm['X_umap_wnn_3d']
    
    for mdkey in mudat.mod.keys():
        if 'X_umap_3d' in mudat[mdkey].obsm_keys():
            new_obsm[f'X_umap_{mdkey}_3d'] = mudat[mdkey].obsm['X_umap_3d']
    
    
    return new_obsm

def collapse_var(mudat: mu.MuData):

    var_prot = mudat['prot'].var
    #print(var_prot)
    var_rna = mudat['rna'].var
    new_var = pd.concat([var_prot,var_rna],axis = 0, join = 'outer')

    return new_var

def cleanup_obs(adat: AnnData) -> pd.DataFrame:
    """cleanup obs for combined anndata

    Parameters
    ----------
    adat : AnnData
        anndata to be cleaned up

    Returns
    -------
    pd.DataFrame
        cleaned up obs
    """    
    ## remove mod_weight columns
    mod_weight_cols = [col for col in adat.obs.columns if 'mod_weight' in col]
    new_obs = adat.obs.drop(columns = mod_weight_cols, inplace = False)

    new_obs.columns = [col.replace('rna:','').replace('prot:','') for col in new_obs.columns]
    
    return new_obs

def collapse_X_raw(mdat: mu.MuData):
    """Generate new sparse matrix from the transformed/normalized data"""
    if isspmatrix(mdat['prot'].raw.X):
        spm_prot = mdat['prot'].raw.X
    else:
        spm_prot = csr_matrix(mdat['prot'].raw.X,dtype= np.float32)
    
    if isspmatrix(mdat['rna'].raw.X):
        spm_rna = mdat['rna'].raw.X
    else:
        spm_rna = csr_matrix(mdat['rna'].raw.X)

    return hstack([spm_prot,spm_rna], format='csr')

def collapse_raw_mu(mdat: mu.MuData) -> AnnData:
    ## generate X from the transformed/normalized data
    
    new_ann = AnnData(X = collapse_X_raw(mdat), 
                         var = collapse_var(mdat))
    
    return new_ann
# def collapse_mu(mdat: mu.MuData, prot_var:pd.DataFrame):
#     ## generate X from the transformed/normalized data
    
#     new_ann = ad.AnnData(X = collapse_X(mdat), 
#                          obsm = collapse_obsm(mdat),obs = mdat.obs,var = collapse_var(mdat))
    
#     return new_ann

def collapse_mu(mdat: mu.MuData) -> AnnData:
    ## generate X from the transformed/normalized data
    
    new_ann = ad.AnnData(X = collapse_X(mdat), 
                         obsm = collapse_obsm(mdat),
                         obs = mdat.obs,
                         var = collapse_var(mdat))
    
    new_ann.obs = cleanup_obs(new_ann)
    
    return new_ann
    

def mu_to_ann(mudat: MuData, include_raw= False) -> AnnData:
    """
    create an AnnData object from a MuData object. 

    Parameters
    ----------
    mudat : MuData
        input mudata object
    include_raw : bool, optional
        whether to include the raw slot

    Returns
    -------
    AnnData
        anndata object with combined data
    """    

    check_mu_req(mudat)
    
    # collapsing mudata to anndata
    comb_ann = collapse_mu(mudat)
    if include_raw:
        comb_ann.raw = collapse_raw_mu(mudat)

    return comb_ann
