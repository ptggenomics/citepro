from scipy.sparse import issparse
from muon import MuData
from anndata import AnnData
import numpy as np
from statistics import median
#import scanpy as sc

from typing import Union

from .utils import find_feature

"""pp module contains utility codes for data preprocessing
"""



def arcsinh_transform(adata: AnnData, densify = True, to_layers = True, cofactor=1, noise_mean = 0.5, noise_sd = 0.5) -> AnnData:
    """Implements the jittered hyperbolic arcsin transformation

    Parameters
    ----------
    adata : AnnData
        input anndata object, typically the prot modal of a mudata object mdata['prot']
    densify : bool, optional
        whether to convert to dense matrix, [default: True]
    to_layers : bool, optional
        whether to insert transformed value to .layer slot. data will be in the .X slot if this parameter is set to False [default: True]
    cofactor : int, optional
        denominator after asinh transformation, [default: 1]
    noise_mean : float, optional
        mean of gaussian distribution for jittering, [default: 0.5]
    noise_sd : float, optional
        standard deviation of gaussian distribution for jittering [default: 0.5]

    Returns
    -------
    AnnData
        the result anndata object, with transformed values in .X or .layers['arcsinh']
    """    
    
    x = adata.X
    if densify and issparse(x):
        x = np.array(x.todense())

    if issparse(x):
        ni = np.abs(np.random.normal(noise_mean, noise_sd, len(x.data)))
        transformed = x.copy()
        transformed.data = np.arcsinh( np.add(transformed.data, ni) / cofactor ) 

    else:
        ni = np.abs(np.random.normal(noise_mean, noise_sd, x.shape))
        transformed = np.arcsinh( np.add(x, ni) / cofactor )

    if to_layers:
        if adata.layers is None:
            adata.layers = dict()
        adata.layers['arcsinh'] = transformed
    else:
        adata.X = transformed

    return adata

def isotype_norm(mdata: MuData, isotype_dict, to_layers=False, inplace=True, default_isotype='66360-1') -> Union[np.ndarray, None]:
    """Implements the isotype normalization similar to visium

    Returns
    -------
    Union[np.ndarray, None]
        returns the isotype normalized counts, or None if inplace is True
    """    
    count_mat = np.asarray(mdata['prot'].X.todense())
    result_arr = np.zeros(count_mat.shape)
    for ab_idx, ab_count in enumerate(count_mat.T):
        ab_count = ab_count.reshape(-1)
        med_original = median(ab_count)
        
        isotype_id = isotype_dict.get(mdata['prot'].var['gene_ids'][ab_idx], None)
        if isotype_id:
            isotype_count = find_feature(mdata['prot'], id = isotype_id[:7]) 
            # [:7] is to strip off the dash number form the isotype name
            # the 'Notavailable' of isotype_name is handled in the following if block
            if isotype_count is None:
                isotype_count = find_feature(mdata['prot'], id = default_isotype)


        isotype_count = np.asarray(isotype_count.todense()).reshape(-1)

        divided = np.divide(ab_count +1 , isotype_count +1)
        scale_factor = med_original / median(divided)
        result_arr[:, ab_idx] = np.arcsinh(divided * scale_factor)
    
    if inplace:
        mdata['prot'].X = result_arr

        if to_layers:
            if mdata['prot'].layers is None:
                mdata['prot'].layers = dict()
            mdata['prot'].layers['isonorm'] = result_arr
        else:
            mdata['prot'].X = result_arr
        return
        
    else:
        return result_arr


def gen_adata_celltypist(mdata_raw: MuData, ct_model='Immune_All_Low.pkl', target_cpm=1e4, inplace = True) -> Union[AnnData,None]:
    """_summary_

    Parameters
    ----------
    mdata_raw : MuData
        _description_
    ct_model : str, optional
        _description_, by default 'Immune_All_Low.pkl'
    target_cpm : _type_, optional
        _description_, by default 1e4
    inplace : bool, optional
        _description_, by default True

    Returns
    -------
    Union[AnnData,None]
        _description_
    """    

    import celltypist as ct
    import scanpy as sc
    
    adata = mdata_raw['rna'].copy()
    adata.var.index = [gn.replace("rna:","") for gn in adata.var.index]
    sc.pp.normalize_total(adata, target_cpm)
    sc.pp.log1p(adata)
    ct.models.download_if_required()
    ct_model = ct.models.Model.load(model = ct_model)
    predictions = ct.annotate(adata, model = ct_model, majority_voting = True)
    adata_result = predictions.to_adata(prefix='ct_')
    adata_result.obs.rename(columns={'ct_predicted_labels':'celltype_ct_labels', 
                                     'ct_majority_voting':'celltype_ct_majvote',
                                     'ct_conf_score': 'celltype_ct_confscore'}, inplace = True)

    if inplace:
        mdata_raw.obs['celltype_ct_labels'] = adata_result.obs['celltype_ct_labels']
        mdata_raw.obs['celltype_ct_majvote'] = adata_result.obs['celltype_ct_majvote']
        mdata_raw.obs['celltype_ct_confscore'] = adata_result.obs['celltype_ct_confscore']
        return
    else:
        return adata_result

def calc_qc_prot_var(adata_prot: AnnData):
    """Calculate the qc metrics for protein features. The adata_prot.var metadata will be updated directly.

    Parameters
    ----------
    adata_prot : AnnData
        input anndata object or the prot modal of a mudata object mudat['prot'].
    """    
    prt_sum = np.asarray(adata_prot.X.sum(axis =0),dtype = np.int64).flatten()
    prt_percent = ['{:.5f}'.format(i) for i in np.divide(prt_sum, np.sum(prt_sum))*100]
    prt_count = adata_prot.X.todense()
    prt_medians = np.asarray(np.median(prt_count, axis=0),dtype = np.int64).flatten()
    prt_75th = np.asarray([np.quantile(np.asarray(prt_count[:,col]).reshape(-1), q = 0.75) for col in range(prt_count.shape[1])], dtype = np.int64)
    prt_95th = np.asarray([np.quantile(np.asarray(prt_count[:,col]).reshape(-1), q = 0.95) for col in range(prt_count.shape[1])], dtype = np.int64)
    adata_prot.var['sum'] = prt_sum
    adata_prot.var['percent'] = prt_percent
    adata_prot.var['median'] = prt_medians
    adata_prot.var['75th'] = prt_75th
    adata_prot.var['95th'] = prt_95th

def _summarize_by_cells(adata: AnnData, np_func = np.sum):
    return np.squeeze(np.asarray(np_func(adata.X, axis = 1), dtype=np.int64))


def calc_cell_qc(mudat: MuData):
    """_summary_

    Parameters
    ----------
    mudat : MuData
        _description_
    """    
    
    mt_idx = [idx for idx, vn in enumerate(mudat['rna'].var_names) if 'MT-' in vn]
    isotype_idx = [idx for idx, vn in enumerate(mudat['prot'].var_names) if 'IsoCtl' in vn]
    mudat.obs['totalumi_rna'] = _summarize_by_cells(mudat['rna'])
    mudat.obs['totalumi_mtgene'] = _summarize_by_cells(mudat['rna'][:,mt_idx])
    mudat.obs['pct_mt_umi'] = mudat.obs['totalumi_mtgene']/mudat.obs['totalumi_rna']

    mudat.obs['totalumi_prot'] = _summarize_by_cells(mudat['prot'])
    mudat.obs['totalumi_isotype'] = _summarize_by_cells(mudat['prot'][:,isotype_idx])
    mudat.obs['mean_isotype'] = _summarize_by_cells(mudat['prot'][:,isotype_idx], np_func= np.mean)