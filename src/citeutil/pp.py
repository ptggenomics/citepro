from scipy.sparse import isspmatrix, issparse
from muon import MuData
import muon as mu
import anndata as ad
import numpy as np
from statistics import median

def arcsinh_transform(adata: ad.AnnData, densify = True, to_layers = True, cofactor=1, noise_mean = 0.5, noise_sd = 0.5):
    '''
    implements the hyperbolic arcsin transformation
    '''
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

def isotype_norm(mdata: mu.MuData, isotype_dict, to_layers=False, inplace=True, default_isotype='66360-1'):
    count_mat = np.asarray(mdata['prot'].X.todense())
    result_arr = np.zeros(count_mat.shape)
    for ab_idx, ab_count in enumerate(count_mat.T):
        ab_count = ab_count.reshape(-1)
        med_original = median(ab_count)
        
        isotype_name = isotype_dict.get(mdata['prot'].var['gene_ids'][ab_idx], None)
        if isotype_name:
            isotype_count = find_feature_id(mdata, isotype_name[:7]) 
            # [:7] is to strip off the dash number form the isotype name
            # the 'Notavailable' of isotype_name is handled in the following if block
            if isotype_count is None:
                isotype_count = find_feature_id(mdata, default_isotype)


        isotype_count = np.asarray(isotype_count.todense()).reshape(-1)

        divided = np.divide(ab_count +1 , isotype_count +1)
        scale_factor = med_original / median(divided)
        #print(f"sf: {scale_factor}")
        #np.log1p(np.divide(ab_count+1 , isotype_count +1))
        result_arr[:, ab_idx] = np.arcsinh(divided * scale_factor)
    
    if inplace:
        mdata['prot'].X = result_arr

        if to_layers:
            if mdata['prot'].layers is None:
                mdata['prot'].layers = dict()
            mdata['prot'].layers['isonorm'] = result_arr
        else:
            mdata['prot'].X = result_arr
        
    else:
        return result_arr
    
def gen_isotype_dict(annd:ad.AnnData):
    ## select only proteins
    prot_var = annd.var[annd.var['feature_types'] == 'Antibody Capture'][['gene_ids', 'IsoCtrl']] 

    ## map isotype names to isotype ids, like '66360-1': 'prot:Mu IgG1 IsoCtl.66360.1'
    isotype_names = {ab_id: ab_name for ab_name, (ab_id, _) in prot_var.iterrows() if 'IsoCtl' in ab_name} 

    ## create dictionary that will return a tuple of isotype id and isotype name for a given protein id
    isotype_dict = {row[0]: (row[1][:7], isotype_names.get(row[1][:7], 'NotAvail')) for _, row in prot_var.iterrows()} 

    return isotype_dict