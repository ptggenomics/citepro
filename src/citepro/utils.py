
import numpy as np
from anndata import AnnData
from scipy.sparse import issparse

from typing import Optional
from logging import getLogger

logger = getLogger('citepro')


def find_feature(adata: AnnData, id:str = None, name:str = None) -> Optional[np.ndarray]:
    """return the counts of the specified feature in the .X slot of supplied anndata object, or None if such id or name is not found.

    Parameters
    ----------
    adata : AnnData
        input anndata object
    id : str, optional
        id to search for. Typically is 'rna:CD8A' or 'prot:CD8A.65146.1'. must sepcify either id or name. [default: None]
    name : str, optional
        name to search for. Typically 'ENSGxxxxx' or proteintech id '12345-1'. must sepcify either id or name. [default: None]

    Returns
    -------
    Union[np.ndarray, None]
        return the counts of the specified feature, or None if not found

    Raises
    ------
    AttributeError
        raised if both id and name are specified.
    """    
    if not bool(id) ^ bool(name):
        raise AttributeError("Specify either id or name")

    if id:
        idx = np.where(adata.var_names == id)
    else:
        idx = np.where(adata.var['gene_ids'] == name)


    if len(idx[0])>0:
        ridx = idx[0][0]
        res = adata.X[:, ridx]
        if issparse(res):
            res = np.array(res.todense()).reshape(-1)
        return res
    else:
        return None
    
def find_encoding_rna(adata, prot:str) -> Optional[str]:
    if 'map_rna' not in adata.var_keys():
        logger.warning('No map_rna column found in adata.var')
        return None

    if prot not in adata.var.index:
        logger.warning(f'Protein {prot} not found in provided anndata')
        return None

    rna = adata.var.loc[prot,'map_rna']
    if rna == 'n.n':
        #logger.warning(f'No encoding RNA found for {prot}')
        return None

    if f"rna:{rna}" not in adata.var.index:
        logger.warning(f'Encoding RNA {rna} not found in provided anndata')
        return None

    return f"rna:{rna}"
