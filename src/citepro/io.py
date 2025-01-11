from muon import read_10x_mtx, read_10x_h5, MuData
import pandas as pd

from pathlib import Path
from typing import Optional

import logging
logger = logging.getLogger('citepro')
logger.setLevel(logging.INFO)
"""the io module contains utility codes for reading data
"""


def _read_maprna(path_map_rna:str):
    logger.info('Reading cellranger feature reference file with map_rna column')
    if not Path(path_map_rna).exists():
        raise FileNotFoundError(f"File not found: {path_map_rna}")
    
    map_rna = pd.read_csv(path_map_rna)

    if 'map_rna' not in map_rna.columns:
        raise AttributeError('No "map_rna" column in the supplied feature reference file')
    
    if 'map_rna' not in map_rna.columns:
        raise AttributeError('No "map_rna" column in the supplied feature reference file')

    if 'IsoCtrl' in map_rna.columns:
        logger.info('IsoCtrl detected in the supplied reference file')
        map_rna = map_rna[['id','map_rna','ENS.GENE','IsoCtrl']]
    else:
        logger.info('No IsoCtrl detected in the supplied reference file')
        map_rna = map_rna[['id','map_rna','ENS.GENE']]
    
    return map_rna


def read_10x_filter(path_count, allow_file:str=None, block_file:str=None, path_map_rna:Optional[str] = None,  *args, **kwargs) -> MuData:
    if allow_file is not None and block_file is not None:
        raise AttributeError("Cannot use bc_allow and bc_block together")

    if not Path(path_count).exists():
        raise FileNotFoundError(f"File not found: {path_count}")
    
    input_is_h5 = path_count.endswith('.h5')

     
    if input_is_h5:
        logger.info('Reading h5 format of count matrix')
        #mdata = read_10x_h5_filter(path_count, allow_file=bc_allow, block_file=bc_block)
        mudat = read_10x_h5(path_count, *args, **kwargs)
    
    else:
        logger.info('Reading mtx format of count matrix')

        if path_map_rna:
            map_rna = _read_maprna(path_map_rna)
        else:
            raise AttributeError("Mtx format count matrix must be supplied with feature reference csv.")
        
        mudat = read_10x_mtx(path_count, *args, **kwargs)
        var_mod=mudat['prot'].var.reset_index().merge(map_rna, left_on = 'gene_ids',right_on = 'id',how = 'inner').set_index('index').copy()
        var_mod.drop(columns = ['id'],inplace = True)
        mudat['prot'].var = var_mod.copy()
        del var_mod
        #mdata = read_10x_mtx_filter(path_count, allow_file=bc_allow, block_file=bc_block, map_rna=map_rna)
    
    mudat.var_names_make_unique()
    
    if 'map_rna' not in mudat['prot'].var_keys():
        logger.warning('')

    if allow_file:
        with open(allow_file, 'r') as f:
            bc_allow = set([l.split(',')[0] for l in f.read().splitlines()])
        bc_keep = [bc in bc_allow for bc in mudat.obs_names]
        logger.info(f'Supplied {len(bc_allow)} allowed barcodes, kept {sum(bc_keep)} / {mudat.shape[0]} barcodes')
        mudat = mudat[bc_keep, :]

    if block_file:
        with open(block_file, 'r') as f:
            bc_block = set([l.split(',')[0] for l in f.read().splitlines()])
        bc_keep = [bc not in bc_block for bc in mudat.obs_names]
        logger.info(f'Supplied {len(bc_block)} blocked barcodes, kept {sum(bc_keep)} / {mudat.shape[0]} barcodes')
        mudat = mudat[bc_keep, :]

    mudat['rna'].var_names  = [f"{'rna'}:{idx}" for idx in mudat['rna'].var.index]
    
    mudat.raw = mudat.copy()
    mudat['prot'].raw = mudat['prot'].copy()
    mudat['rna'].raw = mudat['rna'].copy()
    
    mudat.update()

    return mudat