import muon as mu
import numpy as np

def read_10x_mtx_filter(path, bc_allow=None, bc_block=None, *args, **kwargs) -> mu.MuData:
    mudat = mu.read_10x_mtx(path, *args, **kwargs)

    if bc_allow is not None and bc_block is not None:
        raise ValueError("Both bc_allow and bc_block cannot be used together")

    if bc_allow:
        bc_keep = [bc in bc_allow for bc in mudat.obs_names]
        print(f'Supplied {len(bc_allow)} allowed barcodes, kept {sum(bc_keep)} / {mudat.shape[0]} barcodes')
        mudat = mudat[bc_keep, :]

    if bc_block:
        bc_keep = [bc in bc_block for bc in mudat.obs_names]
        print(f'Supplied {len(bc_block)} blocked barcodes, removed {sum(bc_keep)} / {mudat.shape[0]} barcodes')
        mudat = mudat[np.logical_not(bc_keep), :]

    mudat.var_names_make_unique()
    mudat.raw = mudat.copy()
    return mudat