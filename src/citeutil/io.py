from muon import read_10x_mtx, MuData

def read_10x_mtx_filter(path, allow_file:str=None, block_file:str=None, *args, **kwargs) -> MuData:
    if allow_file is not None and block_file is not None:
        raise AttributeError("Cannot use bc_allow and bc_block together")

    mudat = read_10x_mtx(path, *args, **kwargs)

    if allow_file:
        with open(allow_file, 'r') as f:
            bc_allow = set([l.split(',')[0] for l in f.read().splitlines()])
        bc_keep = [bc in bc_allow for bc in mudat.obs_names]
        print(f'Supplied {len(bc_allow)} allowed barcodes, kept {sum(bc_keep)} / {mudat.shape[0]} barcodes')
        mudat = mudat[bc_keep, :]

    if block_file:
        with open(block_file, 'r') as f:
            bc_block = set([l.split(',')[0] for l in f.read().splitlines()])
        bc_keep = [bc not in bc_block for bc in mudat.obs_names]
        print(f'Supplied {len(bc_block)} blocked barcodes, kept {sum(bc_keep)} / {mudat.shape[0]} barcodes')
        mudat = mudat[bc_keep, :]

    mudat.var_names_make_unique()
    mudat.raw = mudat.copy()
    return mudat