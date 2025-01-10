# CITE-seq analysis python toolbox

<img src="https://www.ptglab.com/img/logos/PTG_Genomics_logo.png" width = 300/>

A collection of code that we think is useful for analyzing Multiomic single cell RNA-seq data

# Prerequisite - python package manager

We recommend using (uv)[https://docs.astral.sh/uv/] for dependency management. Please follow the instruction on the (installation)[https://docs.astral.sh/uv/getting-started/installation/] for installing uv.

After installing uv, create a virtual environment is recommended.
```sh
uv venv --python 3.11
```

And then to activate the virtual environment, run the following code on Linux
```sh
source activate .venv/bin/activate
```
or on Windows

```
.venv\Script\activate
```

----
# Install/Upgrade

To install this package, use the following pip command
```
uv pip install -u git+https://github.com/ptggenomics/citepro.git
```


# Example Usage
This set of code are can work in jupyter/colab environment or being incorporated into your processing script

```python

import citepro 

## read filtered barcode matrix into a mudata object
mdata = citepro.recipe.create_mudata("/path/to/sample_filtered_feature_barcode_matrix.h5", allow_file="/path/to/blocked_barcodes.txt")

## generate anndata object from the mudata
adata = citepro.recipe.mu_to_ann(mdata)

```

Please visit the (citepro-notebook)[http://github.com/ptggenomics/citepro_notebook] repository for more tutorial code


----


This repository is under active development
