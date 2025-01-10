# CITE-seq analysis python toolbox

<img src="https://www.ptglab.com/img/logos/PTG_Genomics_logo.png" width = 300/>

A collection of code that we think is useful for analyzing Multiomic single cell RNA-seq data

----
# Install/Upgrade

To install this package, use the following pip command
```
pip install git+https://https://github.com/ptggenomics/citepro.git
```

To install this package, use the following pip command
```
pip install -u git+https://https://github.com/ptggenomics/citepro.git
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


This file is under active development
