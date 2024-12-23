# CITE-seq analysis python toolbox

<img src="https://www.ptglab.com/img/logos/PTG_Genomics_logo.png" width = 300/>

A collection of code that we think is useful for analyzing Multiomic single cell RNA-seq data

----
# Install/Upgrade

To install this package, use the following pip command
```
pip install git+https://https://github.com/ptggenomics/citeutil.git
```

To install this package, use the following pip command
```
pip install -u git+https://https://github.com/ptggenomics/citeutil.git
```

# Example Usage
This set of code are can work in jupyter/colab environment or being incorporated into your processing script

```python

import citeutil

## read filtered barcode matrix into a mudata object
mudat = read_10x_mtx_filter("/path/to/sample_filtered_feature_barcode_matrix.h5", bc_block="/path/to/blocked_barcodes.txt")

## perform arcsinh transformation on the protein modality
arcsinh_transformation(mudat['prot'], inplace=True, densify= True)

```

Please visit the notebook repository for more tutorial code


----


This file is under active development

KitD made this change at feat branch