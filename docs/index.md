# Introduction
Welcome!

!!! tip "Try citepro on Colab" 

    Google Colab is a ready to use, setup-free notkbook environment that you can launch directly. Please sign up for an google account, and click the launch button on the right -> \( <a target="_blank" href="https://colab.research.google.com/github/ptggenomics/citepro_notebook/blob/28d28b9cbcac2ceb2c5cf09ea750af3335e22dfa/Basic_usage_GUI.ipynb"> :simple-googlecolab:  </a> \) to try citepro with GUI directly.


Citepro is a collection of wrapper codes that is aimed to provide streamlined data processing experience to convert Multiomic RNA+Protein cellranger count matrix to h5ad, a data container format that is compatible with downstream analysis code written in python or R  
 
Citepro is built on the popular and well known python packages in the [scverse](https://scverse.org/) ecosystem, such as the following:

 * [scanpy](https://scanpy.readthedocs.io/en/stable/), the ***de facto*** analysis tool suite for single cell data analysis.

 * [anndata](https://anndata.readthedocs.io/en/stable/), the annotated data matrices container format.

 * [muon](https://muon.readthedocs.io/en/latest/), multiomic data container format the is built on top of anndata.

 * [celltypist](https://www.celltypist.org/), an automated cell type annotating tool that use RNA signature to predict the celltypes \(such as T or B cells\). 

The purpose of citepro is to provide a suggestion of using the publicly available libraries in a research-only manner, the user is ultimately responsible of determine the quality of data generated. 


# Tutorials
Go to the [Tutorial](/tutorials/) page for interactive notebooks that can be launch on [Google Colab](https://research.google.com/colaboratory/faq.html)

# Reference
Go to the [Reference](/reference/) page for a detailed list of options to customize the setting for the functions, such as adding block list or allow list accompanied with the count matrix 

!!! note

    The citepro is under active development, and codes may be changing frequently.