# scover 

This GitHub repository contains the neural network tool **scover** associated with a publication: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03021-9.

<img src="https://github.com/jacobhepkema/scover/raw/master/scover_logo.png" width=300 align=right>

__Q__: What is this? 

__A__: __scover__ is a convolutional neural network (CNN) for *de novo* inference of *cis*-regulatory motifs from single-cell data. 
It finds weights for these motifs across pseudo-bulks and also reports the 'influence' of each motif. The network is written in 
[pytorch](https://pytorch.org/).

## Workflow

I have attached a few example jupyter notebooks for in the `example_notebooks` directory: 1) [data generation using scRNA-seq data](https://github.com/jacobhepkema/scover/blob/master/example_notebooks/01_pool_dataset.ipynb), 2) [extracting promoter sequences](https://github.com/jacobhepkema/scover/blob/master/example_notebooks/02_get_sequences.ipynb), 3) [running scover](https://github.com/jacobhepkema/scover/blob/master/example_notebooks/03_run_scover.ipynb). Note that for the first notebook, the dataset is gzipped (otherwise it would exceed github file size limits), running `gunzip Marrow.h5ad.gz` will give the file necessary for that notebook.

For the pooling step, we have [added an additional notebook with the option to run SEACells instead of our pooling algorithm](https://github.com/jacobhepkema/scover/blob/master/example_notebooks/01b_pool_dataset_SEACells.ipynb), but we have _not_ tried to train the model using data generated with SEACells yet.

## Plotting

For plotting, see [the figures github](https://github.com/jacobhepkema/scoverplots).

## Data availability

For the datasets used in the manuscript, please see our [Zenodo](https://doi.org/10.5281/zenodo.8060659).

## License
MIT
