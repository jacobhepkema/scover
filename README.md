# scover 

[![Documentation Status](https://readthedocs.org/projects/scover/badge/?version=latest)](https://scover.readthedocs.io/en/latest/?badge=latest)

This GitHub repository contains the neural network tool **scover** associated with a Biorxiv pre-print: https://www.biorxiv.org/content/10.1101/2020.11.26.400218v2.

<img src="https://github.com/jacobhepkema/scover/raw/master/scover_logo.png" width=300 align=right>

__Q__: What is this? 

__A__: __scover__ is a convolutional neural network (CNN) for *de novo* inference of *cis*-regulatory motifs from single-cell data. 
It finds weights for these motifs across pseudo-bulks and also reports the 'influence' of each motif. The network is written in 
[pytorch](https://pytorch.org/).

## Workflow

I have attached a few example jupyter notebooks for in the `example_notebooks` directory: 1) [data generation using scRNA-seq data](https://github.com/jacobhepkema/scovernew/blob/master/example_notebooks/01_pool_dataset.ipynb), 2) [extracting promoter sequences](https://github.com/jacobhepkema/scovernew/blob/master/example_notebooks/02_get_sequences.ipynb), 3) [running scover](https://github.com/jacobhepkema/scovernew/blob/master/example_notebooks/03_run_scover.ipynb).

## Plotting

For plotting, see [the figures github](https://github.com/jacobhepkema/scoverplots).

## License
MIT
