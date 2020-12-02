# scover 

[![Documentation Status](https://readthedocs.org/projects/scover/badge/?version=latest)](https://scover.readthedocs.io/en/latest/?badge=latest)

This GitHub repository contains the files needed and instructions to run the Nextflow pipeline of **scover**. The primary developer of 
scover is Jacob Hepkema who is currently a PhD student at the [Wellcome Sanger Insitute](https://www.sanger.ac.uk).

<img src="https://github.com/jacobhepkema/scover/raw/master/scover_logo.png" width=300 align=right>

__Q__: What is this? 

__A__: __scover__ is a convolutional neural network (CNN) for *de novo* inference of *cis*-regulatory motifs from single-cell data. 
It finds weights for these motifs across pseudo-bulks and also reports the 'influence' of each motif. The network is written in 
[pytorch](https://pytorch.org/), with the downstream analyses written in R (using [ggplot2](https://ggplot2.tidyverse.org/) for plotting). 
Running the network and running the downstream analysis is implemented in a Nextflow pipeline. Furthermore, motifs are aligned with 
[Tomtom](http://meme-suite.org/tools/tomtom) from the MEME suite[1]. __scover__ requires that cells are annotated for cell type 
(or other category). For scRNA-seq data TSS annotation is required, as the promoter sequences are obtained directly from 
genomic sequence relative to the TSS. 

The Biorxiv pre-print is now available at [this page](https://www.biorxiv.org/content/10.1101/2020.11.26.400218v1).

__Q__: All right, so when would I use this?

__A__: If you have a clustered scRNA-seq or scATAC-seq dataset, and you want to extract regulatory information, you have come to the right place. 
__scover__ will work as a 'hypothesis generator' to help you identify the most important regulatory elements for your dataset. 
In general, it is hard to pinpoint exactly which transcription factors are the importance since many will bind to the same motif family. 
However, __scover__ will help you by providing a list of transcription factors that it deems to be the most interesting to follow up on. 
Furthermore, __scover__ gives information about the relative impact of motifs in the different cell types in your dataset.
For scRNA-seq, __scover__ will find motifs using the proximal promoter sequences, and for scATAC-seq, __scover__ will find motifs using
the peak regions. In both cases, __scover__ will also find motif influence scores across the cell types in your dataset. 

__Q__: How to install/run __scover__?

__A__: To install __scover__, follow these steps:

1. Install the prerequisites [by following the guide here](https://scover.readthedocs.io/en/latest/how_to_install_scover_prerequisites.html). 
2. Choose the right profile to run with (or create a new profile) [by following the guide here](https://scover.readthedocs.io/en/latest/profiles.html).
3. Run scover by [following the steps on this page](https://scover.readthedocs.io/en/latest/how_to_run_scover.html). If you want to, try to run using an example dataset by [following the steps on this page](https://scover.readthedocs.io/en/latest/how_to_run_scover.html#run-an-example-dataset). I advise to run __scover__ using GPUs, as the run times can increase significantly (~5x in some of our benchmarks) when using CPUs.

__Q__: Can I run this on Windows?

__A__: I have not tried yet, but because __scover__ uses Nextflow which is made for UNIX environments, it might get somewhat tricky. My advice is to first install a virtual machine if you are using a Windows computer. 

## Workflow

After you have installed all the prerequisites (see above), there are three main steps in using __scover__:
1. Preparing your dataset - see [here for scRNA-seq](https://scover.readthedocs.io/en/latest/how_to_prepare_scover_input_using_scRNA_seq.html) and [here for scATAC-seq](https://scover.readthedocs.io/en/latest/how_to_prepare_scover_input_using_scATAC_seq.html). For a very quick scRNA-seq data generation script, see [this link](https://scover.readthedocs.io/en/latest/how_to_prepare_scover_input_using_create_dataset.html)
3. Training __scover__ - see the [guide on this page](https://scover.readthedocs.io/en/latest/how_to_run_scover.html). I advise to run __scover__ using GPUs, as the run times can increase significantly (~5x in some of our benchmarks) when using CPUs.
4. Analysing the output - see [this link](https://scover.readthedocs.io/en/latest/how_to_analyse_scover_output.html)

## Questions and errors
If you have any questions, or want to report an error, please use our [github issues page](https://github.com/jacobhepkema/scover/issues)

## References
[1] Gupta, S., Stamatoyannopoulos, J. A., Bailey, T. L., & Noble, W. S. (2007). Quantifying similarity between motifs. Genome biology, 8(2), R24.

[2] Weirauch, M. T., Yang, A., Albu, M., Cote, A. G., Montenegro-Montero, A., Drewe, P., ... & Zheng, H. (2014). Determination and inference of eukaryotic transcription factor sequence specificity. Cell, 158(6), 1431-1443.

## License
MIT
