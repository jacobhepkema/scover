scanem
======

Welcome to the **scanem** documentation. This page contains guides on data generation for **scanem** and analysis of **scanem** output. For information on how to run **scanem**, `click here to go to the main repository <https://github.com/jacobhepkema/scanem>` 

**scanem** is a convolutional neural network (CNN) for *de novo* inference of cis-regulatory motifs from single-cell data. 
It finds weights for these motifs across pseudo-bulks and also reports the 'impact' of each motif. 
The network is written in pytorch, with the downstream analyses written in R. 
Running the network and running the downstream analysis is implemented in a Nextflow pipeline. 
Furthermore, motifs are aligned with Tomtom from the MEME suite[1]. 
scanem requires that cells are annotated for cell type (or other category). 
For scRNA-seq data TSS annotation is required, as the promoter sequences are obtained directly from genomic sequence relative to the TSS.



.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   How to run scanem <how_to_run_scanem.rst>  
   Profiles <profiles.rst>
   How to generate scanem input using scRNA-seq data <how_to_prepare_scanem_input_using_scRNA_seq.rst>
   How to generate scanem input using scATAC-seq data <how_to_prepare_scanem_input_using_scATAC_seq.rst>
   How to generate scanem input using scRNA-seq data and create_dataset.R <how_to_prepare_scanem_input_using_create_dataset.rst>
   How to analyse scanem output <how_to_analyse_scanem_output.rst>  


References
##########

[1]  Gupta, S., Stamatoyannopoulos, J. A., Bailey, T. L., & Noble, W. S. (2007). Quantifying similarity between motifs. Genome biology, 8(2), R24.
