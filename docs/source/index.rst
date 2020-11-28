scover
======

Welcome to the **scover** documentation. This page contains guides on data generation for **scover**, running **scover**, and analysis of **scover** output. 

**scover** is a convolutional neural network (CNN) for *de novo* inference of cis-regulatory motifs from single-cell data. 
It finds weights for these motifs across pseudo-bulks and also reports the 'impact' of each motif. 
The network is written in pytorch, with the downstream analyses written in R. 
Running the network and running the downstream analysis is implemented in a Nextflow pipeline. 
Furthermore, motifs are aligned with Tomtom from the MEME suite[1]. 
scover requires that cells are annotated for cell type (or other category). 
For scRNA-seq data TSS annotation is required, as the promoter sequences are obtained directly from genomic sequence relative to the TSS.


.. toctree::
   :maxdepth: 2
   :caption: Dataset creation:

   How to generate scover input using scRNA-seq data <how_to_prepare_scover_input_using_scRNA_seq.rst>
   How to generate scover input using scATAC-seq data <how_to_prepare_scover_input_using_scATAC_seq.rst>
   How to generate scover input using scRNA-seq data and create_dataset.R <how_to_prepare_scover_input_using_create_dataset.rst>
   How to extract promoter sequences <how_to_extract_promoter_sequences.rst>

.. toctree::
   :maxdepth: 2
   :caption: Running scover:
   
   How to install scover prerequisites <how_to_install_scover_prerequisites.rst>
   How to run scover <how_to_run_scover.rst>  
   Profiles <profiles.rst>
 
.. toctree::
   :maxdepth: 2
   :caption: Analysing scover output:

   How to analyse scover output <how_to_analyse_scover_output.rst>  



Method overview
###############

**scover** tries to predict expression or accessibility values across multiple samples from either promoter sequences or accessible region sequences. 
It does so by simultaneously inferring the features of the sequences that are predictive of the data, and weighing these features across the different samples.

One downside of scRNA-seq and scATAC-seq datasets is that they are very sparse. This makes it very hard to predict the data. To counter this, we 'pool cells' of the same cell type annotation together, and sum their expression / accessibility scores. This way, the sparsity goes down depending on how many cells you 'pool' together. One caveat of this is that if you pool 100 cells together in each pool, cell types with fewer than 100 cells will be omitted. 

An example of how sparsity goes down with pool size can be seen here:

.. image:: _static/sparsity.png
   :width: 400
   :alt: Sparsity example

After this data preprocessing, the network `can be run <https://scover.readthedocs.io/en/latest/how_to_run_scover.html>`_. The network is a shallow convolutional neural network, with one convolutional layer, one global maximum pooling layer, and one fully connected layer (with the amount of pools as the amount of output channels). The network will start with initially random convolutional kernels (think of these as "motif detectors") and fully connected layer weights (think of these as the "impact scores of motifs in the pools"). The training step is the most computationally expensive and can take quite some hours depending on your dataset size. This is why I advice to run **scover** on GPUs. 

After the training step has been completed, **scover** will align the motifs back to a motif database using Tomtom [1], build a motif alignment graph, and identify reproducible motif clusters. Though **scover** will generate some output plots, for further downstream analysis I recommend the guide for analysing **scover** output `here <https://scover.readthedocs.io/en/latest/how_to_analyse_scover_output.html>`_


References
##########

[1]  Gupta, S., Stamatoyannopoulos, J. A., Bailey, T. L., & Noble, W. S. (2007). Quantifying similarity between motifs. Genome biology, 8(2), R24.
