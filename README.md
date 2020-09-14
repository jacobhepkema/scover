# scanem 

This GitHub repository contains the files needed and instructions to run the Nextflow pipeline of **scanem**. The primary developer of scanem is Jacob Hepkema who is currently a PhD student at the [Wellcome Sanger Insitute](https://www.sanger.ac.uk).

<img src="https://github.com/jacobhepkema/scanem/raw/master/scanem_logo.png" width=300 align=right>

__Q__: What is this? 

__A__: __scanem__ is a convolutional neural network (CNN) for *de novo* inference of *cis*-regulatory motifs from single-cell data. It finds weights for these motifs across pseudo-bulks and also reports the 'impact' of each motif. The network is written in [pytorch](https://pytorch.org/), with the downstream analyses written in R (using [ggplot2](https://ggplot2.tidyverse.org/) for plotting). Running the network and running the downstream analysis is implemented in a Nextflow pipeline. Furthermore, motifs are aligned with [Tomtom](http://meme-suite.org/tools/tomtom) from the MEME suite[1]. __scanem__ requires that cells are annotated for cell type (or other category). For scRNA-seq data TSS annotation is required, as the promoter sequences are obtained directly from genomic sequence relative to the TSS. 

__Q__: All right, so when would I use this?

__A__: If you have a clustered scRNA-seq or scATAC-seq dataset, and you want to extract regulatory information, you have come to the right place. __scanem__ will work as a 'hypothesis generator' to help you identify the most important regulatory elements for your dataset. In general, it is hard to pinpoint exactly which transcription factors are the importance since many will bind to the same motif family. However, __scanem__ will help you by providing a list of transcription factors that it deems to be the most interesting to follow up on. Furthermore, __scanem__ gives information about the relative impact of motifs in the different cell types in your dataset.

__Q__: How to install/run __scanem__?

__A__: To install __scanem__, simply clone the repository to the directory where you want to run it. Before running, there is some data-preprocessing required (see [workflow](#workflow)). To run __scanem__, you will also need to have [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) and [Singularity](https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps) installed. The guide on how to run __scanem__ can be found [further down this page](#training-scanem). I advise to run __scanem__ using GPUs, as the run times can increase significantly (~5x in some of our benchmarks) when using CPUs.

__Q__: How can I run __scanem__ using GPUs?

__A__: Currently, the Singularity image for training the network only detects CPUs, so to use GPUs, you will need to create an anaconda environment with [pytorch](https://pytorch.org/) and other required packages. This is done very easily; if you have [installed anaconda or miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) you can `cd` into the `scanem` folder, which contains a `scanem_env.yml` file that specifies which files should be in the environment. To create the conda environment, simply type 

```{bash}
conda env create -f scanem_env.yml
```

which will create the `scanem` environment with all the required packages. Before running __scanem__, simply activate the environment with `conda activate scanem` so that when you [start training scanem](#training-scanem) with Nextflow, it will have the right packages. Then, run __scanem__ with `profile -local_gpu` if you want to run locally, or with `profile -lsf_gpu` if you want to use the Platform LSF scheduler. We have not yet tried to run __scanem__ using other schedulers.
(For running on CPU, the conda step is not needed, as it will use a Singularity image that has the right packages)

__Q__: Can I run this on Windows?

__A__: I have not tried yet, but because __scanem__ uses Nextflow which is made for UNIX environments, it might get somewhat tricky. My advice is to first install a virtual machine if you are using a Windows computer. 

---------------------------------------------------------------------------------------------------
## Workflow

There are three main steps in using __scanem__:
1. Preparing your dataset - see [here for scRNA-seq](https://scanem.readthedocs.io/en/latest/how_to_prepare_scanem_input_using_scRNA_seq.html) and [here for scATAC-seq](https://htmlpreview.github.io/?https://github.com/jacobhepkema/scanem/blob/master/guides/how_to_prepare_scanem_input_using_scATAC_seq.html). For a very quick scRNA-seq data generation script, see [this link](https://scanem.readthedocs.io/en/latest/how_to_prepare_scanem_input_using_create_dataset.html)
2. Training __scanem__ - see the guide further on this [page](#training-scanem). 
3. Analysing the output - see [this link](https://scanem.readthedocs.io/en/latest/how_to_analyse_scanem_output.html)

## Install/Dependencies

Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) by following the instructions on [their website](https://www.nextflow.io/docs/latest/getstarted.html#installation). 
Move the nextflow file to a directory accessible by your `$PATH` variable so that nextflow can be run by typing `nextflow` rather than the full path to the executable. One way of doing this is with the following line:
```
export PATH="path/to/nextflow:$PATH"
```
Alternatively, manually add `export PATH="path/to/nextflow:$PATH"` to your `~/.bashrc` or `~/.bash_profile` file, depending on which file is used as a source for the command line environment. 

To get __scanem__, clone this repository into your current working directory with the following command:
```
git clone https://github.com/jacobhepkema/scanem
```
Alternatively, download the repository directly and place in a folder of your choice. 
`cd` into the top directory (`/scanem`) to run __scanem__.


---------------------------------------------------------------------------------------------------
## Training scanem

Usage:
```
nextflow run [options] scanem.nf [arg...]
```

Example command:
```
nextflow run -profile lsf_gpu scanem.nf \
  --name test_run /
  --data data/mock_data.tsv \
  --celldata data/mock_data_colData.tsv \
  --tomtom resources/Mus_musculus.meme \
  --num_calibrations 30 \
  --val_factor 10 \
  --epochs 50 \
  --num_errtest 10 \
  --motif_amount 300 \
  --motif_length 12
```

__scanem__ requires a run name (`--name`), a pooled dataset (`--data`), and pool cell type annotations (`--celldata`). 
See the guides for [scRNA-seq](https://htmlpreview.github.io/?https://github.com/jacobhepkema/scanem/blob/master/guides/how_to_prepare_scanem_input_using_scRNA_seq.html) and [scATAC-seq](https://htmlpreview.github.io/?https://github.com/jacobhepkema/scanem/blob/master/guides/how_to_prepare_scanem_input_using_scATAC_seq.html) on how to create the pooled dataset. Another important argument is `--tomtom`, for specifying
the motif `.meme` database file to align found motifs to. I have included `Mus_musculus.meme` and `Homo_sapiens.meme`
from CIS-BP[2] in the `resources` directory. 

Before running __scanem__ you will probably want to identify the best `-profile` to run with. This will define the executor
used by __scanem__. [See below](#profile) for options and customisation. 

Minimal command:
```
nextflow run scanem.nf \
  --name test_run \
  --data data/mock_data.tsv \
  --celldata data/mock_data_colData.tsv \
  --tomtom resources/Mus_musculus.meme
```

## Important:

It might seem like it is stuck during training as the line with 
```
[62/3eb656] process > train_scanem       [  0%] 0 of 1
```
keeps saying `[  0%]`. This happens because the training stage only contains one
neural network step, and it will switch to `[100%]` once the training is complete. 
There is a way to find out how far into training the network is, but it is not entirely straightforward:
Open a second terminal window and `cd` into the `scanem/work` directory. From there, you 
can `cd` into the directory starting with the code shown in front of the specific 
task (in this case, `62/3eb656`, so the directory will be something like `scanem/work/62/3eb656aae5e84c420b7aa267dfeb57`). 
Then, by running `tail .command.log` you can get an idea of how far along the training is. 

---------------------------------------------------------------------------------------------------
## Options

### Profile

`-profile`
Choose a specified profile configuration from the `/conf` directory. These profiles are included in the 
`nextflow.config` file; if you are adding an option do not forget to add it to `nextflow.config`. Options currently include:

* `-profile local_gpu` will run locally, and will use a GPU if available. Note that you need to have the right packages installed (see "How can I run __scanem__ using GPUs?" at the start of this README). 
* `-profile local_cpu` will run locally, and it will also use Singularity for the first step in the workflow. This will only support CPUs. 
* `-profile local_nosingularity` will run locally and will not use Singularity images. 
* `-profile lsf_gpu` will run using the [Platform LSF](https://en.wikipedia.org/wiki/Platform_LSF) scheduler using GPUs. This might require some editing of the `conf/lsf_gpu.conf` file to be compatible with GPU queues on your LSF setup. See [this page](https://www.nextflow.io/docs/latest/executor.html) for more information on how to specify executors.
* `-profile lsf_cpu` will run using the [Platform LSF](https://en.wikipedia.org/wiki/Platform_LSF) scheduler using CPUs. This might require some editing of the `conf/lsf_cpu.conf` file to be compatible with your LSF setup. See [this page](https://www.nextflow.io/docs/latest/executor.html) for more information on how to specify executors.

#### Important: 

Since the Singularity environments can be quite large, it might be helpful to store them in a specific location. In the `.conf` file of your choice, I advice adding a line that includes the `cacheDir` such that the images are stored in that directory:

```
singularity {
  runOptions = '--no-home --cleanenv'
  enabled = true
  autoMounts = true

  cacheDir = "ADD/PATH/TO/CACHE/DIR"      <---  Add this
}
```

### Report

`-with-report name.html`
Generates a Nextflow report webpage with information on task run times, CPU usage, memory usage, and I/O. Note that this does _not_ include information on GPU usage.

<!-- Options include `local` (default option), `lsf`, and `docker`. The `lsf` profile  uses a [Singularity](https://sylabs.io/docs/) image (built from [this git repo](https://github.com/jacobhepkema/scanem-wip) using [Singularity Hub](https://singularity-hub.org/)). For more options, see the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html) or run `nextflow run -h`. The cache directory for the Singularity image can be set in `/conf/lsf.config` by setting `cacheDir = ...` -->

## Arguments

`--name`
The name of the current experiment. 
Default value: `experiment`

`--data`
The path to the dataset (relative to `SCANEM.nf`) to train the network. Add your data in the example tab-separated format into the `/data` folder. 
Default value: `/data/input_data.tsv`.

`--celldata`
The path to the dataset annotation (relative to `SCANEM.nf`) to annotate the cells. 
Add your data in the example format (see default file) into the `/data` folder. 
For more information, see the data generation guide in this repository. 
Default value: `/data/input_data.tsv`.

`--tomtom`
The path to the `.meme` format motif database to align found motifs to. 
Default value: `resources/Mus_musculus.meme`.

`--num_calibrations`
The number of randomly intialized calibrations for hyperparameter optimization. 
Default value: `30`

`--num_candidates`
The number of candidate models with optimal initial parameters that should be run. 
Default value: `10`

`--val_factor`
Sets K for the K-fold cross-validation. 
Default value: `10`.

`--epochs`
Sets the number of epochs (1 epoch = 1 forward and backward cycle through the entire training set) the network goes through. 
Default value: `100`

`--batch_size`
The size of one training batch (amount of sequences and the corresponding outputs). After one batch, the network will update its parameters through back-propagation. Default value: `32`

`--motif_length`
The length of each individual regulatory motif. 
Default value: `12`

`--motif_amount`
The amount of regulatory motifs to look for. 
Default value: `300`

`--sigma_motifs_min`
sigma_motifs will be drawn from loguniform(sigma_motifs_min, sigma_motifs_max).
Subsequently, sigma_motifs acts as the standard deviation of the normal
distribution from which convolutional kernel coefficients are drawn.
Default value: `1e-7`

`--sigma_motifs_max`
sigma_motifs will be drawn from loguniform(sigma_motifs_min, sigma_motifs_max).
Subsequently, sigma_motifs acts as the standard deviation of the normal
distribution from which convolutional kernel coefficients are drawn.
Default value: `1e-3`

`--sigma_net_min`
sigma_net will be drawn from loguniform(sigma_net_min, sigma_net_max).
Subsequently, sigma_net acts as the standard deviation of the normal
distribution from which neural network layer coefficients are drawn.
Default value: `1e-5`

`--sigma_net_max`
sigma_net will be drawn from loguniform(sigma_net_min, sigma_net_max).
Subsequently, sigma_net acts as the standard deviation of the normal
distribution from which neural network layer coefficients are drawn.
Default value: `1e-2`

`--epsilon_min`
The learning rate will be drawn from loguniform(epsilon_min, epsilon_max).
Default value: `5e-4`

`--epsilon_max`
The learning rate will be drawn from loguniform(epsilon_min, epsilon_max).
Default value: `5e-2`

`--opt`
The pytorch optimizer to use. Options include `SGD` and `Adam`. 
Default value: `SGD`

---------------------------------------------------------------------------------------------------
## Example output
```
N E X T F L O W  ~  version 20.04.1
Launching `scanem.nf` [berserk_swirles] - revision: bb8afb69fd
=========================================================================
=========================================================================

  scanem  v0.1 

=========================================================================

  run name             : example_dataset
  data path            : data/example_dataset.tsv.gz
  cell label data path : data/example_dataset_celltypes.tsv
  motif length         : 12
  amount of motifs     : 300
  epochs               : 30
  batch size           : 128
  K in K-fold CV       : 10
  number of cal        : 30 
  number of candidates : 10
  tomtom db file       : resources/Mus_musculus.meme  
  random seed          : 42     

=========================================================================
=========================================================================

         
executor >  lsf (1)
[62/3eb656] process > scanem             [  0%] 0 of 1
[-        ] process > tomtom             -
[-        ] process > motif_analysis     -
```

---------------------------------------------------------------------------------------------------

## Questions and errors
If you have any questions, or want to report an error, please use our [github issues page](https://github.com/jacobhepkema/scanem/issues)

---------------------------------------------------------------------------------------------------

## References
[1] Gupta, S., Stamatoyannopoulos, J. A., Bailey, T. L., & Noble, W. S. (2007). Quantifying similarity between motifs. Genome biology, 8(2), R24.

[2] Weirauch, M. T., Yang, A., Albu, M., Cote, A. G., Montenegro-Montero, A., Drewe, P., ... & Zheng, H. (2014). Determination and inference of eukaryotic transcription factor sequence specificity. Cell, 158(6), 1431-1443.

---------------------------------------------------------------------------------------------------

## License
MIT
