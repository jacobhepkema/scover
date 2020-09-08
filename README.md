# scanem 

This GitHub repository contains the files needed and instructions to run the Nextflow pipeline of **scanem**.

<img src="https://github.com/jacobhepkema/scanem/raw/master/scanem_logo.png" width=300 align=right>

What is **scanem**? **scanem** is a convolutional neural network (CNN) for *de novo* inference of *cis*-regulatory motifs by training on single-cell data. It finds weights for these motifs across pseudo-bulks (weighing their 'impact'). 

The network is written in pytorch, with the downstream analyses written in R (using ggplot for plotting). Running the network and running the downstream analysis is implemented in a Nextflow pipeline. Furthermore, motifs are aligned with Tomtom from the MEME suite[1].

**scanem** requires that the cells are annotated for cell type (or other category). Furthermore, TSS annotation is required, as the promoter sequences are obtained directly from genomic sequence relative to the TSS. 

---------------------------------------------------------------------------------------------------
## Workflow

There are three main steps in using **scanem**:
1. Preparing your dataset - see [this link for scRNA-seq](https://htmlpreview.github.io/?https://github.com/jacobhepkema/scanem/blob/master/guides/how_to_prepare_scanem_input_using_scRNA_seq.html) or [this link for scATAC-seq](https://htmlpreview.github.io/?https://github.com/jacobhepkema/scanem/blob/master/guides/how_to_prepare_scanem_input_using_scATAC_seq.html). For a very quick scRNA-seq data generation script, see [this link](https://htmlpreview.github.io/?https://github.com/jacobhepkema/scanem/blob/master/guides/how_to_prepare_scanem_dataset_using_create_dataset.html)
2. Training **scanem** - see the guide further on this page: [link](#training-scanem). 
3. Analysing the output - see [this link](guides/how_to_analyse_scanem_output.html)

## Install/Dependencies

Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) by following the instructions on [their website](https://www.nextflow.io/docs/latest/getstarted.html#installation). 
Move the nextflow file to a directory accessible by your `$PATH` variable so that nextflow can be run by typing `nextflow` rather than the full path to the executable. One way of doing this is with the following line:
```
export PATH="path/to/nextflow:$PATH"
```
Alternatively, manually add `export PATH="path/to/nextflow:$PATH"` to your `~/.bashrc` or `~/.bash_profile` file, depending on which file is used as a source for the command line environment. 

When running locally (in non-singularity mode), the following dependencies are required:
```
python 3.6 or above

UPDATE

R ...
```

To get **scanem**, clone this repository into your current working directory with the following command:
```
git clone https://github.com/jacobhepkema/scanem
```
Alternatively, download the repository directly and place in a folder of your choice. 
`cd` into the top directory (`/scanem`) to run **scanem**.


---------------------------------------------------------------------------------------------------
## Training scanem

Usage:
```
nextflow run [options] scanem.nf [arg...]
```

Example command:
```
nextflow run -profile singularity scanem.nf \
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

**scanem** requires a run name (`--name`), a pooled dataset (`--data`), and pool cell type annotations (`--celldata`). 
See the guides above on how to create the pooled dataset. Another important argument is `--tomtom`, for specifying
the motif `.meme` database file to align found motifs to. I've included `Mus_musculus.meme` and `Homo_sapiens.meme`
from CIS-BP[2] in the `resources` directory. 

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
keeps saying `[  0%]`. This is however because the training stage only contains one
neural network step, and it will switch to `[100%]` once the training is complete. 
There is a way to find out how far into training the network is, but it is not entirely straightforward:
Open up a second terminal window and `cd` into the `scanem/work` directory. From there, you 
can `cd` into the directory starting with the code shown in front of the specific 
task (in this case, `62/3eb656`). Then, by running `tail .command.log` you might get 
an idea of how far along the training is. 

---------------------------------------------------------------------------------------------------
**Options**:

`-profile`
Choose a specified profile configuration from the `/conf` directory. Options:
`-profile singularity` 

`-with-report name.html`
Generates a Nextflow report webpage with information on task run times, CPU usage, memory usage, and I/O. Note that this does _not_ include information on GPU usage.

<!-- Options include `local` (default option), `lsf`, and `docker`. The `lsf` profile  uses a [Singularity](https://sylabs.io/docs/) image (built from [this git repo](https://github.com/jacobhepkema/scanem-wip) using [Singularity Hub](https://singularity-hub.org/)). For more options, see the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html) or run `nextflow run -h`. The cache directory for the Singularity image can be set in `/conf/lsf.config` by setting `cacheDir = ...` -->

---------------------------------------------------------------------------------------------------
**Arguments**:

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
The amount of randomly intialized calibrations for hyperparameter optimization. 
Default value: `30`

`--num_candidates`
The amount of candidate models with optimal initial parameters that should be run. 
Default value: `10`

`--val_factor`
Sets the K in K-fold cross-validation. 
Default value: `10`.

`--epochs`
Sets the amount of epochs (1 epoch = 1 forward and backward cycle through the entire training set) the network goes through. 
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
[62/3eb656] process > train_scanem       [  0%] 0 of 1
[-        ] process > tomtom             -
[-        ] process > tomtom_allmotifs   -
[-        ] process > motif_analysis     -
```

---------------------------------------------------------------------------------------------------

## Questions and errors
If you have any questions, or want to report an error, please use our [github issues page](https://github.com/jacobhepkema/scanem_pytorch/issues)

---------------------------------------------------------------------------------------------------

## References
[1] Gupta, S., Stamatoyannopoulos, J. A., Bailey, T. L., & Noble, W. S. (2007). Quantifying similarity between motifs. Genome biology, 8(2), R24.

[2] Weirauch, M. T., Yang, A., Albu, M., Cote, A. G., Montenegro-Montero, A., Drewe, P., ... & Zheng, H. (2014). Determination and inference of eukaryotic transcription factor sequence specificity. Cell, 158(6), 1431-1443.

---------------------------------------------------------------------------------------------------

## License
MIT
