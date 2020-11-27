How to run scover
=================

Usage:

.. code-block:: bash

 nextflow run [options] scover.nf [arg...]


Example command:

.. code-block:: bash
 
 nextflow run -profile lsf_gpu scover.nf \
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


**scover** requires a run name (:code:`--name`), a pooled dataset (:code:`--data`), and pool cell type annotations (:code:`--celldata`). 
For creating a pooled dataset, find the right guide in the `documentation <https://scover.readthedocs.io/en/latest/index.html>`_.
Another important argument is :code:`--tomtom`, for specifying
the motif :code:`.meme` database file to align the found motifs to. I have included :code:`Mus_musculus.meme` and :code:`Homo_sapiens.meme`
from CIS-BP[1] in the :code:`resources` directory. 

Before running **scover** you will probably want to identify the best :code:`-profile` to run with. This will define the executor
used by **scover**. `This page <profiles.html>`_ has more information on the profiles specified for **scover**. `See this page <https://scover.readthedocs.io/en/latest/profiles.html>`_ for options and customisation. 

Minimal command:

.. code-block:: bash
 
 nextflow run scover.nf \
  --name test_run \
  --data data/mock_data.tsv \
  --celldata data/mock_data_colData.tsv \
  --tomtom resources/Mus_musculus.meme


Important:
It might seem like it is stuck during training as the line with 

.. code-block:: bash

   [62/3eb656] process > scover             [  0%] 0 of 1

keeps saying :code:`[  0%]`. This happens because the training stage only contains one
neural network step, and it will switch to :code:`[100%]` once the training is complete. 
There is a way to find out how far into training the network is, but it is not entirely straightforward:
Open a second terminal window and :code:`cd` into the :code:`scover/work` directory. From there, you 
can :code:`cd` into the directory starting with the code shown in front of the specific 
task (in this case, :code:`62/3eb656`, so the directory will be something like :code:`scover/work/62/3eb656aae5e84c420b7aa267dfeb57`). 
Then, by running :code:`tail .command.log` you can get an idea of how far along the training is. 


Run an example dataset
######################

Try to run **scover** with an example dataset by `choosing your profile to run with <profiles.html>`_, making sure
you are in an environment that can run both Nextflow and Singularity, and then run the following from the :code:`scover` directory:

.. code-block:: bash

   nextflow run -profile local_gpu scover.nf \
    --data data/small_dataset.tsv \
    --celldata data/small_dataset_colData.tsv \
    --name example_run \
    --epochs 8 \
    --num_candidates 8 \
    --num_calibrations 4 \
    --tomtom resources/Mus_musculus.meme

Note that you should change :code:`local_gpu` to whatever profile you will use. This should create a folder in the :code:`scanem/output` directory named :code:`example_run` with all the output files. Note that most of the results will not make sense here, this is a dataset that contains the right data structure but not enough data for the network to converge properly. 



Arguments
#########

:code:`--name`
The name of the current experiment. 
Default value: :code:`experiment`

:code:`--data`
The path to the dataset (relative to :code:`scover.nf`) to train the network. Add your data in the example tab-separated format into the `/data` folder. 

:code:`--celldata`
The path to the dataset annotation (relative to :code:`scover.nf`) to annotate the cells. 
Add your data in the example format (see default file) into the :code:`/data` folder. 
For more information, see the data generation guides in the documentation. 

:code:`--tomtom`
The path to the :code:`.meme` format motif database to align found motifs to. 
Default value: :code:`resources/Mus_musculus.meme`.

:code:`--num_calibrations`
The number of randomly intialized calibrations for hyperparameter optimization. 
Default value: :code:`30`

:code:`--num_candidates`
The number of candidate models with optimal initial parameters that should be run. 
Default value: :code:`10`

:code:`--val_factor`
Sets K for the K-fold cross-validation. 
Default value: :code:`10`.

:code:`--epochs`
Sets the number of epochs (1 epoch = 1 forward and backward cycle through the entire training set) the network goes through. 
Default value: :code:`100`

:code:`--batch_size`
The size of one training batch (amount of sequences and the corresponding outputs). After one batch, the network will update its parameters through back-propagation. 
Default value: :code:`32`

:code:`--motif_length`
The length of each individual regulatory motif. 
Default value: :code:`12`

:code:`--motif_amount`
The amount of regulatory motifs to look for. 
Default value: :code:`300`

:code:`--sigma_motifs_min`
sigma_motifs will be drawn from loguniform(sigma_motifs_min, sigma_motifs_max).
Subsequently, sigma_motifs acts as the standard deviation of the normal
distribution from which convolutional kernel coefficients are drawn.
Default value: :code:`1e-7`

:code:`--sigma_motifs_max`
sigma_motifs will be drawn from loguniform(sigma_motifs_min, sigma_motifs_max).
Subsequently, sigma_motifs acts as the standard deviation of the normal
distribution from which convolutional kernel coefficients are drawn.
Default value: :code:`1e-3`

:code:`--sigma_net_min`
sigma_net will be drawn from loguniform(sigma_net_min, sigma_net_max).
Subsequently, sigma_net acts as the standard deviation of the normal
distribution from which neural network layer coefficients are drawn.
Default value: :code:`1e-5`

:code:`--sigma_net_max`
sigma_net will be drawn from loguniform(sigma_net_min, sigma_net_max).
Subsequently, sigma_net acts as the standard deviation of the normal
distribution from which neural network layer coefficients are drawn.
Default value: :code:`1e-2`

:code:`--epsilon_min`
The learning rate will be drawn from loguniform(epsilon_min, epsilon_max).
Default value: :code:`5e-4`

:code:`--epsilon_max`
The learning rate will be drawn from loguniform(epsilon_min, epsilon_max).
Default value: :code:`5e-2`

:code:`--opt`
The pytorch optimizer to use. Options include :code:`SGD` and :code:`Adam`. 
Default value: :code:`SGD`



**References**

[1] Weirauch, M. T., Yang, A., Albu, M., Cote, A. G., Montenegro-Montero, A., Drewe, P., ... & Zheng, H. (2014). Determination and inference of eukaryotic transcription factor sequence specificity. Cell, 158(6), 1431-1443.
