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
the motif :code:`.meme` database file to align found motifs to. I have included :code:`Mus_musculus.meme` and :code:`Homo_sapiens.meme`
from CIS-BP[1] in the :code:`resources` directory. 

Before running **scover** you will probably want to identify the best :code:`-profile` to run with. This will define the executor
used by **scover**. `See this page <https://scover.readthedocs.io/en/latest/profiles.html>`_ for options and customisation. 

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



**References**

[1] Weirauch, M. T., Yang, A., Albu, M., Cote, A. G., Montenegro-Montero, A., Drewe, P., ... & Zheng, H. (2014). Determination and inference of eukaryotic transcription factor sequence specificity. Cell, 158(6), 1431-1443.
