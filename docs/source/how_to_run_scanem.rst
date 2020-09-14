How to run scanem
=================

Usage:

.. code-block:: bash

 nextflow run [options] scanem.nf [arg...]


Example command:

.. code-block:: bash
 
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


**scanem** requires a run name (:code:`--name`), a pooled dataset (:code:`--data`), and pool cell type annotations (:code:`--celldata`). 
For creating a pooled dataset, find the right guide in the `documentation <https://scanem.readthedocs.io/en/latest/index.html>`_.
Another important argument is :code:`--tomtom`, for specifying
the motif :code:`.meme` database file to align found motifs to. I have included :code:`Mus_musculus.meme` and :code:`Homo_sapiens.meme`
from CIS-BP[1] in the :code:`resources` directory. 

Before running **scanem** you will probably want to identify the best :code:`-profile` to run with. This will define the executor
used by **scanem**. `See this page <https://scanem.readthedocs.io/en/latest/profiles.html>`_ for options and customisation. 

Minimal command:

.. code-block:: bash
 
 nextflow run scanem.nf \
  --name test_run \
  --data data/mock_data.tsv \
  --celldata data/mock_data_colData.tsv \
  --tomtom resources/Mus_musculus.meme

**References**

[1] Weirauch, M. T., Yang, A., Albu, M., Cote, A. G., Montenegro-Montero, A., Drewe, P., ... & Zheng, H. (2014). Determination and inference of eukaryotic transcription factor sequence specificity. Cell, 158(6), 1431-1443.
