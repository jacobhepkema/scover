How to install scover prerequisites
===================================

**Scover** requires a couple of prerequisite programs. 

Nextflow
########

Install `Nextflow <https://www.nextflow.io/docs/latest/getstarted.html#installation>`_ by following the instructions on `their website <https://www.nextflow.io/docs/latest/getstarted.html#installation>`_.

Move the nextflow file to a directory accessible by your :code:`$PATH` variable so that nextflow can be run by typing :code:`nextflow` rather than the full path to the executable. One way of doing this is with the following line:

.. code-block:: bash
   
   export PATH="path/to/nextflow:$PATH"

Alternatively, manually add :code:`export PATH="path/to/nextflow:$PATH"` to your :code:`~/.bashrc` or :code:`~/.bash_profile` file, depending on which file is used as a source for the command line environment. 


Singularity
###########

Singularity is an alternative to Docker; Singularity images are environments that contain the right software packages/versions and I'm using them here to allow to run the different steps of :code:`scover` reproducibly.

If you don't have it yet, install `Singularity <https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps>`_ by following `the instructions on their website <https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps>`_. Try it out by typing:

.. code-block:: bash

   singularity --version


Scover 
######

To get :code:`scover`, clone the repository into your current working directory with the following command:

.. code-block:: bash

   git clone https://github.com/jacobhepkema/scover

Alternatively, download the repository directly and place in a folder of your choice. 
:code:`cd` into the top directory (:code:`/scover`) before you `run scover <how_to_run_scover.html>`_.
