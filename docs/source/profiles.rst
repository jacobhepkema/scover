Choosing the right profile
==========================

.. code-block:: bash

 nextflow run -profile PROFILE_NAME scover.nf [args]

Nextflow configuration profiles are files that specify the executor (e.g. LSF/SLURM/local/AWS), requirements, and other things for the different parts of the Nextflow workflow. Currently, this repo includes configurations for running locally with GPU (:code:`local_gpu`), on LSF with GPU (:code:`lsf_gpu`), and some other options, including ones that run on CPU. However, given that scover already takes some time to run on GPU, I would advice against running it on CPU.  

Using :code:`-profile`, you can choose a specified profile configuration from the :code:`/conf` directory. These profiles are included in the :code:`nextflow.config` file; if you are adding an option do not forget to add it to :code:`nextflow.config`. Options currently include:

* profile :code:`local_gpu` will run locally, and will use a GPU if available. 
* profile :code:`local_cpu` will run locally, and it will only support CPUs.
* profile :code:`local_nosingularity` will run locally and will not use Singularity images.
* profile :code:`lsf_gpu` will run using the Platform LSF scheduler using GPUs. This might require some editing of the :code:`conf/lsf_gpu.conf` file to be compatible with GPU queues on your LSF setup. See `this page <https://www.nextflow.io/docs/latest/executor.html>`_ for more information on how to specify executors. The most important thing to change is the :code:`queue` option. 
* profile :code:`lsf_cpu` will run using the Platform LSF scheduler using CPUs. This might require some editing of the :code:`conf/lsf_cpu.conf` file to be compatible with your LSF setup. See `this page <https://www.nextflow.io/docs/latest/executor.html>`_ for more information on how to specify executors.


Custom profiles
###############

Depending on your setup, it might be necessary to specify a custom profile (for e.g. SLURM). 
For this, you can add a profile in the :code:`scanem/conf` directory (where the current :code:`.conf` files are located).
In this file, you can specify the executor and other requirements for the different parts of the workflow.
For information on how to specify these, see the `Nextflow documentation on configuration files <https://www.nextflow.io/docs/latest/config.html>`_ and the `Nextflow documentation on executors <https://www.nextflow.io/docs/latest/executor.html>`_.
Afterwards, you will also need to specify the profile in the :code:`scanem/nextflow.profile` file.

For example, if your profile is specified in the file :code:`scanem/conf/custom_profile.conf`, then you will need to add this to your :code:`scanem/nextflow.profile` file:

.. code-block:: bash

  custom_profile {
    includeConfig 'conf/custom_profile.conf'
  }

Then, you can `run scover <how_to_run_scover.html>`_ using :code:`-profile custom_profile`:

.. code-block:: bash

 nextflow run -profile custom_profile scover.nf [args]


