Profiles
========

.. code-block:: bash

 nextflow run -profile PROFILE_NAME scover.nf [args]

Using :code:`-profile`, you can choose a specified profile configuration from the :code:`/conf` directory. These profiles are included in the :code:`nextflow.config` file; if you are adding an option do not forget to add it to :code:`nextflow.config`. Options currently include:

* profile :code:`local_gpu` will run locally, and will use a GPU if available. Note that you need to have the right packages installed (see "How can I run scover using GPUs?" at the start of this README).
* profile :code:`local_cpu` will run locally, and it will also use Singularity for the first step in the workflow. This will only support CPUs.
* profile :code:`local_nosingularity` will run locally and will not use Singularity images.
* profile :code:`lsf_gpu` will run using the Platform LSF scheduler using GPUs. This might require some editing of the :code:`conf/lsf_gpu.conf` file to be compatible with GPU queues on your LSF setup. See `this page <https://www.nextflow.io/docs/latest/executor.html>`_ for more information on how to specify executors.
* profile :code:`lsf_cpu` will run using the Platform LSF scheduler using CPUs. This might require some editing of the :code:`conf/lsf_cpu.conf` file to be compatible with your LSF setup. See `this page <https://www.nextflow.io/docs/latest/executor.html>`_ for more information on how to specify executors.
