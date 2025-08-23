*stctm* retrievals: basic orienteering
======================================


Two types of retrievals
-----------------------

With ``stctm``, you can either fit **transmission spectra** to obtain constraints on the TLSE (assuming it can explain all of the spectral variations), or fit **stellar spectra** to infer the mixtures of spectral components that can reproduce an observed flux-calibrated stellar spectrum (``exotune`` sub-module).

* The basic setup of ``stctm`` allows you to obtain posterior distributions on stellar surface parameters and ranges of transmission spectra that best match an observed transmission spectrum from the effect of unocculted stellar surface heterogeneities **alone** (TLS retrieval). To run such a retrieval, use the examples in ``stellar_contamination_analysis/`` (results to be populated in a ``stellar_contamination_results/`` folder).
* For stellar spectrum retrievals with ``exotune``, use the examples in ``exotune_analysis/`` (results to be populated in a ``exotune_results/`` folder).

I provide below the basics of the folder structure that should be adopted for all retrieval types, and custom instructions are provided in other derivative pages.

Detailed instructions on how to run TLS retrievals are available :ref:`here <runningTLS>`, while the information to run stellar spectrum retrievals with *exotune* can be consulted :ref:`here <running_exotune>`.


Input structure
---------------

* The inputs required to run a retrieval are structured into a file with a ``.toml`` extension. Examples are provided in each of ``stellar_contamination_analysis/`` and ``exotune_analysis/`` for TLS and stellar spectrum retrievals respectively. You can refer to the relevant page for details on how to structure inputs for each type of retrievals to create your own TOML file.

* The relative paths are set up such that your run files and toml input files are stored in a folder of your choice (say, ``your_analysis_name/``), which will be a subfolder of either of the analysis directories depending on the type of retrieval: either ``stellar_contamination_results/your_analysis_name/`` or ``exotune_results/your_analysis_name/``. Before starting a retrieval of your own, I recommend creating a copy of the folders as they are set up in the ``example/`` directory.

* You will likely create one "master" TOML file per retrieval, which can be saved under the ``stellar_contamination_analysis/your_analysis_name/`` or the ``exotune_analysis/your_analysis_name/`` subfolder, depending on the retrieval type, with the possibility to make minor variations without creating a brand new file each time for different retrieval flavors! Again, this is documented in the section on TOML files specific to each type of retrievals.

Running a retrieval
-------------------

Typically, running a retrieval should only entail the following steps:

1. Navigate via the command line to the analysis directory where you have your TOML file.
2. Run using the neat ``stctm_TLS``/``stctm_exotune`` CLI utility and with the TOML file as an argument,

so in the command line, this looks like e.g.::

   stctm_TLS your_toml_file.toml

or for the stellar spectrum retrievals::

   stctm_exotune your_toml_file.toml

Typically, you should not need to edit these run scripts to run standard retrievals as intended by the current code structure.

3. Inspect the output files generated in the corresponding results folder (under ``stellar_contamination_results/`` or ``exotune_results/``).

Retrieval outputs
-----------------

The outputs will be saved in a custom folder created by the run file in ``stellar_contamination_results/`` or ``exotune_results/``. You will find more details in the following pages on the outputs produced for each retrieval type.
