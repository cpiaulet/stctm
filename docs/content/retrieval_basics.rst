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

* The inputs required to run a retrieval are structured into a file with a ``.ini`` extension. Examples are provided in each of ``stellar_contamination_analysis/`` and ``exotune_analysis/`` for TLS and stellar spectrum retrievals respectively. You can refer to the relevant page for details on how to structure inputs for each type of retrievals to create your own ini file.

* The relative paths are set up such that your run files and ini files are stored in a folder of your choice (say, ``your_analysis_name/``), which will be a subfolder of either of the analysis directories depending on the type of retrieval: either ``stellar_contamination_results/your_analysis_name/`` or ``exotune_results/your_analysis_name/``. Before starting a retrieval of your own, I recommend creating a copy of the folders as they are set up in the ``example/`` directory.

* You will likely create one "master" ini file per retrieval, which can be saved under the ``stellar_contamination_analysis/your_analysis_name/`` or the ``exotune_analysis/your_analysis_name/`` subfolder, depending on the retrieval type, with the possibility to make minor variations without creating a brand new file each time for different retrieval flavors! Again, this is documented in the section on ini files.

The ini files are provided as inputs to the template run files. I provide one template run file per retrieval type (``stellar_retrieval_v15_generic_runfile.py`` and ``exotune_runscript_v5_clean_20250422.py`` respectively. You should not need to edit these files as they serve three main purposes:
* Read the ini files and set up the retrievals.
* Set up the retrievals to run in serial or in parallel, if you have multiple cores available.
* Post-process the retrieval outputs to create results files you can use for post-processing, as well as diagnostics and paper-ready plots.

Running a retrieval
-------------------

Typically, running a retrieval should only entail the following steps:

1. Navigate via the command line to the analysis directory where you have your ini file and the template run file.
2. Run the template run file with the ini file as an argument,

e.g.::

   python stellar_retrieval_v15_generic_runfile.py your_ini_file.ini

or for the stellar spectrum retrievals::

   python exotune_runscript_v5_clean_20250422.py your_ini_file.ini

Typically, you should not need to edit these run scripts to run standard retrievals as intended by the current code structure.

3. Inspect the output files generated in the corresponding results folder (under ``stellar_contamination_results/`` or ``exotune_results/``).

Retrieval outputs
-----------------

The outputs will be saved in a custom folder created by the run file in ``stellar_contamination_results/`` or ``exotune_results/``. You will find more details in the following pages on the outputs produced for each retrieval type.
