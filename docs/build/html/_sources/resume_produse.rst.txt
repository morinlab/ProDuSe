Resume ProDuSe
==============

Description
^^^^^^^^^^^

Restarts a previously terminated instance of `run_produse`_

.. _run_produse: run_produse.html

Run Using
^^^^^^^^^

::

	produse resume_produse

or

::

	/path/to/ProDuSe/Clone/ProDuSe/ResumePipeline.py

Parameters
^^^^^^^^^^

	:-d --produse_dir:
		Path to the base ProDuSe analysis directory (usually named produse_analysis_directory)
	:-j --jobs:
		Number of samples to process in parallel

Additional Information
^^^^^^^^^^^^^^^^^^^^^^

`run_produse`_ automatically generates <task>_Complete file when each pipeline
component is completed for each sample. These files are placed in the "config"
directory of the coresponding sample. This script identifies samples in which not all <task>_Complete files have been generated, and resumes the analysis from there.

.. note:: If you wish to re-run a stage of the pipeline, simply remove the coresponding <task>_Complete file

.. _run_produse: run_produse.html