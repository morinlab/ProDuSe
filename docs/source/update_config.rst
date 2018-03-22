Update_Config
=============

Description
^^^^^^^^^^^

Coverts an older produse config file to a format compatible with this version of ProDuSe


Run Using
^^^^^^^^^

::

	produse update_config

or 
::

	python /path/to/ProDuSe/ProDuSe/UpdateConfig.py


Parameters
^^^^^^^^^^

	:-i/--input:
		The configuration file to be reformatted
	:-o/--output:
		Path to the new, reformatted configuration file

Additional Info
^^^^^^^^^^^^^^^

Duplicate parameters are not permitted in config files. If two parameters specify the same argument, one will be removed
automatically. If they specify different arguments, they will both be retained. One of the arguments must be removed manually.

Any parameters with no arguments are removed automatically

