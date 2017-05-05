Adapter Predict
===============

Predicts the barcoded adapter sequence range used in a given file.

Run Using
^^^^^^^^^

::

    produse adapter_predict

or

::

    python /path/to/ProDuSe/ProDuSe/adapter_predict.py


Parameters
^^^^^^^^^^

    :-i --input:
        Paired fastq files. Two files must be specified in this argument.
    :-m --max_adapter_length:
        Limit the adapter sequence prediction to this length

Additional Considerations
^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning:: Ensure your samples are de-multiplexed prior to running Adapter Predict. The adapter sequence is estimated from ALL reads in the fastq files.



    
