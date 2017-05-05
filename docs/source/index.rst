ProDuSe: Process Duplex Sequence data
=====================================

    :Authors:
        | Marco Albuquerque
        | Christopher Rushton
    :Contact:
        ckrushto@sfu.ca.
    :Source Code:
        Github_
    :License:
        `GNU General Public License v3.0`_      

    .. _Github: https://github.com/morinlab/ProDuSe
    .. _GNU General Public License v3.0: License.html

ProDuSe is a Python command-line variant caller designed for use on Illumina sequencing data obtained from libraries constructing using barcoded adapters. Using these adapters,
read duplicates are identified and collapsed into a single consensus sequence, correcting PCR and sequencing errors. In samples with a high duplication rate, this allows extremely rare
variants be be called confidently.

Useful Links
============
.. toctree::
   :maxdepth: 1

   Quick Start
   The ProDuSe Pipeline

ProDuSe Commands
================
.. toctree::
   :maxdepth: 1

   Config
   Trim
   bwa
   Collapse
   SplitMerge
   snv
   Filter

Other Commands
==============

.. toctree::
   :maxdepth: 1

   Adapter Predict

Additional Links
================

.. toctree::
   :maxdepth: 1

   License
