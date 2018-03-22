ProDuSe: Process Duplex Sequence data
=====================================

    :Authors:
        | Christopher Rushton
        | Marco Albuquerque
    :Contact:
        ckrushto@sfu.ca.
    :Source Code:
        Github_
    :License:
        `GNU General Public License v3.0`_      

    .. _Github: https://github.com/morinlab/ProDuSe
    .. _GNU General Public License v3.0: License.html

ProDuSe is a variant caller designed for use with libraries prepared and sequenced using `semi-degenerate barcoded adapters`_. These
barcodes allow PCR and optical duplicates which were derived from the same starting molecule to be flagged and merged into a single
consensus sequence. The addition of a strand-specific tag also allows molecules which are derived from different parental strands
to be flagged, allowing variants at extremely low frequency (VAF~=0.1%) to be called confidently.

.. _semi-degenerate barcoded adapters: How_ProDuSe_Works.html

Quick Links
============
.. toctree::
   :maxdepth: 1

   Quick Start <Quick_Start>
   The ProDuSe Pipeline <The_ProDuSe_Pipeline>

Main Commands
=============
.. toctree::
   :maxdepth: 1

   run_produse
   adapter_predict
   Trim
   Collapse
   ClipOverlap
   Call

Additional Commands
===================
.. toctree::
  :maxdepth: 1

  update_config
  Train

Additional Links
================

.. toctree::
   :maxdepth: 1

   How ProDuSe Works <How_ProDuSe_Works>
   Config Files <Config_Files>
   License

