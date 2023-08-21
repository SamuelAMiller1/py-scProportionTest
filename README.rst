.. raw:: html

   <picture>
       <source srcset="https://raw.githubusercontent.com/SamuelAMiller1/py-scProportionTest/main/images/proptest_webimage_dark.png" media="(prefers-color-scheme: dark)">
       <img src="https://raw.githubusercontent.com/SamuelAMiller1/py-scProportionTest/main/images/proptest_webimg_light.png" alt="Image">
   </picture>

**Single Cell Proportion Test in Python** is designed to analyze the difference between the proportion of cells in clusters between two scRNA-seq samples. A permutation test is used to calculate a p-value for each cluster, and a confidence interval for the magnitude difference is returned via bootstrapping. Results can be visualized using a point range plot. View the `tutorial <https://github.com/SamuelAMiller1/py-scProportionTest/blob/main/tutorials/scPropTest_tutorial.ipynb>`_ to get started. If you are working with an Seurat object, you may consider the `implementation in R <https://github.com/rpolicastro/scProportionTest/tree/master>`_.

Installation
============

Current release:

.. code-block:: bash

   pip install py-scProportionTest

Development version:

.. code-block:: bash

   pip install git+https://github.com/SamuelAMiller1/py-scProportionTest.git

Citing
======

Miller SA, Policastro RA, Sriramkumar S, Lai T, Huntington TD, Ladaika CA, Kim D, Hao C, Zentner GE, O'Hagan HM. LSD1 and Aberrant DNA Methylation Mediate Persistence of Enteroendocrine Progenitors That Support BRAF-Mutant Colorectal Cancer. Cancer Res. 2021 Jul 15;81(14):3791-3805. doi: 10.1158/0008-5472.CAN-20-3562. Epub 2021 May 25. PMID: 34035083; PMCID: PMC8513805.

