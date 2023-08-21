# ![example_plot](/images/dot_movestext.png) 

Single Cell Proportion Test in Python is designed to analyze the difference between the proprotion of cells in clusters between two scRNA-seq samples. A permutation test is used to calculate a p-value for each cluster, and a confidence interval for the magnitude difference is returned via bootstrapping.
Results can be visualized using a point range plot. View the [tutorial](https://github.com/SamuelAMiller1/py-scProportionTest/blob/main/tutorials/scPropTest_tutorial.ipynb) to get started. If you are working with an seurat object you may consider the [implementation in R](https://github.com/rpolicastro/scProportionTest/tree/master).

## Installation

Current release.
```
pip install py-scProportionTest
```
Development version.
```
pip install git+https://github.com/SamuelAMiller1/py-scProportionTest.git
```

## Citing

Miller SA, Policastro RA, Sriramkumar S, Lai T, Huntington TD, Ladaika CA, Kim D, Hao C, Zentner GE, O'Hagan HM. LSD1 and Aberrant DNA Methylation Mediate Persistence of Enteroendocrine Progenitors That Support BRAF-Mutant Colorectal Cancer. Cancer Res. 2021 Jul 15;81(14):3791-3805. doi: 10.1158/0008-5472.CAN-20-3562. Epub 2021 May 25. PMID: 34035083; PMCID: PMC8513805.
