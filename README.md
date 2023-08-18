# Single Cell Proportion Test in Python

Analyze the difference between the proprotion of cells
in clusters between two scRNA-seq samples.
A permutation test is used to calculate a p-value for each cluster,
and a confidence interval for the magnitude difference is returned via bootstrapping.
There is also a function to generate a point range plot to display the results.

## Installation

Current release.
```
pip install py-scProportionTest
```
Development version.
```
pip install git+https://github.com/SamuelAMiller1/py-scProportionTest.git
```
