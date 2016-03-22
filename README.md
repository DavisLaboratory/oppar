#oppar
**Outlier Profile and Pathway Analysis in R**

OPPAR is the R implementation of the modified Cancer Outlier Profile Analysis [(mCOPA)](http://jclinbioinformatics.biomedcentral.com/articles/10.1186/2043-9113-2-22) method, published by Wang et al. in 2012.


### Download Instructions

Enter the following lines in R/Rstudio to download and install the source code for the package from Github.
```{r}
# install.packages("devtools") # if not already installed
library(devtools)
install_github("DavisLaboratory/oppar", build_vignettes = TRUE)
library(oppar)
browseVignettes(package = "oppar") # opens up package vignette 
```


For a complete example workflow, and to learn about how `oppar` works, please refer to the package vignette or [Wiki](https://github.com/DavisLaboratory/oppar/wiki).

To report bugs or for any other enquiries, please email hediyehzadeh.s@wehi.edu.au

**Release notes:**
The modified `gsva` function for Gene Set Enrichment Analysis is now available. Several Examples are provided in the [Wiki](https://github.com/DavisLaboratory/oppar/wiki) for GSEA. It has also been demonstrated that how `oppar` can be used in conjuction with other packages such as `limma` for gene set testing. In addition, examples on how to apply `oppar` tools are provided in man/help pages. 

