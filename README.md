#oppar
##Outlier Profile and Pathway Analysis in R

OPPAR is the R implementation of the modified Cancer Outlier Profile Analysis [(mCOPA)](http://jclinbioinformatics.biomedcentral.com/articles/10.1186/2043-9113-2-22) method, published by Wang et al. in 2012.


### Download Instructions

Enter the following lines in R/Rstudio to download and install the source code for the package from Github.
```{r}
# install.packages("devtools") # if not already installed
library(devtools)
install_github("DavisLaboratory/oppar", build_vignettes = TRUE)
library(oppar)
```


For a complete example workflow, and to learn about how `oppar` works, please refer to the [Wiki](https://github.com/DavisLaboratory/oppar/wiki).

To report bugs or for any other enquiries, please email hediyehzadeh.s@wehi.edu.au

#### Release notes:
The current version of the package does not include the `gsva` function for Gene Set Enrichment Analysis. However, examples are provided in the [Wiki](https://github.com/DavisLaboratory/oppar/wiki) which suggest how `oppar` can be used in conjuction with other packages such as `limma` for gene set testing.

