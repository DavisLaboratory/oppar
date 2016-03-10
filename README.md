# oppar
Outlier Profile and Pathway Analysis in R

The Cancer Outlier Profile Analysis (COPA) is a common analysis to identify genes that might be down-regulated or up-regulated only in a proportion of samples with the condition of interest. OPPAR is the R implementation of modified COPA (mCOPA) method, originally published by Chenwei Wang et la. in 2012. The aim is to identify genes that are outliers in samples with condition of interest, compared to normal samples. The methods implemented in OPPAR enable the users to perform the analysis in various ways, namely detecting outlier features in control versus condition samples (whether or not there is a information on subtypes), and detecting genes that are outlier in one subtype compared to the other samples, if the subtypes are known. 

OPPAR can also be used for Gene Set Enrichment Analysis (GSEA). Here, a modified version of [GSVA](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7) method is implemented. GSVA can be used to determine which samples in the study are enriched for gene expression signatures that are of interest. The `gsva` function in GSVA package returns an enrichment score for each sample, for the given signatures/gene sets. With the current implementation of the method, samples that strongly show enrichment for down(-regulated) gene expression signatures will receive negative scores.
However, Often it is in the interest of the biologists and researchers to get positive scores for samples that are enriched in both up and down signatures. Therefore, the `gsva` function has been modified to assign positive scores to samples that are enriched for the up-regulated and down-regulated gene expression signatures.

See Github Wiki for an example workflow.
