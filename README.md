# Calculating-mitochondria-percentage-for-different-species-
The function built up for calculating mitochondria percentage for any species used for scRNA seq and analyzed with Seurat packages. The R code not only includes the main R function but also includes examples explaining how to use this function. 
## Please make sure your Seurat version is v4 to use this function
## Background: 
Nowadays, single-cell RNA seq (scRNA seq) has become a common tool to solve biological problems, enabling people to evaluate the genes' expression at single-cell resolution.
Given its powerful nature, many packages were developed for scRNA seq analysis. Among them, Seurat is one of the most widely used R packages, which covers almost all the necessary functions for analyzing.
In the Seurat pipeline, analysis starts with quality control(QC), which will exclude doublets and dead cells. 
The count of unique molecular identifiers (UMI) will be used to filter the doublets. It is assumed that if there are two/multiple cells with single-cell barcodes, the UMI count for that cell barcode will be high.
The percentage of mitochondria(mt) will be used to filter the dead cells. Cause the dead cells or the cells with bad quality are assumed to have a high mt percentage.

## Problem to solve:
The UMI and gene count will be originally stored in the metadata, which can be directly used for QC filtering. However, the mt percentage is not provided. 
There is a function in the Seurat package that can calculate the mt percentage and add it to the metadata of the Seurat Object, but this function only works for mice and human samples.
For the people who have samples from other species (i.e., Monkeys, zebrafish, guinea pigs, etc.), there will be a problem.

# Function Introduction:
This function will calculate the raw count and percentage of mitochondria genes for the QC filtering of scRNA seq data.
This function would be a comprehensive method for calculating the mitochondria percentage compared with the function implemented by the Seurat pipeline.

## Function body with annotations
Five arguments are used in these functions, but providing all of them is unnecessary.
The object is mandatory, which should be the SeuratObject. 
Species are optional and only used if the genes are annotated using the Ensemble reference gene resources. The two resources for reference genes are Refseq and Ensemble, which will have different gene IDs for mitochondria genes.
If the 'RefSeq' is chosen for resources, gtf and seqname must be specified. The gif file is downloadable on the NCBI website, and seqname for mt genes can also be found in NCBI.
