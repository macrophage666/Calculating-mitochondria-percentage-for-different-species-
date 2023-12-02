#Background: 
#Nowadays, single cell RNA seq (scRNA seq) has become a common tool to solve biological problems, enabling people to evaluate the genes' expression at single cell resolution.
#Given its powerful nature, tons of packages were developed for scRNA seq analysis.Among them, Seurat is one of the most widely used R package which cover almost all the necessary functions for analyzing.
#In Seurat pipeline, analyzing starts from quality control(QC), which will exclude the doublets and dead cells. 
#To filter the doublets, the count of unique molecular indentifier(UMI) will be used. Cause it is assumed that if there are two/multiple cells with single cell barcode, the UMI count for that cell barcode will be high.
#To filter the dead cells, the percentage of mitochondria(mt) will be used. Cause the dead cells or the cells with bad quality are assumed to have high mt percentage.

#Problem to solve:
#The UMI count will be originally stored in the metadata which can be directly used for QC filtering, however the mt percentage is not provided. 
#Therefore, there is a function in Seurat package, which can calculate the mt percentage and add it to the metadata of Seurat Object,but this function only work for mice and human samples.
#For the people who have the samples from other species (i.e.Monkey, zebrafish, guinea pig....), there will be a problem.

#Function Introduction:
#This function will calculate the raw count and percentage of mitochondria genes for the QC filtering of scRNA seq data.
#This function would be a comprehensive method for calculating the mitochondria percentage compared with the function implemented by Seurat pipeline,
#since it can be worked for any species provided in the argument.

#Function body with annotations
#Five arguments are used in these function, but it is not necessary to provide all of them.
#Object is mandatory, which should be the SeuratObject; 
#Species is optional, which is only used if the genes are annotated by using Ensemble reference gene
#Resources has to be specified. The two resources for reference gene are Refseq and Ensemble, which will have different gene id for mitochondria genes.
#gtf and seqname have to be specified, if the 'Refseq' is chosen for resources. The gtf file is downloadable in NCBI website and seqname for mt genes can also be found in NCBI.

library(AnnotationHub)#Annotationhub is used to query the database in ensemble
library(ensembldb)#To use the function genes(),which can retrieve the queried annotation data from ensemble  


CalcMTPercent <- function(object, species, resources, gtf, seqname){
  if (resources == 'Ensemble'){
    ah <- AnnotationHub()
    ahq <- query(ah,pattern = c(species, "EnsDb"), ignore.case = TRUE)#query by using provided info for species
    ahq %>% 
      mcols() #access the metadata columns. This will show all the databases that fit the provided patterns
    id <- ahq %>% #For each database, there will be an id assigned
      mcols() %>%
      rownames() %>% 
      tail(n = 1) #The latest updated database will at the end 
    edb <- ah[[id]] #Get the database by providing the index id
    annotations <- genes(edb, 
                         return.type = "data.frame")#This function will transform the ensemble object to a dataframe
    annotations <- annotations %>%
      dplyr::select(gene_id, gene_name, seq_name)#Pull out the columns that are needed
    mt <- annotations %>%
      dplyr::filter(seq_name == "MT") %>% #The mitochondria genes will have their seq_name as 'MT'
      dplyr::pull(gene_name) #Grab the gene_name for the gene with 'MT' as seq_name
    mt1 <- annotations %>%
      dplyr::filter(seq_name == "MT") %>%
      dplyr::pull(gene_id) #There are some mt genes without gene_name,instead, the gene_id will be used and presented in the metadata.
    mt <- c(mt, mt1) #Combine the gene_name and gene_id for all mt genes
  } else if (resources == 'Refseq') { #The case that the reference gene are from 'Refseq' instead of 'Ensemble'
      gtf <- as.data.frame(gtf) #Refseq does not develop package to process their genome database, so .gtf file has to be manually downloaded from NCBI
      mt <- gtf %>% #After converting the .gtf file to a data frame, it is able to pull out the mt gene by provided seqname
        dplyr::filter(seqnames == seqname) %>% #The seqnames can also be found in NCBI website, the Refseq will assign an unique seqname for mt genes
        dplyr::pull(gene_id)
      mt1 <- gtf %>%
        dplyr::filter(seqnames == seqname) %>%
        dplyr::pull(gene)
      mt <- c(mt,mt1)
      mt <- gsub('_','-',mt) #the underscore will be automatically transformed to "-" in Seurat, so this conversion is done here
  } 
  object@meta.data$mtUMI <- Matrix::colSums(object[which(rownames(object) %in% mt),], na.rm = T) #This is to calculate the UMI count for mt gene in each cell and add them to the metadata
  object@meta.data$mitoPercent <- object@meta.data$mtUMI*100/object@meta.data$nCount_RNA #By dividing the total UMI count by mt UMI count, we can get mt percentage, which also can be added to the metadata of the SeuratObject
  object <- object #update the changed SeuratObject
}

#Example of using this function
#The URL for dataset used in example: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195574

#Example1: The scRNA seq data from Rhesus Macaque, which used Ensemble to annotate genes
#The packages that are used for loading the data
library(Seurat)
library(hdf5r)
setwd('D:/data/Xiaoke Mid Term Project')
test_Ensmbl <- Read10X_h5('./test_Ensmbl.h5')##Load scRNA seq data
test_Ensmbl <- CreateSeuratObject(test_Ensmbl)#Create a SeuratObject
head(test_Ensmbl@meta.data)
test_Ensmbl <- CalcMTPercent(object = test_Ensmbl, species = 'Macaca mulatta', resources = 'Ensemble')
head(test_Ensmbl@meta.data)#Check if mt percentage is added to the metadata.See 'Result for test_Ensmbl.png'.

#Plot the QC
win.graph()
VlnPlot(test_Ensmbl, features = c('nCount_RNA','nFeature_RNA','mitoPercent'))

#Example2: The scRNA seq data from Rhesus Macaque, which used Refseq to annotate genes
library(Seurat)
library(hdf5r)
library(rtracklayer) #The package to read .gtf file
gtf <- rtracklayer::import("./genomic.gtf") #This .gtf file is downloaded from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003339765.1/
test_Refseq <- Read10X_h5('./test_Refseq.h5')##Load scRNA seq data
test_Refseq <- CreateSeuratObject(test_Refseq)#Create a SeuratObject
test_Refseq <- CalcMTPercent(object = test_Refseq, resources = 'Refseq', gtf = gtf, seqname = 'NC_005943.1')
head(test_Refseq)#Check if mt percentage is added to the metadata.See 'Result for test_Refseq.png'.

#Example3: The scRNA seq data from Human PBMC, which used Ensmble to annotate genes
test2_Ensmbl <- Read10X_h5('./test2_Ensmbl.h5')#This dataset was downloaded from: https://cf.10xgenomics.com/samples/cell-exp/7.0.1/SC3pv3_GEX_Human_PBMC/SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix.h5
test2_Ensmbl <- CreateSeuratObject(test2_Ensmbl)
#By querying Ensemble dataset
test2_Ensmbl <- CalcMTPercent(object = test2_Ensmbl, resources = 'Ensemble', species = 'Homo sapiens')
head(test2_Ensmbl)
