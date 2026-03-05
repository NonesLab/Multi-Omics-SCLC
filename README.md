# Multi-Omics-SCLC
Code to re-create results and figures from "Multi-omics reveals key molecular and cellular features of advanced small cell lung cancers associated with distinct therapeutic opportunities"
The scripts presented here are related to the article “Multi-omics reveals key molecular and cellular features of advanced small cell lung cancers associated with distinct therapeutic opportunities” Nones et al. 2026 Genome Medicine, xxx.


 ![Figure 6 overview](/../../blob/main/Figure6.png)
 


This figure was created in BioRender (www.biorender.com).
Summary:
This repository contains the workflow used to generate the data presented in the manuscript. This workflow uses idat.files (methylation arrays raw data) and outputs from other Tools and analysis workflows referenced below.
Requirements and inputs:
1)	R and suit of publicly available packages.
2)	idat.files (raw data from Illumina EPIC arrays) - GSE307029
3)	Clinical Information – Supplementary Table 1
4)	RNASeq counts – TMM counts from EGAS00001007832
5)	Gene sets from GSEA-MSigDB database (https://www.gsea-msigdb.org/gsea/msigdb)
6)	WGS somatic vcf files (EGAS00001007832) and copy number (output) estimated by ascatNgs ( Raine et al. Curr Protoc Bioinformatics. 2016;56:15 9 1- 9 7)
7)	Telomere length (output) estimated by qmotif tool (Holmes et al Bioinform Adv. 2022;2(1):vbac005)
8)	Focal amplifications and ecDNA  identified by AmpliconSuite v1.3.4 pipeline (Deshpande et al, Nat Commun. 2019;10(1):392) https://github.com/AmpliconSuite/AmpliconSuite-pipeline/
9)	cfDNA somatic mutations as reported by DRAGEN secondary analysis software, v. 21.10.6, build 5661 (Illumina).
Workflow Overview:
Idat.files were imported, filtered and normalised using ChAMP package (v 2.36.0) (Morris  et al Bioinformatics. 2014;30(3):428-30).
Unsupervised hierarchical clustering using Euclidean distance and ward.D2 method
cell-type deconvolution was estimated using the MethylCIBERSORT package (Chakravarthy, et al Nat Commun. 2018;9(1):3220.)
GSEA was executed with the R package fgsea v 1.32.0.( https://github.com/alserglab/fgsea)
MSI was estimated with MSIsensor v0.2 (Niu et al. Bioinformatics. 2014;30(7):1015-6.)
Mutational Signatures - R package YAPSA v1.34.0 (Hübschmann et al Genes Chromosomes Cancer. 2021;60(5):314-31) (https://bioconductor.org/packages/YAPSA)  

Outputs:
The Workflow script produces tables and figures presented in the manuscript.
