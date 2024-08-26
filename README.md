# Microbiome_Analysis_Tutorials
Introductory Tutorials For Analysis of Microbiome Data

This is a tutorial to analyze microbiome data with R. 

QIIME 2 data files are defined in qza format , which are actually just zipped files. Below are data from the Moving Pictures in QIIME 2.
1. rooted-tree.qza
2. sample-metadata.tsv
3. table.qza
4. taxonomy.qza

The demo data-set comes from the QIIME 2 tutorial - Moving Pictures (https://docs.qiime2.org/2024.5/tutorials/moving-pictures/).

The repository contains convert.sh (to transform these qza files to a R compatible format) as well as microbiome_analysis.R (analysis of microbiome data).

The R script starts from the processed output from metagenomic sequencing, i.e. feature matrix. 


The R script has all the processes and steps defined with a "#", for the user's convenience. 

All the image files (heatmaps, plots, and others) generated from the R script have been provided here.
