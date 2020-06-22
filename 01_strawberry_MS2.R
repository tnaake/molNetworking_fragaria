## source functions for ms2 similarity
source("~/Projects/molNetworking_swM/01_ms2_similarity_functions.R")

## load libraries
library(MetCirc)
library(MSnbase)

## read files 
setwd("~/Projects/molNetworking_fragaria/")
aln_pos <- read.csv("MS2_data/pos/PeakID_0_20202201438.txt", skip = 4, 
    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
aln_neg <- read.csv("MS2_data/neg/PeakID_0_2020220167.txt", skip = 4,
    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## create reference spectra
aln_pos_ref <- createRefSpectra(aln_pos)
aln_neg_ref <- createRefSpectra(aln_neg)

## remove these MS2 that contain 5 or less mass peak
aln_pos_ref <- aln_pos_ref[unlist(lapply(aln_pos_ref, length)) > 5]
aln_neg_ref <- aln_neg_ref[unlist(lapply(aln_neg_ref, length)) > 5]

## construct Spectrum2
spl_pos_ref <- construct_Spectrum2(aln_pos_ref)
spl_neg_ref <- construct_Spectrum2(aln_neg_ref)

## create Spectra
spectra_pos_ref <- create_Spectra(spl_pos_ref)
spectra_neg_ref <- create_Spectra(spl_neg_ref)

## compare Spectra
simMat_pos_ref <- compare_Spectra(spectra_pos_ref, 
    fun = normalizeddotproduct, binSize = 0.01)
simMat_neg_ref <- compare_Spectra(spectra_neg_ref, 
    fun = normalizeddotproduct, binSize = 0.01)

## rename cols and rows of the similarity matrices
colnames(simMat_pos_ref) <- rownames(simMat_pos_ref) <- 
    paste(spectra_pos_ref@elementMetadata$precursorMz, 
          spectra_pos_ref@elementMetadata$rt, sep = "_")
colnames(simMat_neg_ref) <- rownames(simMat_neg_ref) <- 
    paste(spectra_neg_ref@elementMetadata$precursorMz, 
          spectra_neg_ref@elementMetadata$rt, sep = "_")

## save the objects
save(aln_pos_ref, aln_neg_ref, 
     file = "MS2_data/01_MS2_reference_spectrum_pos_neg.RData")
save(spl_pos_ref, spl_neg_ref, 
     file = "MS2_data/01_MS2_spectrum_list_reference_pos_neg.RData")
save(spectra_pos_ref, spectra_neg_ref,
     file = "MS2_data/01_MS2_spectra_reference_pos_neg.RData")
save(simMat_pos_ref, simMat_neg_ref, 
     file = "MS2_data/01_MS2_similarityMatrix_reference_pos_neg.RData")

## create the network

library(igraph)
g_pos <- graph_from_adjacency_matrix(simMat_pos_ref, weighted = TRUE)
g_neg <- graph_from_adjacency_matrix(simMat_neg_ref, weighted = TRUE)
