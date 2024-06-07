#' @title Consensus OPLS for Multi-Block Data Fusion
#' @description
#' Merging data from multiple sources is a relevant approach for 
#' comprehensively evaluating complex systems. However, the inherent problems 
#' encountered when analyzing single tables are amplified with the generation 
#' of multi-block datasets, and finding the relationships between data layers 
#' of increasing complexity constitutes a challenging task. For that purpose, 
#' a generic methodology is proposed by combining the strengths of established 
#' data analysis strategies, i.e. multi-block approaches and the OPLS framework
#' to provide an efficient tool for the fusion of data obtained from multiple 
#' sources. The package enables quick and efficient implementation of the 
#' consensus OPLS model for any horizontal multi-block data structure 
#' (observation-based matching). Moreover, it offers an interesting range of 
#' metrics and graphics to help to determine the optimal number of components 
#' and check the validity of the model through permutation tests. 
#' Interpretation tools include scores and loadings plots, as well as 
#' Variable Importance in Projection (VIP), and performance coefficients such 
#' as R2, Q2 and DQ2 coefficients.
#' 
#' This package uses functions from the K-OPLS package, developed by Max
#' Bylesjo, University of Umea, Judy Fonville and Mattias Rantalainen, Imperial
#' College.
#' 
#' Copyright (c) 2007-2010 Max Bylesjo, Judy Fonville and Mattias Rantalainen
#' 
#' This code has been extended and adapted under the terms of the GNU General
#' Public License version 2 as published by the Free Software Foundation.
#'
#' @aliases ConsensusOPLS-package
#' @docType package
#' 
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' Three-block omics data
#'
#' @name demo_3_Omics
#' @description A demonstration case study available from a public repository
#' of the National Cancer Institute, namely the NCI-60 data set, was used to 
#' illustrate the method's potential for omics data fusion. A subset of NCI-60
#' data (transcriptomics, proteomics and metabolomics) involving experimental 
#' data from 14 cancer cell lines from two tissue origins, i.e. colon and 
#' ovary, was used. The object proposed in this package contains, in a list, 
#' all the information needed to make a model: the three data blocks, a list 
#' of observation names (samples) and the binary response matrix Y.
#' @docType data
#' @author Boccard & Rutledge
#' @references J. Boccard and D.N. Rutledge. A consensus OPLS-DA strategy for
#' multiblock Omics data fusion. Analytica Chimica Acta, 769, 30-39, 2013.
#' 
NULL
