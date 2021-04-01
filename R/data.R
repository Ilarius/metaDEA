#' Example of a differential expression analysis generated with `DESeq2`.
#'
#' mRNA-seq profiling of iPS cells (4 euploid and 3 trisomy 21) derived from fibroblasts of monozygotic twins discordant for trisomy 21
#'
#' @format A data frame with 58336 rows and 6 variables:
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52249}
"ipsc_deseq2"


#' Example of a differential expression analysis generated with `limma`.
#'
#' Thymus transcriptomic profiles of patients with and without Down syndrome were compared in order to identify differentially expressed transcripts.
#'
#' @format A data frame with 20994 rows and 16 variables:
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69210}
"thymus_limma"


#' Example of a dataset generated from 124 different differential expression analyses
#' 
#' The dataset include microarray and RNA-seq experiments comparing trisomic samples versus euploid samples coming from differnt tissue and organs both in Homo sapines and Mus musculus, with the final goal to detect the genes most consistently differentially expressed in Down syndrome. Only genes with an adjusted p-value <0.1 are reported.
#'
#' @format list of 3 lists called `names`,  `adjpval`, and `log2FC` containing the same sublistis each with the same number of ordered elements:
#'   * _names_ list of 124 lists of gene names, each corresponding to a differnt comparison
#'   * _adjpval_ list of 124 lists of adjusted p-values, each corresponding to a differnt comparison
#'   * _log2FC_ list of 124 lists of log2 fold changes, each corresponding to a differnt comparison
"list_array"