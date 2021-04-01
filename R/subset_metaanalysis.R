#' Generating the dataset for the meta-analysis
#'
#' Given a number of differential expression analyses it creates a dataset with the right format for the downstream analyses.
#' @param ... dataframes with the differential expression analyses to be included in the final dataset.
#' @param col_names character vector with all the possible column names that indicate the gene names in the input dataframes. If this column is absent the rownames will be used as gene names.
#' @param col_adjpvals character vector with all the possible column names that indicate the adjusted p-value in the input dataframes.
#' @param col_log2FC character vector with all the possible column names that indicate the log2 fold changes in the input dataframes.
#' @return  dataset list of lists, including a list called `names` containing sublists, each corresponding to a specific comparison from a study, of differentially expressed genes (each sublist is a character vector). The sublists are named with the names given to each of the objects given as differential expression analyses.
#' @examples
#' # import dataset with DESeq2 differential expression analysis
#' data("ipsc_deseq2")
#' # import dataset with limma differential expression analysis
#' data("thymus_limma")
#' list_array_2studies=create_dataset(thymus_limma, ipsc_deseq2,col_names="GeneName", col_adjpvals=c("adj.P.Val", "padj"), col_log2FC=c("logFC", "log2FoldChange"))
#' @export
create_dataset <- function(...,col_names="GeneName", col_adjpvals=c("adj.P.Val", "padj"), col_log2FC=c("logFC", "log2FoldChange")){
  
names_dataset <- as.character(substitute(list(...)))[-1L]
deas=list(...)
# ###function to get variable name as string
# get_var_name=function(object){
#   deparse(substitute(object))
# }

##initialize lists
names <- list()
adjpval <- list()
log2FC <- list()

for(i in seq_along(deas)){
  df <- deas[[i]]
  name_col <- intersect(col_names, colnames(df))
  log2FC_col <- intersect(col_log2FC, colnames(df))
  adjp_col <- intersect(col_adjpvals, colnames(df))
  
  if(length(name_col)==0){ #if there is no column for the gene name use rownames
    names <- append(names, list(rownames(df)))
  } else {names <- append(names, list(df[,name_col]))}
  
  log2FC <- append(log2FC, list(df[,log2FC_col]))
  adjpval <- append(adjpval, list(df[,adjp_col]))
  
}

names(names) <- names_dataset
names(adjpval) <- names_dataset
names(log2FC) <- names_dataset

list_array <- list(names, adjpval, log2FC)

names(list_array) <- c("names", "adjpval", "log2FC")
return(list_array)

}


#' Subsetting the dataset
#'
#' Subset DE analysis data for a bunch of comparisons based on fold change and p-values
#' @param dataset list of 3 lists, the first with gene IDs (called `names`), the second with the adjusted p-values for each genes (called `adjpval`) and the third with log2 fold changes (called `log2FC`). Each one of these three lists have the number of sublists corresponding to the results of the differential expression analyses for each contrasts corresponding respectively to the geneIDs, adjusted p-values, and log2FCs.
#' @param adjpval desired adjusted p-value threshold for subsetting the data
#' @param max_n_genes maximum number of genes to be kept. It is used only if for a given contrasts there are more than `max_n_genes` genes that pass the chosen filtered. Ranking is based on the adjusted p-values and ties are included so there might be more than `max_n_genes` included.
#' @param abslog2FC desired absolute log2 fold change threshold for subsetting the data
#' @return  A list of lists structured identically to `dataset` but where the genes that do not pass the filters are excluded.
#' @examples
#' # create a list of lists were only the first 500 most significant genes with adjusted p-value < 0.05 and fold change >1.5 or < -1.5 are included
#' data(list_array) #load data
#' list_array.05_fc1.5_max500 <- subset_metanalysis(dataset=list_array, adjpval = 0.05, abslog2FC = log2(1.5), max_n_genes = 500 )
#' @export
subset_metanalysis <- function(dataset, adjpval=0.1, max_n_genes=Inf, abslog2FC=0){
  if( !identical(names(dataset), c("names", "adjpval", "log2FC") )) stop("dataset  should have 'names', 'adjpval' and 'log2FC' as names")
  if( length(dataset$names) !=length(dataset$adjpval) | length(dataset$names) !=length(dataset$log2FC))  stop("dataset$names,   dataset$adjpval, and dataset$log2FC should have the same length")
  indexes1 <- lapply(dataset$adjpval, FUN=function(x) which(x < adjpval))
  nomi <- names(dataset$names)
  dataset$names <-  lapply(seq_along(dataset$names), FUN=function(x) (dataset$names[[x]][(indexes1[[x]])]))
  dataset$adjpval <- lapply(seq_along(dataset$adjpval), FUN=function(x) (dataset$adjpval[[x]][(indexes1[[x]])]))
  dataset$log2FC <-  lapply(seq_along(dataset$log2FC), FUN=function(x) (dataset$log2FC[[x]][(indexes1[[x]])]))
  names(dataset$names) <- nomi
  names(dataset$adjpval) <- nomi
  names(dataset$log2FC) <- nomi


  subset <- which(lapply(dataset$names, length)>0)
  dataset$names <- dataset$names[subset]
  dataset$adjpval <- dataset$adjpval[subset]
  dataset$log2FC <- dataset$log2FC[subset]
  nomi <- names(dataset$names)

  indexes1 <- lapply(dataset$log2FC, FUN=function(x) which((x > abslog2FC )| (x < -abslog2FC )))

  dataset$names <-  lapply(seq_along(dataset$names), FUN=function(x) (dataset$names[[x]][(indexes1[[x]])]))
  dataset$adjpval <-  lapply(seq_along(dataset$adjpval), FUN=function(x) (dataset$adjpval[[x]][(indexes1[[x]])]))
  dataset$log2FC <-  lapply(seq_along(dataset$log2FC), FUN=function(x) (dataset$log2FC[[x]][(indexes1[[x]])]))
  names(dataset$names) <- nomi
  names(dataset$adjpval) <- nomi
  names(dataset$log2FC) <- nomi

  subset <- which(lapply(dataset$names, length)>0)
  dataset$names <- dataset$names[subset]
  dataset$adjpval <- dataset$adjpval[subset]
  dataset$log2FC <- dataset$log2FC[subset]


  for(x in seq_along(dataset$adjpval)){
    if(!is.na(dataset$adjpval[[x]][max_n_genes])){
      max_genes <- max_n_genes
      a <- TRUE
      while(a){
        if(dataset$adjpval[[x]][max_genes] == dataset$adjpval[[x]][max_genes+1]){
          max_genes <- max_genes+1

        } else{a <- FALSE}}

      dataset$names[[x]] <-  dataset$names[[x]][1:max_genes][!is.na(dataset$names[[x]][1:max_genes])]
      dataset$adjpval[[x]] <- dataset$adjpval[[x]][1:max_genes][!is.na(dataset$adjpval[[x]][1:max_genes])]
      dataset$log2FC[[x]] <- dataset$log2FC[[x]][1:max_genes][!is.na(dataset$log2FC[[x]][1:max_genes])]

    }
  }
  return(dataset)
}



#' Find co-differentially expressed genes (co-DE genes)
#'
#' Find all genes that are found differentially expressed in the same comparison(s).
#' @param gene a gene (character) included in one of the sublists in dataset. Can also be a vector of more than one gene, in that case only co-DE events including the genes in the vector will be reported.
#' @param dataset list of lists, including a list called `names` containing sublists, each corresponding to a specific comparison from a study, of differentially expressed genes (each sublist is a character vector). The sublists should be named to identify the comparisons.
#' @return  a dataframe with 4 columns: `gene1` and `gene2` represent a co-DE gene pair; `pair_occurrency` indicates in how many comparisons a given co-DE gene pair is reported; `comparisons` lists these comparisons separated by ";". Results are sorted per `pair_occurrency`.
#' @examples
#' data(list_array) #load data
#' # create a list of lists were only the first 500 most significant genes with adjusted p-value < 0.05 and fold change >1.5 or < -1.5 are included
#' data(list_array) #load data
#' list_array.05_fc1.5_max500 <- subset_metanalysis(dataset=list_array, adjpval = 0.05, abslog2FC = log2(1.5), max_n_genes = 500 )
#' #find co-DE genes for the DYRK1A gene (ensemblID="ENSG00000157540")
#' head(find_coDE(gene="ENSG00000157540",dataset=list_array.05_fc1.5_max500))
#' @export
find_coDE <- function(gene, dataset){

  suppressWarnings(if(!(gene %in% (unlist(dataset)))) stop("the provided gene ID is not present in your datasets"))
  if(!("names" %in% names(dataset))) stop("there is no sublist called 'names'")
  if(is.null(names(dataset$names))) stop("your sublists should be named")
  
  
  ### subset only for dataset containing the gene queried
  indeces <- rep(NA, length(dataset$names))
  for(i in seq_along(dataset$names)){ 
    indeces[i] <- sum(gene %in% dataset$names[[i]])
    }
  dataset$names <- dataset$names[!is.na(indeces)]
  dataset$adjpval <- dataset$adjpval[!is.na(indeces)]
  dataset$log2FC <- dataset$log2FC[!is.na(indeces)]
  
  cat("calculating combinations...\n")
  combinazioni <- lapply(dataset$names, FUN = function(x)(if(length(x)>1)(t(utils::combn(x,2))))) ## create all pairwise combination
  combinazioni <- combinazioni[!(unlist(lapply(combinazioni, is.null)))]
  dimensioni <- unlist(lapply(combinazioni, FUN=function(x) dim(x)[1]))

  int_matrix <- as.data.frame(matrix(rep(NA, 3*sum(dimensioni)), ncol=3))
  cat("creating gene co-DE matrix...\n")
  contatore <- 1
  for(i in seq_along(dimensioni)){

    int_matrix[contatore:(contatore+dimensioni[i]-1),c(1,2)] <- combinazioni[[i]]
    int_matrix[contatore:(contatore+dimensioni[i]-1),3] <- names(dimensioni)[i]
    contatore <- contatore+dimensioni[i]
  }

  index1 <- which(int_matrix[,1] %in% gene) ## rows with the given gene as first iteractor
  index2 <- which(int_matrix[,2] %in% gene) ## rows with the given gene as second interactor

  ### only subset matrix with the rows containing the gene of interests
  int_matrix_subset <- int_matrix[sort(unique(c(index1, index2))),]

  int_matrix_nometadata <- int_matrix_subset[,c(1,2)]

  #order couple of interaction alphabetically
  int_matrix_nometadata <- as.list(as.data.frame(t(int_matrix_nometadata)))
  int_matrix_nometadata <- lapply(int_matrix_nometadata, sort)
  unlisted <- unlist(int_matrix_nometadata)
  unlisted <- as.character(unlisted)
  matrix_int_nomet <- matrix(unlisted, ncol=2, byrow=TRUE)
  matrix_int_nomet <- as.data.frame(matrix_int_nomet)
  colnames(matrix_int_nomet) <- c("gene1", "gene2")
  matrix_int_nomet$comparisons <- int_matrix_subset$V3
  matrix_int_nomet$pair <- paste(matrix_int_nomet$gene1, matrix_int_nomet$gene2, sep=".")

  ranking <- sort(table(matrix_int_nomet$pair), decreasing=TRUE)

  summary_df <- data.frame("pair"=names(ranking),"pair_occurrency"=as.character(ranking))
  summary_df <- merge(summary_df, matrix_int_nomet)
  
  if(length(gene)>1){
    indici1 <- which(summary_df$gene1 %in% gene)
    indici2 <- which(summary_df$gene2 %in% gene)
    if(length(unique(intersect(indici1, indici2)))==0) stop("No interaction among the selected genes. Try to query genes individually.")
    summary_df <- summary_df[unique(intersect(indici1, indici2)),]
  }
  
  df1=stats::aggregate(comparisons~pair, summary_df, FUN=function(x)(paste(x, collapse=";")))
  
  summary_df=summary_df[,1:4]
  summary_df=summary_df[which(!duplicated(summary_df$pair)),]
  
  summary_df=merge(summary_df, df1)
  
  summary_df <- summary_df[order(summary_df$pair_occurrency, decreasing=TRUE),]
  
  summary_df <- summary_df[, c("gene1", "gene2", "pair_occurrency", "comparisons")]

  return(summary_df)

}






#' Retrieve fold changes and adjusted-pvalues
#'
#' Retrieve for a given genes all the log2FCs and adjusted-p-values across the different comparison(s)
#' @param gene a gene (character) included in one of the sublists in `dataset`.
#' @param dataset list of 3 lists, the first with gene IDs (called `names`), the second with the adjusted p-values for each genes (called `adjpval`) and the third with log2 fold changes (called `log2FC`). Each one of these three lists have a number of sublists corresponding to the gene IDs, adjusted p-values, and log2FCs, respectively from the differential expression analyses of each comparison.
#' @return  a dataframe with 2 columns: `log2FC` and `adjpval`. Each row is a different comparison where the given gene is found DE.
#' @examples
#' # create a list of lists were only the first 500 most significant genes with adjusted p-value < 0.05 and fold change >1.5 or < -1.5 are included
#' data(list_array) #load data
#' list_array.05_fc1.5_max500 <- subset_metanalysis(dataset=list_array, adjpval = 0.05, abslog2FC = log2(1.5), max_n_genes = 500 )
#' #retrieve log2FCs and adjusted-p-values for the DYRK1A gene (ensemblID="ENSG00000157540")
#' find_fold_changes_and_pvalues(gene="ENSG00000157540",dataset=list_array.05_fc1.5_max500)
#' @export
find_fold_changes_and_pvalues <- function(gene, dataset){
  if(!(gene %in% (unlist(dataset)))) stop("the provided gene ID is not present in your datasets")
  if( !identical(names(dataset), c("names", "adjpval", "log2FC") )) stop("dataset  should have 'names', 'adjpval' and 'log2FC' as names")
  if( length(dataset$names) !=length(dataset$adjpval) | length(dataset$names) !=length(dataset$log2FC))  stop("dataset$names,   dataset$adjpval, and dataset$log2FC should have the same length")
  
  dataset_true <- dataset
  dataset_true$names <- dataset$names[unlist(lapply(dataset$names, FUN=function(x) (gene %in% x)))]
  dataset_true$adjpval <- dataset$adjpval[unlist(lapply(dataset$names, FUN=function(x) (gene %in% x)))]
  dataset_true$log2FC <- dataset$log2FC[unlist(lapply(dataset$names, FUN=function(x) (gene %in% x)))]


  output <- unlist(lapply(seq_along(dataset_true$log2FC), FUN=function(x)(dataset_true$log2FC[[x]][which(dataset_true$names[[x]]==gene)])))

  output2 <- unlist(lapply(seq_along(dataset_true$adjpval), FUN=function(x)(dataset_true$adjpval[[x]][which(dataset_true$names[[x]]==gene)])))


  results <- as.data.frame(cbind(output, output2))
  rownames(results) <- names(dataset_true$names)
  colnames(results) <- c("log2FC", "adjpval")

  results <- results[order(results$adjpval, decreasing = FALSE),]
  return(results)
}




#' Retrieve statistics for all or a subset of genes present in the meta-analysis
#'
#' For each gene retrieve mean, median, and standard deviation of the log2FC and the number of times it is ecountered differentially expressed in the various comparisons.
#' @param dataset list of 3 lists, the first with gene IDs (called `names`), the second with the adjusted p-values for each genes (called `adjpval`) and the third with log2 fold changes (called `log2FC`). Each one of these three lists have a number of sublists corresponding to the gene IDs, adjusted p-values, and log2FCs, respectively from the differential expression analyses of each comparison.
#' @param genes vector of gene names. Limit the results only to the elements of `genes`. Cannot be set together with `top`.
#' @param top integer. Limit the results to the top number of genes with the highest occurence in the dataset. Cannot be set together with `genes`.
#' @return  a dataframe with 6 columns and as many rows as genes. The column `names` reports the gene names used in the dataset. `occurences` indicates the number of times a gene is found differentially expressed. `median_log2FC`, `mean_log2FC` and `sd_log2FC`: median, mean, and standard deviations of all the log2 fold changes in which the gene is encountered across the comparisons. `pseudotscore`: this statistics is calculated as `mean_log2FC/(sd_log2FC/sqrt(occurences))`. Genes with a high pseudo-t-score have high absolute fold change that does not vary across comparisons and are found differentially expressed in a relatively high number of comparisons.
#' @examples
#' # create a list of lists were only the first 500 most significant genes with adjusted p-value < 0.05 and fold change >1.5 or < -1.5 are included
#' data(list_array) #load data
#' list_array.05_fc1.5_max500 <- subset_metanalysis(dataset=list_array, adjpval = 0.05, abslog2FC = log2(1.5), max_n_genes = 500 )
#' #retrieve the statistics for the most frequent 500 DE genes
#' find_gene_statistics(dataset=list_array.05_fc1.5_max500, top=500)
#' @export
find_gene_statistics <- function(dataset, genes=NULL, top=NULL){
  if( !identical(names(dataset), c("names", "adjpval", "log2FC") )) stop("dataset  should have 'names', 'adjpval' and 'log2FC' as names")
  if( length(dataset$names) !=length(dataset$adjpval) | length(dataset$names) !=length(dataset$log2FC))  stop("dataset$names,   dataset$adjpval, and dataset$log2FC should have the same length")
  
  gene_frequency <- sort(table(unlist(dataset$names)), decreasing=TRUE)
  gene_matrix <- data.frame(names=names(gene_frequency), occurences=as.numeric(gene_frequency))

suppressWarnings(if(!(is.null(genes))){
  if(!(is.null(top))) stop("You can set whether 'top' or 'genes', not both")
  if(sum(genes %in% gene_matrix$names)==0) stop("None of your genes is present in the dataset")
  genes <- intersect(genes, gene_matrix$names)
  gene_matrix <- gene_matrix[which(gene_matrix$names %in% genes),]
})

suppressWarnings(if(!(is.null(top))){
  if(!(is.null(genes))) stop("You can set whether 'top' or 'genes', not both")
  if(top>dim(gene_matrix)[1]) stop("top cannot be higher than the total number of genes in the dataset")
  gene_matrix <- gene_matrix[1:top,]
})

gene_matrix$median_log2FC <- NA
gene_matrix$sd_log2FC <- NA
gene_matrix$mean_log2FC <- NA

cat("calculating statistics...\n")
notify <- seq(0, dim(gene_matrix)[1], 100)
for(i in seq_len(dim(gene_matrix)[1])){
  gene <- gene_matrix$names[i]
  fcs <- find_fold_changes_and_pvalues(gene, dataset)
  gene_matrix$median_log2FC[i] <- stats::median(fcs$log2FC, na.rm=TRUE)
  gene_matrix$sd_log2FC[i] <- stats::sd(fcs$log2FC, na.rm=TRUE)
  gene_matrix$mean_log2FC[i] <- mean(fcs$log2FC, na.rm=TRUE)
  if(i %in% notify){cat(i, 'of', dim(gene_matrix)[1],'completed...', '\n')}
}

gene_matrix$pseudotscore <- gene_matrix$mean_log2FC/(gene_matrix$sd_log2FC/sqrt(gene_matrix$occurences))
cat("Done!\n")

return(gene_matrix)
}




#' Plot a network of co-DE genes
#'
#' Plot a network where each node is a gene and each edge represent an event of co-differential expression.
#' @param dataset list of 3 lists, the first with gene IDs (called `names`), the second with the adjusted p-values for each genes (called `adjpval`) and the third with log2 fold changes (called `log2FC`). Each one of these three lists have a number of sublists corresponding to the gene IDs, adjusted p-values, and log2FCs, respectively from the differential expression analyses of each comparison.
#' @param genes it is taken into consideration only when `dataset` is not `NULL`. Limit the search only to the desired genes. Cannot be set together with `top`.
#' @param top it is taken into consideration only when `dataset` is not `NULL`. Integer that limits the results to the top number of genes with the highest frequence in the dataset. Cannot be set together with `genes`.
#' @param gene_stats_matrix dataframe of statistics only used when `dataset=NULL` and together with `code_matrix`. Can be produced using the function `find_gene_statistics`.
#' @param code_matrix dataframe of co-DE events only used when `dataset=NULL` and together with `gene_stats_matrix`. Can be produced using the function `find_coDE`. The gene should match the one present in `gene_stats_matrix`.
#' @param node_color color nodes based on the genes' mean (`"mean"`) or median (`"median"`) log2 fold changes or based on their pseudo-t-scores (`"pseudotscore"`).
#' @param node.names whether to print (`=TRUE`) or not (`=FALSE`) node names.
#' @return  plot a network where each node is a gene (node size corresponds to the number of comparisons in which a given gene is found co-DE) and each edge represent an event of co-differential expression (edge thickeness corresponds to the number of comparison where the genes are found co-DE). Nodes are colored based on log2FC (median or mean across comparisons) or pseudo-t-score (for more information check documentation for `find_gene_statistics`). Also returns the co-DE graph object produced by the igraph package.
#' @examples
#' # create a list of lists were only the first 500 most significant genes with adjusted p-value < 0.05 and fold change >1.5 or < -1.5 are included
#' data(list_array) #load data
#' list_array.05_fc1.5_max500 <- subset_metanalysis(dataset=list_array, adjpval = 0.05, abslog2FC = log2(1.5), max_n_genes = 500)
#' #plot the co-DE network for the most frequent 10 genes and color the nodes based on their mean log2FC
#' plot_gene_coDE_newtork(dataset=list_array.05_fc1.5_max500, top = 10, node_color="mean")
#' @export
plot_gene_coDE_newtork <- function(dataset=NULL, genes=NULL, top=NULL, code_matrix=NULL, gene_stats_matrix=NULL,node_color=c("mean", "median", "pseudotscore"), node.names=TRUE){
  if(!(node_color %in% c("mean", "median", "pseudotscore"))) stop("node_color should be one o 'mean', 'median' or 'pseudotscore'")
  if (is.null(code_matrix) & is.null(gene_stats_matrix)){
    if(is.null(dataset)) stop("You should define whether a code_matrix and a gene_stats_matrix or you should provide a dataset")
    gene_stats <- find_gene_statistics(dataset=dataset, genes=genes, top=top)
    codes <- find_coDE(gene = gene_stats$names, dataset=dataset)
  } else{
    gene_stats <- gene_stats_matrix
    codes <- code_matrix
  }

  rownames(gene_stats) <- gene_stats$names
  graph_code <- igraph::graph_from_data_frame(codes[,1:2], directed = FALSE, vertices = NULL)
  igraph::V(graph_code)$size <- gene_stats[igraph::V(graph_code)$name,2]
  igraph::V(graph_code)$color <- "white"
  
  igraph::E(graph_code)$weight <- as.numeric(codes$pair_occurrency)
  igraph::E(graph_code)$width <- log2(igraph::E(graph_code)$weight)
  
  
  whitered <- grDevices::colorRampPalette(c("white", "red"))
  bluewhite <- grDevices::colorRampPalette(c("blue", "white"))
  
  if(node_color=="pseudotscore"){
  cutcolors_pos <- whitered(100)[as.numeric(cut(gene_stats$pseudotscore[which(gene_stats$mean_log2FC>=0)],breaks = 100))]
  
  cutcolors_neg <- bluewhite(100)[as.numeric(cut(gene_stats$pseudotscore[which(gene_stats$mean_log2FC<0)],breaks = 100))]
  
  clim1 <- max(abs(min(gene_stats$pseudotscore)), abs(max(gene_stats$pseudotscore)))
  clim2 <- clim1*(-1)
  
  clabtag <- "pseudo t-score"
  }

  if(node_color=="mean"){
    cutcolors_pos <- whitered(100)[as.numeric(cut(gene_stats$mean_log2FC[which(gene_stats$mean_log2FC>=0)],breaks = 100))]
    
    cutcolors_neg <- bluewhite(100)[as.numeric(cut(gene_stats$mean_log2FC[which(gene_stats$mean_log2FC<0)],breaks = 100))]
    
    clim1 <- max(abs(min(gene_stats$mean_log2FC)), abs(max(gene_stats$mean_log2FC)))
    clim2 <- clim1*(-1)
    clabtag <- "mean log2FC"
    
  }
  
  if(node_color=="median"){
    cutcolors_pos <- whitered(100)[as.numeric(cut(gene_stats$median_log2FC[which(gene_stats$mean_log2FC>=0)],breaks = 100))]
    
    cutcolors_neg <- bluewhite(100)[as.numeric(cut(gene_stats$median_log2FC[which(gene_stats$mean_log2FC<0)],breaks = 100))]

    clim1 <- max(abs(min(gene_stats$median_log2FC)), abs(max(gene_stats$median_log2FC)))
    clim2 <- clim1*(-1)    
    clabtag <- "median log2FC"
    
  }
    
  cutcolors_df_pos <- data.frame("name"=gene_stats$names[which(gene_stats$mean_log2FC>=0)],"colors"=cutcolors_pos)
  cutcolors_df_neg <- data.frame("name"=gene_stats$names[which(gene_stats$mean_log2FC<0)],"colors"=cutcolors_neg)
  
  
  rownames(cutcolors_df_neg) <- cutcolors_df_neg$name
  rownames(cutcolors_df_pos) <- cutcolors_df_pos$name
  
  cutcolors_df <- rbind(cutcolors_df_neg,cutcolors_df_pos)
  cutcolors_df <- cutcolors_df[igraph::V(graph_code)$name,]
  cat("plotting the network...")
  
  layout <- igraph::layout_with_fr(graph_code, weights=igraph::E(graph_code)$weight)
  if(node.names==FALSE){
    graphics::par(mar=c(0,0,0,4), xpd=TRUE)
    plot(graph_code,layout=layout, vertex.color=as.character(cutcolors_df$colors), vertex.label=NA)
    plot3D::colkey(c(bluewhite(100), whitered(100)), add=TRUE, clab=clabtag, width=.7, length=.3, cex.clab=.75, line.clab=1
                   , clim=c(clim1,clim2)
                   # , shift=(-1)*0.2, dist=0.05 
    )
    
    graphics::legend("bottomright", 
           c("low", "high"), bty="n",pt.cex=c(1,3),
           pch=1, title="occurrences:",cex=1)
    graphics::legend("topright", 
           c("low", "high"), bty="n",lwd=c(1,3),
            cex=1, title="number of co-DE events:",pch="_")
  } else {
  graphics::par(mar=c(0,0,0,4),xpd=TRUE)
  plot(graph_code,layout=layout, vertex.color=as.character(cutcolors_df$colors))
  plot3D::colkey(c(bluewhite(100), whitered(100)), add=TRUE, clab=clabtag, width=.7, length=.3, cex.clab=.75, line.clab=1
                 , clim=c(clim1,clim2)
                 # , shift=(-1)*0.2, dist=0.05 
                 )
  
  graphics::legend("bottomright", 
         c("low", "high"), bty="n",pt.cex=c(1,3),
         pch=1, title="occurrences:",cex=1)
  graphics::legend("topright", 
         c("low", "high"), bty="n",lwd=c(1,3)
         , cex=1, title="number of co-DE events:",pch="_")
  }
  
  return(graph_code)
}




#library(devtools)
#devtools::install_github(build_vignettes = TRUE)
# load_all(".")

#library(roxygen2)
#roxygenise()
#devtools::build_manual()
#usethis::use_vignette("my-vignette")
#devtools::install_github(build_vignettes = TRUE)
# R CMD build --resave-data
# R CMD check /Users/ilario/OneDrive\ -\ CRG\ -\ Centre\ de\ Regulacio\ Genomica/Dropbox\ \(CRG\)/metaDEA
# usethis::use_testthat()
#use_git()
# devtools::build() to create a package bundle with the vignettes included.
#devtools::install_github(build_vignettes = TRUE)
# $ git config --global user.name "Ilarius"
# $ git config --global user.email "ilario.detoma@gmail.com"
# $ cd /Users/ilario/OneDrive\ -\ CRG\ -\ Centre\ de\ Regulacio\ Genomica/Dropbox\ \(CRG\)/metaDEA 
# $ git remote add origin https://github.com/Ilarius/metaDEA.git
# $ git pull origin master
# fatal: Couldn't find remote ref master
# $ git push -u origin master