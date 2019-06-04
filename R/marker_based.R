#' Estimate cell type proportions using first PC of expression matrix
#' 
#' @param x A sample by gene bulk expression matrix. Genes should be marker
#'   genes
#' @param weighted Boolean. If weighted=TRUE, multiply scaled gene expression by
#'   gene weights
#' @param w Numeric vector. Weights of genes 
#' @return ret List. Attribute \strong{pcs} contains matrix of PCs, where PC1
#'   should be used as estimates for cell type abundances
#'   Attribute \strong{sdev} contains eigenvalues of eigendecomposition of
#'   var-covar matrix. The 1st eigenvalue should explain most of the variance.
#'   Attribute \strong{genes} contains names of genes.
EstimatePCACellTypeProportions <- function(x, weighted=FALSE, w=NULL){
  x <- base::scale(x)
  if (weighted) {
    # Intersect gene names of weights and column names of x
    common.markers <- base::intersect( base::colnames(x), base::names(w) )
    # not sure if this will ever happen, since we check this at line 268
    #if ( length(common.markers) == 0 ) {
    #  base::stop(base::paste0("Genes from weights w do not match with column ",
    #                         "names of expression matrix x."))
    #}
    x <- x[,common.markers]
    w <- w[common.markers]
    wd <- base::diag(w)
    xw <- x %*% wd
    varcov <- t(xw) %*% xw
  }
  else {
    varcov <- t(x) %*% x
  }
  varcov.ed <- base::eigen(varcov)
  rot <- varcov.ed$vectors
  wpcs <- x %*% rot
  sds <- varcov.ed$values
  # x contains PCs, sdev contains eigenvalues of eigendecomposition
  return(list( pcs = wpcs, sdev = sds, genes=colnames(x)))
}

#' Get number of genes to use with weighted PCA
#'
#' @param x Numeric Matrix. A sample by gene expression matrix containing the
#'   marker genes.
#' @param w Numeric Vector. The weights of the genes that correspond to the
#'   columns of x.
#' @param min.gene Numeric. Minimum number of genes to consider as markers.
#' @param max.gene Numeric. Maximum number of genes to consider as markers.
#' @return best.n Numeric. Number of genes to use
GetNumGenesWeighted = function(x, w, min.gene = 25, max.gene = 200){
  max.gene = base::min(max.gene, base::ncol(x))
  ratios = base::lapply(min.gene:max.gene,
                        function(i) {
                          ret = EstimatePCACellTypeProportions(x[,1:i],
                                        weighted=TRUE,
                                        w=w[1:i])
                          vars = ret$sdev
                          vars[1] / vars[2]
                        })
  best.n = base::which.max(ratios) + min.gene - 1
  return(best.n)
}

#' Get number of genes to use with no weighted information
#'
#' @param x Numeric Matrix. A sample by gene expression matrix containing the
#'   marker genes.
#' @param min.gene Numeric. Minimum number of genes to consider as markers.
#' @param max.gene Numeric. Maximum number of genes to consider as markers.
#' @return best.n Numeric. Number of genes to use
GetNumGenes = function(x, min.gene = 25, max.gene = 200){
  max.gene = base::min(max.gene, base::ncol(x))
  ratios = base::lapply(min.gene:max.gene,
                        function(i) {
                          ret = EstimatePCACellTypeProportions(x[,1:i])
                          vars = ret$sdev
                          vars[1] / vars[2]
                        })
  best.n = base::which.max(ratios) + min.gene - 1
  return(best.n)
}

#' Get unique markers present in only 1 cell type
#'
#' Given a data frame of marker genes for cell types, 
#' returns a new data frame with non-unique markers removed.
#'
#' @param x Data frame. Contains column with marker gene names
#' @param gene_col Character string. Name of the column that contains
#'   the marker genes
#' @return x Data frame. Markers with non-unique markers removed
#'
GetUniqueMarkers <- function(x, gene_col="gene"){
  keep <- ! (duplicated(x[,gene_col], fromLast = FALSE) |
             duplicated(x[,gene_col], fromLast = TRUE) )
  return(x[keep,])
}

#' Correlate columns of data frame
#'
#' This function runs correlation between markers of a data frame or matrix, 
#' returning the values of the lower/upper triangular of the correlation matrix
#' in a vector.
#'
#' @param x Data frame or matrix. Column vectors are correlated
#' @param method Character string. Name of method passed to cor. 
#'   Pearson by default.
#' @return cors Numeric vector. Correlation coefficients of pairs
#'
CorTri <- function(x, method="pearson"){
  cors <- stats::cor(x, method=method)
  cors <- cors[base::lower.tri(cors)]
  return(cors)
}

#' Return cell type proportions from bulk
#'
#' Calculate cell type proportions from a data frame containing bulk expression
#' values. Uses PCA (weighted or regular) to estimate relative proportions
#' within each cell type.
#'
#' @param bulk Expression Set containing bulk data
#' @param cell_types Character vector. Names of cell types.
#' @param markers Data frame with columns specifying cluster and gene,
#'   and optionally a column for weights, typically the fold-change of the gene.
#'   Important that the genes for each cell type are row-sorted by signficance.
#' @param ct_col Character string. Column name specifying cluster/cell type
#'   corresponding to each marker gene in \strong{markers}. 
#' @param gene_col Character string. Column name specifying gene names in
#'   \strong{markers}.
#' @param min_gene Numeric. Min number of genes to use for each cell type.
#' @param max_gene Numeric. Max number of genes to use for each cell type.
#' @param weighted Boolean. Whether to use weights for gene prioritization
#' @param w_col Character string. Column name for weights, such as "avg_logFC",
#'  in \strong{markers}
#' @param verbose Boolean. Whether to print log info during decomposition.
#'   Errors will be printed regardless. 
#'
#' @return A List. Slot \strong{cors} contains list of vectors with correlation
#'   coefficients. Slot \strong{ctps} contains list of CTP objects returned by
#'   GetCTP
GetCTP <- function(bulk,
                   cell_types,
                   markers,
                   ct_col,
                   gene_col,
                   min_gene,
                   max_gene,
                   weighted,
                   w_col,
                   verbose){
  ctp = base::lapply(cell_types, function(ct){
                     # Get marker genes
                     markers_ct = markers[ markers[,ct_col] == ct , , drop=FALSE]
                     ctm = base::make.unique(as.character(markers_ct[, gene_col]))
                     # Get markers in common between bulk and markers data frame
                     common_markers = base::intersect(ctm, Biobase::featureNames(bulk))
                     if ( base::length(common_markers) == 0 ){
                       base::stop("No marker genes found in bulk expression data")
                     }
                     expr = base::t(Biobase::exprs(bulk)[common_markers,])
                     if ( base::ncol(expr) < min_gene ){
                       base::stop(base::paste0(base::sprintf("For cell type %s, There are less marker genes in ", ct),
                                               base::sprintf("the bulk expression set (%i) than the ", base::ncol(expr)),
                                               base::sprintf("minimum number of genes set (%i) ", min_gene),
                                               "for PCA-based decomposition\nSet the min_gene parameter to a lower integer."))
                     }
                     if (weighted){
                       # Get gene weights
                       ctw = markers_ct[, w_col]; names(ctw) = ctm; ctw = ctw[common_markers]
                       ng = GetNumGenesWeighted(expr, ctw, min_gene, max_gene) # Number of markers for PCA
                       expr = expr[,1:ng,drop=FALSE]
                       if (verbose){
                         base::cat(base::sprintf("Using %i genes for cell type %s; ", ng, ct))
                       }
                       ret = EstimatePCACellTypeProportions(expr, weighted=TRUE, w=ctw[1:ng])
                     }
                     else{
                       ng = GetNumGenes(expr, min_gene, max_gene)
                       expr = expr[,1:ng,drop=FALSE]
                       if (verbose){
                         base::cat(base::sprintf("Using %i genes for cell type %s; ", ng, ct))
                       }
                       ret = EstimatePCACellTypeProportions(expr)
                     }
                     # Flip the sign of the first PC if negatively correlated with most genes
                     cors = stats::cor(expr, ret$pcs[,1])
                     n_pos = sum(cors[,1] > 0)
                     if (n_pos/base::length(cors[,1]) < (0.5)) {
                       ret$pcs[,1] = ret$pcs[,1] * -1
                     }
                     if (verbose){
                       cors = stats::cor(expr, ret$pcs[,1]); n_pos = sum(cors[,1] > 0)
                       pct <- as.character(as.integer(100 * round(n_pos/base::length(cors[,1]), 2)))
                       clen <- as.character(base::length(cors[,1]))
                       base::cat(base::paste0(pct, "% of ", clen, " marker genes correlate positively with PC1 for cell type ", ct, "\n"))
                     }
                     return(ret)
                   })
  return(ctp)
}

#' Performs marker-based decomposition of bulk expression using marker genes
#' 
#' Estimates relative abundances of cell types from PCA-based decomposition.
#' Uses a list of marker genes to subset the expression data, and returns the 
#' first PC of each sub-matrix as the cell type fraction estimates.
#' Optionally, weights for each marker gene can be used to prioritize genes
#' that are highly expressed in the given cell type.
#'
#' Note that this method expects the input bulk data to be normalized, unlike
#' the reference-based method.
#' 
#' @param bulk.eset Expression Set. Normalized bulk expression data.
#' @param markers Data frame with columns specifying cluster and gene,
#'   and optionally a column for weights, typically the fold-change of the gene.
#'   Important that the genes for each cell type are row-sorted by signficance.
#' @param ct_col Character string. Column name specifying cluster/cell type
#'   corresponding to each marker gene in \strong{markers}. 
#' @param gene_col Character string. Column name specifying gene names in
#'   \strong{markers}.
#' @param min_gene Numeric. Min number of genes to use for each cell type.
#' @param max_gene Numeric. Max number of genes to use for each cell type.
#' @param weighted Boolean. Whether to use weights for gene prioritization
#' @param w_col Character string. Column name for weights, such as "avg_logFC",
#'  in \strong{markers}
#' @param unique_markers Boolean. If TRUE, subset markers to include only genes 
#'   that are markers for only one cell type
#' @param verbose Boolean. Whether to print log info during decomposition.
#'   Errors will be printed regardless. 
#' @return A List. Slot \strong{bulk.props} contains estimated relative cell
#'   type abundances. Slot \strong{var.explained} contains variance explained by
#'   first 20 PCs for cell type marker genes. Slot \strong{genes.used} contains
#'   vector of genes used for decomposition.
#' @export
MarkerBasedDecomposition <- function(bulk.eset, 
                                       markers,
                                       ct_col="cluster",
                                       gene_col="gene",
                                       min_gene = 5,
                                       max_gene = 200,
                                       weighted=FALSE,
                                       w_col = "avg_logFC",
                                       unique_markers = TRUE,
                                       verbose=TRUE){
  # Check input
  if ( ! methods::is(bulk.eset, "ExpressionSet") ) {
    base::stop("Expression data should be in ExpressionSet")
  }
  if (min_gene > max_gene){
    base::stop(base::paste0(base::sprintf("min_gene (set at %i) ", min_gene),
                            "must be less than or equal to max_gene ",
                            base::sprintf("(set at %i)\n", max_gene)))
  }

  # Get unique markers if applicable
  if (unique_markers){
    if (verbose) cat("Getting unique markers\n")
    markers <- GetUniqueMarkers(markers, gene_col=gene_col)
    n_genes_clust <- table(markers[,ct_col])
    if (any(n_genes_clust < min_gene)) {
      base::stop(paste0("Clusters must have a minimum of ",
                        base::as.character(min_gene),
                        " unique marker genes"))
    }
  }
  cg <- intersect(unique(markers[,gene_col]), rownames(bulk.eset))
  if (length(cg) == 0) {
    base::stop(cat("No overlapping genes between markers and bulk.eset\n"))
  }
  markers <- markers[markers[,gene_col] %in% cg,]

  # Set variables
  markers[,ct_col] <- as.character(markers[,ct_col])
  cell_types <- sort(unique(markers[,ct_col]))
  n_ct <- length(cell_types)
  n_s <- base::ncol(bulk.eset)
  if (verbose){
    base::cat(base::paste0("Estimating proportions for ",
                           base::sprintf("%i cell types in %i samples\n",
                                         n_ct, n_s)))
  }

  # Remove zero-variance genes
  bulk.eset <- FilterZeroVarianceGenes(bulk.eset, verbose)

  # Get cell type proportions
  ctp <- GetCTP(bulk.eset, cell_types, markers, ct_col, gene_col,
                min_gene, max_gene, weighted, w_col, verbose)

  names(ctp) <- cell_types
  ctp_pc1 <- base::lapply(ctp, function(x) x$pcs[,1])
  ctp_pc1 <- base::do.call(cbind, ctp_pc1)
  ctp_cors <- stats::cor(ctp_pc1)
  ctp_cors <- CorTri(ctp_cors)
  if (verbose & all(ctp_cors > 0)){
    base::cat("***\n")
    base::cat(paste0("WARN: All cell type proportion estimates are correlated ",
                     "positively with each other. Check to make sure the ",
                     "expression data is properly normalized\n"))
    base::cat("***\n")
  }

  ctp_pc1 <- base::t(ctp_pc1)
  ctp_varexpl <- base::sapply(ctp, function(x) x$sdev[1:20])
  rownames(ctp_varexpl) <- base::paste0("PC", base::as.character(1:20))
  genes_used <- base::lapply(ctp, function(x) x$genes)
  if (verbose) cat("Finished estimating cell type proportions using PCA\n")
  return(list(bulk.props=ctp_pc1,
              var.explained=ctp_varexpl,
              genes.used=genes_used))
}

