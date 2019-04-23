#' Find overlapping samples in single-cell and bulk data
#'
#' @param sc.eset Expression Set with single-cell data
#' @param bulk.eset Expression Set with bulk data
#' @param subject.names A character string. Name of phenoData attribute in
#'   sc.eset that indicates individual ID (that would be found in bulk.eset
#'   if overlapping)
#' @param verbose Boolean. Print logging info
#' @return samples A list with attributes \emph{overlapping} and
#'   \emph{remaining}. Each attribute refers to a character vector that lists
#'   the samples found in both datasets and samples found only in bulk,
#'   respectively
#' @examples
#' GetOverlappingSamples(sc.eset, bulk.eset, "SubjectName", TRUE)
GetOverlappingSamples <- function(sc.eset, bulk.eset, subject.names, verbose) {
  bulk.samples <- Biobase::sampleNames(bulk.eset)
  sc.samples <- base::levels(base::factor(sc.eset[[subject.names]]))
  overlapping.samples <- base::intersect(bulk.samples, sc.samples)
  if (base::length(overlapping.samples) == 0) {
    base::stop("No overlapping samples in bulk and single-cell expression.")
  }
  remaining.samples <- base::setdiff(Biobase::sampleNames(bulk.eset),
                                     overlapping.samples)
  if (base::length(remaining.samples) == 0) {
    base::stop("All samples have single-cell data, nothing to process.")
  }
  samples <- base::list("overlapping"=overlapping.samples,
                        "remaining"=remaining.samples)
  if (verbose) {
    n.overlapping <- base::length(samples$overlapping)
    n.remaining <- base::length(samples$remaining)
    base::cat(base::sprintf("Found %i samples", n.overlapping),
                            "with bulk and single-cell expression.\n",
                            sep=" ")
    base::cat(base::sprintf("Remaining %i", n.remaining),
                            "bulk expression samples will be deconvoluted.\n",
                            sep=" ")
  }
  return(samples)
}

#' Find overlapping genes in single-cell data, bulk data, and marker genes
#'
#' @param sc.eset Expression Set with single-cell data
#' @param bulk.eset Expression Set with bulk data
#' @param markers Character vector. List of relevant marker genes
#' @param verbose Boolean. Print logging info
#' @return overlapping.genes Character vector. List of genes found in markers
#'   and both datasets. 
#' @examples
#' GetOverlappingGenes(sc.eset,bulk.eset, c("NRGN", "GFAP", "SLC1A2"), TRUE)
GetOverlappingGenes <- function(sc.eset, bulk.eset, markers, verbose) {
  bulk.genes <- Biobase::featureNames(bulk.eset)
  sc.genes <- Biobase::featureNames(sc.eset)
  overlapping.genes <- base::intersect(bulk.genes, sc.genes)
  if (base::length(overlapping.genes) == 0) {
    base::stop(base::paste0("No overlapping genes found between bulk and ",
                           "single-cell expression."))
  }
  overlapping.genes <- base::intersect(overlapping.genes, markers)
  if (base::length(overlapping.genes) == 0) {
    base::stop(base::paste0("No marker genes found in both bulk and ",
                           "single-cell expression."))
  }
  if (verbose) {
    n.genes <- base::length(overlapping.genes)
    base::cat(base::sprintf("Using %i genes in both", n.genes),
                            "bulk and single-cell expression.\n",
                            sep=" ")
  }
  return(overlapping.genes)
}

#' Estimate cell type proportions using first PC of expression matrix
#' 
#' @param x A sample by gene bulk expression matrix. Genes should be marker genes
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
    if ( length(common.markers) == 0 ) {
      base::stop(base::paste0("Genes from weights w do not match with column ",
                             "names of expression matrix x."))
    }
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
#' @param gene_col Character string. Name of the column that contains the marker genes
#'
#' @return x Data frame. Markers with non-unique markers removed
#'
GetUniqueMarkers <- function(x, gene_col="gene"){
  keep <- ! ( duplicated(x[,gene_col], fromLast = FALSE) | duplicated(x[,gene_col], fromLast = TRUE) )
  return(x[keep,])
}

#' Correlate columns of data frame
#'
#' This function runs correlation between markers of a data frame or matrix, 
#' returning the values of the lower/upper triangular of the correlation matrix
#' in a vector.
#'
#' @param x Data frame or matrix. Column vectors are correlated
#' @param method Character string. Name of method passed to cor. Pearson by default
#'
#' @return cors Numeric vector. Correlation coefficients of pairs
#'
CorTri <- function(x, method="pearson"){
  cors <- cor(x, method=method)
  cors <- cors[lower.tri(cors)]
  return(cors)
}

#' Return cell type proportions from bulk
#'
#' Calculate cell type proportions from a data frame containing bulk expression values.
#' Uses PCA (weighted or regular) to estimate relative proportions within each cell type.
#'
#' @param bulk.eset Expression Set containing bulk data
#' @param cell_types Character vector. Names of cell types.
#' @param markers Data frame with columns specifying cluster and gene,
#'   and optionally a column for weights, typically the fold-change of the gene.
#'   Important that the genes for each cell type are row-sorted by signficance.
#' @param ct.col Character string. Column name specifying cluster/cell type
#'   corresponding to each marker gene in \strong{markers}. 
#' @param gene.col Character string. Column name specifying gene names in
#'   \strong{markers}.
#' @param min.gene Numeric. Min number of genes to use for each cell type.
#' @param max.gene Numeric. Max number of genes to use for each cell type.
#' @param weighted Boolean. Whether to use weights for gene prioritization
#' @param w.col Character string. Column name for weights, such as "avg_logFC",
#'  in \strong{markers}
#' @param verbose Boolean. Whether to print log info during deconvolution.
#'   Errors will be printed regardless. 
#'
#' @return A List. Slot \strong{cors} contains list of vectors with correlation
#'   coefficients. Slot \strong{ctps} contains list of CTP objects returned by GetCTP
#'
#'
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
                     ctm = base::make.unique(markers_ct[, gene_col])
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
                                               "for PCA-based deconvolution\nSet the min_gene parameter to a lower integer."))
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
                     cors = cor(expr, ret$pcs[,1])
                     n_pos = sum(cors[,1] > 0)
                     if (n_pos/base::length(cors[,1]) < (0.5)) ret$pcs[,1] = ret$pcs[,1] * -1
                     if (verbose){
                       cors = cor(expr, ret$pcs[,1]); n_pos = sum(cors[,1] > 0)
                       pct <- as.character(as.integer(100 * round(n_pos/base::length(cors[,1]), 2)))
                       clen <- as.character(base::length(cors[,1]))
                       base::cat(base::paste0(pct, "% of ", clen, " marker genes correlate positively with PC1 for cell type ", ct, "\n"))
                     }
                     return(ret)
                   })
  return(ctp)
}

#' Regress out global expression PCs from bulk.
#' 
#' This function normalizes bulk expression data by first calculating
#' global PCs and then regressing them out sequentially. If genes_subset 
#' is given, only normalizes and returns expression data for those listed 
#' genes, speeding up computation.
#' 
#' @param bulk.eset ExpressionSet. Contains bulk expression data to normalize
#' @param max_pc Numeric. Maximum number of PCs to regress out
#' @param genes_subset Character vector. Vector of genes to subset the ExpressionSet bulk object
#'
#' @return List with each element containing an ExpressionSet object of the normalized
#' bulk gene expression data, starting with 0 (non-normalized) up to and including 
#' max_pc
#'
NormalizeByPC <- function(bulk.eset, max_pc, genes_subset=NULL){
  bulk.eset <- FilterZeroVarianceGenes(bulk.eset, verbose=FALSE)
  bulk.eset <- FilterUnexpressedGenes(bulk.eset, verbose=FALSE)
  bulk <- Biobase::exprs(bulk.eset);
  global_pcs <- prcomp(t(bulk), retx=TRUE, scale.=TRUE)$x
  if (!is.null(genes_subset)){
    bulk <- bulk[genes_subset,]
  }
  bulks <- list()
  # names(bulks) <- paste0("PC", as.character(0:max_pc))
  bulks[["0"]] <- ExpressionSet(bulk)
  for (i in 1:max_pc){
    # pc_name <- paste0("PC", as.character(i))
    bulk <- t(apply(bulk, 1, function(x) resid(lm(x ~ global_pcs[,i,drop=FALSE]))))
    bulks[[as.character(i)]] <- ExpressionSet(bulk)
  }
  return(bulks)
}

#' Return number of PCs to use for PCA CTP estimation
#'
#' Return the minimum number of PCs to regress out based on 
#' CTP-CTP correlations. Given a list of vectors, return 
#' the number of PCs with less than pct_thresh percent of 
#' coefficients that are greater than 0.
#'
#' @param cors List. Contains correlation coefficients.
#' @param pct_thresh Numeric. Threshold of percent of coefficients that are greater than 0.
#' @param verbose Boolean. Print log info.
#'
#' @return best_pc_n Numeric. The number of PCs to regress out to remove global covariates.
GetNumPC <- function(cors, pct_thresh=0.5, verbose=TRUE){
  if ( pct_thresh < 0 | pct_thresh > 1) base::stop("pct_thresh should be between 0 and 1")
  n_gt_0 <- sapply(cors, function(x) table(x > 0))
  lengths <- sapply(cors, function(x) length(x))
  pct_gt_0 <- n_gt_0 / lengths
  good <- pct_gt_0 < pct_thresh
  best_pc_n <- which(good)[1] # Return first number of PCs with > 0 correlations < pct_thresh
  if (is.na(best_pc_n)) base::stop("No PC-normalized bulk expression data passes threshold criteria")
  if (verbose) cat(paste0("Regressed out ", as.character(best_pc_n), " global PCs for cell type proportion estimation\n"))
  return(best_pc_n)
}

#' Return cell type proportions from bulk by first normalizing for PCs
#'
#' Calculate cell type proportions from a data frame containing bulk expression values.
#' Uses PCA (weighted or regular) to estimate relative proportions within each cell type.
#' This function first calculates the number of global PCs to regress out from the data. 
#' It then uses the normalized expression data as input for PCA.
#'
#' @param bulk.eset Expression Set containing bulk data
#' @param max_pc Numeric. If normalized=FALSE, the maximum number of PCs to consider
#'   regressing out for normalization.
#' @param cell_types Character vector. Names of cell types.
#' @param markers Data frame with columns specifying cluster and gene,
#'   and optionally a column for weights, typically the fold-change of the gene.
#'   Important that the genes for each cell type are row-sorted by signficance.
#' @param ct.col Character string. Column name specifying cluster/cell type
#'   corresponding to each marker gene in \strong{markers}. 
#' @param gene.col Character string. Column name specifying gene names in
#'   \strong{markers}.
#' @param min.gene Numeric. Min number of genes to use for each cell type.
#' @param max.gene Numeric. Max number of genes to use for each cell type.
#' @param weighted Boolean. Whether to use weights for gene prioritization
#' @param w.col Character string. Column name for weights, such as "avg_logFC",
#'  in \strong{markers}
#' @param verbose Boolean. Whether to print log info during deconvolution.
#'   Errors will be printed regardless. 
#'
#' @return A List. Slot \strong{cors} contains list of vectors with correlation
#'   coefficients. Slot \strong{ctps} contains list of CTP objects returned by GetCTP
#'
GetCTPpcnorm <- function(bulk.eset, 
                         max_pc, 
                         cell_types, 
                         markers, 
                         ct_col, 
                         gene_col, 
                         min_gene, 
                         max_gene, 
                         weighted, 
                         w_col, 
                         verbose){
  marker_genes <- unique(markers[,gene_col])
  if (verbose) cat("Regressing out PCs from expression...\n")
  bulks <- NormalizeByPC(bulk.eset, max_pc=max_pc, genes_subset=marker_genes)
  if (verbose) cat("finished\n")
  ctps <- list()
  cors <- list()
  # names(ctps) <- names(cors) <- paste0("PC", as.character(0:max_pc))
  for (i in 0:max_pc){
    # pc_name <- paste0("PC", as.character(i))
    if (verbose) cat(paste0("Getting proportions from expression with ", as.character(i), " PCs out...\n"))
    ctp <- GetCTP(bulks[[as.character(i)]], cell_types, markers, ct_col, gene_col, min_gene, max_gene, weighted, w_col, verbose)
    if (verbose) cat("finished\n")
    ctp_pc1 = base::lapply(ctp, function(x) x$pcs[,1])
    ctp_pc1 = base::do.call(cbind, ctp_pc1)
    tc <- CorTri(ctp_pc1)
    ctps[[as.character(i)]] <- ctp
    cors[[as.character(i)]] <- tc
  }
  ret <- list(cors = cors, ctps = ctps)
  return(ret)
}

#' Performs reference-free deconvolution of bulk expression using marker genes
#' 
#' Estimates relative abundances of cell types from PCA-based deconvolution.
#' Uses a list of marker genes to subset the expression data, and returns the 
#' first PC of each sub-matrix as the cell type fraction estimates.
#' Optionally, weights for each marker gene can be used to prioritize genes
#' that are highly expressed in the given cell type.
#'
#' Note that this method expects the input bulk data to be normalized, unlike
#' the reference-based method.
#' 
#' @param bulk.eset Expression Set containing bulk data, 
#' @param markers Data frame with columns specifying cluster and gene,
#'   and optionally a column for weights, typically the fold-change of the gene.
#'   Important that the genes for each cell type are row-sorted by signficance.
#' @param counts Boolean. TRUE if raw counts are given.
#' @param normalized Boolean. TRUE if bulk data is to be treated as normalized. 
#' @param max_pc Numeric. If normalized=FALSE, the maximum number of PCs to consider
#'   regressing out for normalization.
#' @param ct.col Character string. Column name specifying cluster/cell type
#'   corresponding to each marker gene in \strong{markers}. 
#' @param gene.col Character string. Column name specifying gene names in
#'   \strong{markers}.
#' @param min.gene Numeric. Min number of genes to use for each cell type.
#' @param max.gene Numeric. Max number of genes to use for each cell type.
#' @param weighted Boolean. Whether to use weights for gene prioritization
#' @param w.col Character string. Column name for weights, such as "avg_logFC",
#'  in \strong{markers}
#' @param unique_markers Boolean. If TRUE, subset markers to include only genes 
#'   that are markers for only one cell type
#' @param pct_thresh Numeric. When normalized=FALSE. A fraction between 0 and 1 
#'   that specifies the maximum fraction of CTP-CTP correlation coefficients 
#'   that are positive. Expression data with all positive CTP-CTP correlations 
#'   are likely influenced by global covariation.
#' @param verbose Boolean. Whether to print log info during deconvolution.
#'   Errors will be printed regardless. 
#' @return A List. Slot \strong{bulk.props} contains estimated relative cell
#'   type abundances. Slot \strong{var.explained} contains variance explained by
#'   first 20 PCs for cell type marker genes. Slot \strong{genes.used} contains
#'   vector of genes used for deconvolution.
#' @export
ReferenceFreeDeconvolution <- function(bulk.eset, 
                                       markers,
                                       counts=TRUE,
                                       normalized=FALSE,
                                       max_pc=5,
                                       ct_col="cluster",
                                       gene_col="gene",
                                       min_gene = 5,
                                       max_gene = 200,
                                       weighted=FALSE,
                                       w_col = "avg_logFC",
                                       unique_markers = TRUE,
                                       pct_thresh = 0.5,
                                       verbose=TRUE){
  # Check input
  if ( ! methods::is(bulk.eset, "ExpressionSet") ) base::stop("Expression data should be in ExpressionSet")
  if (min_gene > max_gene){
    base::stop(base::paste0(base::sprintf("min_gene (set at %i) ", min_gene),
                            base::sprintf("must be less than or equal to max_gene (set at %i)\n", max_gene)))
  }

  # Get unique markers if applicable
  if (unique_markers){
    if (verbose) cat("Getting unique markers\n")
    markers <- GetUniqueMarkers(markers, gene_col=gene_col)
    n_genes_clust <- table(markers[,ct_col])
    if (any(n_genes_clust < min_gene)) {
      base::stop(paste0("Clusters must have a minimum of ", str(min_gene), " unique marker genes"))
    }
  }
  cg <- intersect(unique(markers[,gene_col]), rownames(bulk.eset))
  if (length(cg) == 0) base::stop(cat("No overlapping genes between markers and bulk.eset\n"))
  markers <- markers[markers[,gene_col] %in% cg,]

  # Set variables
  markers[,ct_col] <- as.character(markers[,ct_col])
  cell_types <- sort(unique(markers[,ct_col]))
  n_ct <- length(cell_types)
  n_s <- base::ncol(bulk.eset)
  if (verbose){
    base::cat(base::sprintf("Estimating proportions for %i cell types in %i samples in bulk\n", n_ct, n_s))
  }

  # Remove zero-variance genes, and un-expressed genes (if applicable)
  bulk.eset <- FilterZeroVarianceGenes(bulk.eset, verbose)
  bulk.eset <- FilterUnexpressedGenes(bulk.eset, verbose)

  # Get cell type proportions
  if (counts) bulk.eset <- CountsToCPM(bulk.eset)
  n_pc=NA
  if (!normalized){
    cors_ctps <- GetCTPpcnorm(bulk.eset, max_pc, cell_types, markers, ct_col, gene_col, min_gene, max_gene, weighted, w_col, verbose)
    # cors_ctps$cors has list of correlation coeffs and cors_ctps$ctps has list of ctp objects
    n_pc <- GetNumPC(cors_ctps[["cors"]], pct_thresh, verbose)
    ctp <- cors_ctps[["ctps"]][[n_pc]]
  } else {
    ctp <- GetCTP(bulk.eset, cell_types, markers, ct_col, gene_col, min_gene, max_gene, weighted, w_col, verbose)
  }

  names(ctp) = cell_types
  ctp_pc1 = base::lapply(ctp, function(x) x$pcs[,1])
  ctp_pc1 = base::do.call(cbind, ctp_pc1)
  ctp_pc1 = base::t(ctp_pc1)
  ctp_varexpl = base::sapply(ctp, function(x) x$sdev[1:20])
  rownames(ctp_varexpl) = base::paste0("PC", base::as.character(1:20))
  genes_used = base::lapply(ctp, function(x) x$genes)
  if (verbose) cat("Finished estimating cell type proportions using PCA\n")
  return(list(CTP=ctp_pc1, VarExpl=ctp_varexpl, markers=genes_used, n_pc=n_pc))
}

