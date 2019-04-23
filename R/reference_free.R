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
#' @return ret List. Attribute \strong{x} contains matrix of PCs, where PC1
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
  return(list( x = wpcs, sdev = sds, genes=colnames(x)))
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
#' @param bulk.eset Expression Set containing bulk data.
#'   Expression data should be normalized.
#' @param markers Data frame with columns specifying cluster and gene,
#'   and optionally a column for weights, typically the fold-change of the gene.
#'   Important that the genes in each row are sorted by signficance for each
#'   cell type.
#' @param ct.col Character string. Column name specifying cluster/cell type
#'   corresponding to each marker gene in \strong{markers}. 
#' @param gene.col Character string. Column name specifying gene names in
#'   \strong{markers}.
#' @param min.gene Numeric. Min number of genes to use for each cell type.
#' @param max.gene Numeric. Max number of genes to use for each cell type.
#' @param weighted Boolean. Whether to use weights for gene prioritization
#' @param w.col Character string. Column name for weights, such as "avg.logFC",
#'  in \strong{markers}
#' @param verbose Boolean. Whether to print log info during deconvolution.
#'   Errors will be printed regardless. 
#' @return A List. Slot \strong{bulk.props} contains estimated relative cell
#'   type abundances. Slot \strong{var.explained} contains variance explained by
#'   first 20 PCs for cell type marker genes. Slot \strong{genes.used} contains
#'   vector of genes used for deconvolution.
#' @export
ReferenceFreeDeconvolution <- function(bulk.eset, 
                                       markers, 
                                       ct.col="cluster", 
                                       gene.col="gene", 
                                       min.gene = 25, 
                                       max.gene = 200, 
                                       weighted=FALSE, 
                                       w.col = "avg.logFC", 
                                       verbose=TRUE) {
  if ( ! methods::is(bulk.eset, "ExpressionSet") ) {
    base::stop("Expression data should be in ExpressionSet")
  }
  if (min.gene > max.gene){
    base::stop(base::paste0(base::sprintf("min.gene (set at %i) ", min.gene),
                            "must be less than or equal to ",
                            base::sprintf("max.gene (set at %i)\n", max.gene)))
  }
  cell.types = sort(unique(markers[,ct.col]))
  n.ct = length(cell.types)
  n.s = base::ncol(bulk.eset)
  if (verbose){
    base::cat(base::sprintf("Estimating proportions for %i cell types", n.ct),
              base::sprintf("in %i samples in bulk\n", n.s), sep=" ")
  }
  # Remove unexpressed and zero-variance genes
  bulk.eset <- FilterZeroVarianceGenes(bulk.eset, verbose)
  bulk.eset <- FilterUnexpressedGenes(bulk.eset, verbose)
  ctp = base::lapply(cell.types,
           function(ct){
             markers.ct = markers[ markers[,ct.col] == ct , , drop=FALSE]
             ctm = base::make.unique(markers.ct[, gene.col])
             # Get markers in common between bulk and markers data frame
             common.markers = base::intersect(ctm,
                                              Biobase::featureNames(bulk.eset))
             if ( base::length(common.markers) == 0 ){
               base::stop("No marker genes found in bulk expression data")
             }
             expr = Biobase::exprs(bulk.eset)[common.markers,]
             expr = base::t(expr)
             if ( base::ncol(expr) < min.gene ){
               base::stop(base::paste0(base::sprintf("For cell type %s,", ct),
                                       " there are less marker genes in the ",
                                       base::sprintf("bulk expression set (%i)",
                                                     base::ncol(expr)),
                                       " than the minimum number of genes set",
                                       base::sprintf(" (%i) ", min.gene),
                                       "for PCA-based deconvolution\nSet the ",
                                       "min.gene parameter to a lower integer"))
             }
             if (weighted) {
                 # Get gene weights
                 ctw = markers.ct[, w.col]
                 base::names(ctw) = ctm
                 ctw = ctw[common.markers]
                 # Number of markers for PCA
                 ng = GetNumGenesWeighted(expr, ctw, min.gene, max.gene) 
                 expr = expr[,1:ng,drop=FALSE]
                 if (verbose) {
                   base::cat(base::sprintf("Using %i genes for cell type %s. ",
                                           ng, ct))
                 }
                 ret = EstimatePCACellTypeProportions(expr, weighted=TRUE,
                                                      w=ctw[1:ng])
             }
             else {
                 ng = GetNumGenes(expr, min.gene, max.gene)
                 expr[,1:ng,drop=FALSE]
                 if (verbose){
                   base::cat(base::sprintf("Using %i genes for cell type %s\n",
                                           ng, ct))
                 }
                 ret = EstimatePCACellTypeProportions(expr)
             }
             # Flip the sign of the first PC if negatively
             # correlated with most genes
             cors = stats::cor(expr, ret$x[,1])
             n.pos = sum(cors[,1] > 0)
             if (n.pos/base::length(cors[,1]) < (0.5)) {
               ret$x[,1] = ret$x[,1] * -1
             }
             if (verbose) {
               cors = stats::cor(expr, ret$x[,1])
               n.pos = base::sum(cors[,1] > 0)
               base::cat(base::sprintf("%i/%i", n.pos, base::length(cors[,1])),
                         "marker genes correlate positively",
                         base::sprintf("with PC1 for cell type %s\n", ct),
                         sep=" ")
             }
             return(ret)
           })
  base::names(ctp) = cell.types
  ctp.pc1 = base::lapply(ctp, function(x) x$x[,1])
  ctp.varexpl = base::sapply(ctp, function(x) x$sdev[1:20])
  base::rownames(ctp.varexpl) = base::paste0("PC", base::as.character(1:20))
  ctp.pc1 = base::do.call(cbind, ctp.pc1)
  ctp.pc1 = base::t(ctp.pc1)
  markers = base::lapply(ctp, function(x) x$genes)
  result <- list(bulk.props=ctp.pc1,
                 var.explained=ctp.varexpl,
                 genes.used=markers)
  return(result)
}
