#TODO: Check that seurat object has ident, cell.names
#' @export
SeuratToExpressionSet <- function(seurat.object, delimiter="-", position=2) {
  # Extracts data from seurat object for deconvolution
  #
  # Note that only cells with labels are used.
  #
  # Args:
  #   seurat.object: Seurat object with raw.data (dgTMatrix) and ident
  #        (from Seurat::FindClusters)
  #   delimiter: Character to split cell.names with
  #   position: After splitting cell.name with delimiter, position of individual
  #             ID. 1-indexed.
  #
  # Returns:
  #   sc@counts: Raw count matrix
  #   sc@cells.labels: Cell type label for all cells
  #   sc@individual.ids: Sample label for all cells
  if (! "Seurat" %in% base::.packages()) {
    base::stop("Seurat package is not attached")
  }
  individual.ids <- base::sapply(base::strsplit(seurat.object@cell.names,
                                                delimiter),
                                 `[[`, position)
  base::names(individual.ids) <- seurat.object@cell.names
  individual.ids <- base::factor(individual.ids)
  n.individuals <- base::length(base::levels(individual.ids))
  base::cat(base::sprintf("Split sample names by \"%s\"", delimiter),
            base::sprintf("and checked position %i.", position),
            base::sprintf("Found %i individuals.\n", n.individuals),
            sep=" ")
  base::cat(base::sprintf("Example: \"%s\" corresponds to individual \"%s\".\n",
                          seurat.object@cell.names[1], individual.ids[1]))
  sample.ids <- base::names(seurat.object@ident)
  sc.pheno <- base::data.frame(check.names=F, check.rows=F,
                               stringsAsFactors=F,
                               row.names=sample.ids,
                               SubjectName=individual.ids,
                               cellType=seurat.object@ident)
  sc.meta <- base::data.frame(labelDescription=base::c("SubjectName",
                                                       "cellType"),
                              row.names=base::c("SubjectName",
                                                "cellType"))
  sc.pdata <- methods::new("AnnotatedDataFrame",
                           data=sc.pheno,
                           varMetadata=sc.meta)
  sc.data <- base::as.matrix(seurat.object@raw.data[,sample.ids])
  sc.eset <- Biobase::ExpressionSet(assayData=sc.data,
                                    phenoData=sc.pdata)
  return(sc.eset)
}

GetOverlappingSamples <- function(sc.eset, bulk.eset, subject.names, verbose) {
  # Returns vector of samples found in both sc and bulk data
  # subject.names is name of phenoData column for individual id 
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

GetOverlappingGenes <- function(sc.eset, bulk.eset, markers, verbose) {
  #Returns vector of genes found in single-cell and bulk data
  bulk.genes <- Biobase::featureNames(bulk.eset)
  sc.genes <- Biobase::featureNames(sc.eset)
  overlapping.genes <- base::intersect(bulk.genes, sc.genes)
  if (base::length(overlapping.genes) == 0) {
    base::stop(base::paste("No overlapping genes found between bulk and ",
                           "single-cell expression.", sep=""))
  }
  overlapping.genes <- base::intersect(overlapping.genes, markers)
  if (base::length(overlapping.genes) == 0) {
    base::stop(base::paste("No marker genes found in both bulk and ",
                           "single-cell expression.", sep=""))
  }
  if (verbose) {
    n.genes <- base::length(overlapping.genes)
    base::cat(base::sprintf("Using %i genes in both", n.genes),
                            "bulk and single-cell expression.\n",
                            sep=" ")
  }
  return(overlapping.genes)
}

CountsToCPM <- function(eset) {
  Biobase::exprs(eset) <- base::sweep(Biobase::exprs(eset),
                                      2, base::colSums(Biobase::exprs(eset)),
                                      `/`) * 1000000
  return(eset)
}

FilterZeroVarianceGenes <- function(eset, verbose=TRUE) {
  indices <- (base::apply(Biobase::exprs(eset), 1, stats::var) != 0)
  if (sum(indices) > 0) {
    eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(eset)[indices,,drop=F],
                             phenoData=eset@phenoData)
  }
  if (verbose) {
    genes.filtered <- base::length(indices) - base::nrow(eset)
    base::cat(sprintf("Filtered %i zero variance genes.\n", genes.filtered))
  }
  return(eset)
}

FilterUnexpressedGenes <- function(eset, verbose=TRUE) {
  indices <- (base::apply(Biobase::exprs(eset), 1, base::sum) != 0)
  if (sum(indices) > 0) {
    eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(eset)[indices,,drop=F],
                             phenoData=eset@phenoData)
  }
  if (verbose) {
    genes.filtered <- base::length(indices) - base::nrow(eset)
    base::cat(sprintf("Filtered %i unexpressed genes.\n", genes.filtered))
  }
  return(eset)
}

GenerateSCReference <- function(sc.eset, cell.types) {
  # Averages expression within each cell type cluster, returning a n.genes
  # by n.celltype dense matrix
  # cell.types is name of Phenodata column for cell types
  cell.labels <- base::factor(sc.eset[[cell.types]])
  all.cell.types <- base::levels(cell.labels)
  aggr.fn <- function(cell.type) {
    base::rowMeans(Biobase::exprs(sc.eset)[,cell.labels == cell.type, drop=F])
  }
  template <- base::numeric(base::nrow(sc.eset))
  sc.ref <- base::vapply(all.cell.types, aggr.fn, template)
  return(sc.ref)
}

CalculateSCCellProportions <- function(sc.eset, subject.names, cell.types) {
  # Counts cells belonging to each individual, returning a n.celltype by
  # n.individuals matrix of proportions
  # subject.names is name of phenoData column for individual id 
  individual.labels <- base::factor(sc.eset[[subject.names]])
  individuals <- base::levels(individual.labels)
  cell.labels <- sc.eset[[cell.types]]
  aggr.fn <- function(individual) {
     base::table(cell.labels[individual.labels == individual]) /
       base::length(cell.labels[individual.labels == individual])
  }
  sc.props <- base::sapply(individuals, aggr.fn)
  return(sc.props)
}

# NOTE
# I would put a check on the input that gene isn't a vector
# If gene is a vector, then base::scale will scale the columns of the Y.train and X.train matrices.
TransformBulk <- function(gene, Y.train, X.train, X.pred) {
  # Learns linear transformation between observed bulk expression and linear
  # combination of sc-based reference and known proportions. Applies this
  # transformation to bulk samples to be deconvoluted. 
  #
  # Used with vapply, processes one gene
  Y.train.scaled <- base::scale(Y.train[gene,,drop=T])
  Y.center <- base::attr(Y.train.scaled, "scaled:center")
  Y.scale <- base::attr(Y.train.scaled, "scaled:scale")
  X.train.scaled <- base::scale(X.train[gene,,drop=T])
  X.center <- base::attr(X.train.scaled, "scaled:center")
  X.scale <- base::attr(X.train.scaled, "scaled:scale")
  # If zero variance in both X and Y train, just solve coefficient directly
  # for one individual.
  if (base::anyNA(X.train.scaled) & base::anyNA(Y.train.scaled)) {
    coeff <- Y.train[gene,,drop=T][1]/X.train[gene,,drop=T][1]
    if (coeff == 0 || ! is.finite(coeff)) {
      coeff = NaN
    }
    Y.pred <- base::matrix(X.pred[gene,,drop=T] * coeff,
                           dimnames=base::list(base::colnames(X.pred), gene))
  }
  # If only one of X or Y has zero variance, return NaN. We shouldn't use this 
  # gene for deconvolution.
  else if (anyNA(X.train.scaled) || anyNA(Y.train.scaled)) {
    Y.pred <- base::matrix(X.pred[gene,,drop=T] * NaN,
                           dimnames=base::list(base::colnames(X.pred), gene))
  }
  # Otherwise, do standard linear model on scaled data, then unscale.
  else {
    X.pred.scaled <- base::scale(X.pred[gene,,drop=T],
                                 center=X.center,
                                 scale=X.scale)
    model <- stats::lm(Y.train.scaled ~ X.train.scaled +0)
    coeff <- base::as.numeric(stats::coefficients(model))
    Y.pred.scaled <- X.pred.scaled * coeff
    Y.pred <- base::matrix((Y.pred.scaled * Y.scale) + Y.center,
                           dimnames=base::list(base::colnames(X.pred), gene))
  }
  return(Y.pred)
}

#=====================================================================
# PCA + Weighted PCA functions
#=====================================================================
# Everything below here is what I added for the PCA-based deconvolution

#' Estimate cell type proportions using first PC of expression matrix
#' 
#' x is a sample x gene bulk expression matrix. Genes should be marker genes
#' If weighted = TRUE, multiply scaled gene expression by gene weights
#'
#' Returns a list:
#'	x matrix of PCs, first column contains PC1, which should be used as 
#'		as the estimate for cell type fractions
#'	sdev eigenvalues of eigendecomposition of var-covar matrix. The 
#'		1st eigenvalue should explain most of the variance.
pca_ctp <- function(x, weighted=FALSE, w=NULL){
	x <- base::scale(x)
	if (weighted){
		# Intersect gene names of weights and column names of x
		common_markers <- base::intersect( base::colnames(x), base::names(w) )
		if ( length(common_markers) == 0 ) base::stop("Genes from weights w do not match with column names of expression matrix x")
		x <- x[,common_markers]; w <- w[common_markers]
		wd <- base::diag(w)
		xw <- x %*% wd
		varcov <- t(xw) %*% xw
	}
	else varcov <- t(x) %*% x
	varcov_ed <- base::eigen(varcov)
	rot <- varcov_ed$vectors
	wpcs <- x %*% rot
	sds <- varcov_ed$values
	# x contains PCs, sdev contains eigenvalues of eigendecomposition
	return(list( x = wpcs, sdev = sds))
}

#' Get number of genes to use with weighted PCA
#' x is a sample by gene expression matrix containing the marker genes
#' w are the weights of the genes that correspond to the columns of x
#' min_gene and max_gene are the minimum and maximum number of genes to 
#'	consider as markers.
get_n_genes_w = function(x, w, min_gene = 25, max_gene = 200){
	max_gene = base::min(max_gene, base::ncol(x))
	ratios = base::rep(-Inf, max_gene)
	for (i in min_gene:max_gene){
		ret = pca_ctp(x[,1:i], weighted=TRUE, w=w[1:i])
		vars = ret$sdev
		ratios[i] = vars[1] / vars[2]
	}
	best_n = base::which.max(ratios)
	return(best_n)
}

#' Get number of genes to use, no weighted information
get_n_genes = function(x, min_gene = 25, max_gene = 200){
	max_gene = base::min(max_gene, base::ncol(x))
	ratios = base::rep(-Inf, max_gene)
	for (i in min_gene:max_gene){
		ret = pca_ctp(x[,1:i])
		vars = ret$sdev
		ratios[i] = vars[1] / vars[2]
	}
	best_n = base::which.max(ratios)
	return(best_n)
}	


#' Returns estimated "proportions" from PCA-based deconvolution
#' Uses a list of marker genes to subset the expression data, and returns the 
#' first PC as the cell type fraction estimates. Optionally, weights for each marker
#' can be used to prioritize genes that are highly expressed in the given cell type.
#' 
#' bulk.eset is ExpressionSet of bulk data (gene x sample). Expression data should be normalized
#' markers is a data frame with columns specifying cluster and gene, and optionally a 
#'	a column for weights, typically the fold-change of the gene. Important that the genes
#'	in each row are sorted by signficance for each cell type.
#' ct_col, gene_col, and w_col are the column names in markers for cell type, gene, and weight
#' 
#' Returns a list of 2 matrices. The first under the slot name CTP contains the estimated cell type proportions
#'	The slot VarExpl contains the variance explained by the first 20 PCs for the cell type marker bulk expression.
DeconvolutePCA <- function(bulk.eset, 
						   markers, 
						   ct_col="cluster", 
						   gene_col="gene", 
						   min_gene = 25, 
						   max_gene = 200, 
						   weighted=FALSE, 
						   w_col = "avg_logFC", 
						   verbose=TRUE){
	if ( ! methods::is(bulk.eset, "ExpressionSet") ) base::stop("Expression data should be in ExpressionSet")
	if (min_gene > max_gene){
		base::stop(base::paste0(base::sprintf("min_gene (set at %i) ", min_gene),
							 base::sprintf("must be less than or equal to max_gene (set at %i)\n", max_gene)))
	}
	cell_types = sort(unique(markers[,ct_col]))
	n_ct = length(cell_types)
	n_s = base::ncol(bulk.eset)
	if (verbose){
		base::cat(base::sprintf("Estimating proportions for %i cell types in %i samples in bulk\n", n_ct, n_s))
	}
	# Remove unexpressed and zero-variance genes
	bulk.eset <- FilterZeroVarianceGenes(bulk.eset, verbose)
	bulk.eset <- FilterUnexpressedGenes(bulk.eset, verbose)
	ctp = base::lapply(cell_types, function(ct){
					   markers_ct= markers[ markers[,ct_col] == ct , , drop=FALSE]
					   ctm = base::make.unique(markers_ct[, gene_col])
					   # Get markers in common between bulk and markers data frame
					   common_markers = base::intersect(ctm, Biobase::featureNames(bulk.eset))
					   if ( base::length(common_markers) == 0 ){
						   base::stop("No marker genes found in bulk expression data")
					   }
					   expr = Biobase::exprs(bulk.eset)[common_markers,]
					   expr = base::t(expr)
					   if ( base::ncol(expr) < min_gene ){
						   base::stop(base::paste0(base::sprintf("For cell type %s, There are less marker genes in ", ct), 
						   						   base::sprintf("the bulk expression set (%i) than the ", base::ncol(expr)),
						   						   base::sprintf("minimum number of genes set (%i) ", min_gene),
						   						   "for PCA-based deconvolution\nSet the min_gene parameter to a lower integer."))
					   }
					   if (weighted){
					   	   # Get gene weights
					   	   ctw = markers_ct[, w_col]; names(ctw) = ctm; ctw = ctw[common_markers]
					   	   ng = get_n_genes_w(expr, ctw, min_gene, max_gene) # Number of markers for PCA
					   	   if (verbose){
							   base::cat(base::sprintf("Using %i genes for cell type %s\n", ng, ct))
						   }
					   	   ret = pca_ctp(expr[,1:ng,drop=FALSE], weighted=TRUE, w=ctw[1:ng])
					   }
					   else{
					   	   ng = get_n_genes(expr, min_gene, max_gene)
					   	   if (verbose){
							   base::cat(base::sprintf("Using %i genes for cell type %s\n", ng, ct))
						   }
					   	   ret = pca_ctp(expr)
					   }
					   return(ret)
						   })
	names(ctp) = cell_types
	ctp_pc1 = base::lapply(ctp, function(x) x$x[,1])
	ctp_varexpl = base::sapply(ctp, function(x) x$sdev[1:20])
	rownames(ctp_varexpl) = base::paste0("PC", base::as.character(1:20))
	ctp_pc1 = base::do.call(cbind, ctp_pc1)
	return(list(CTP=ctp_pc1, VarExpl=ctp_varexpl))
}



#=====================================================================
#=====================================================================


# TODO: Catch division by 0 in CPM conversion (this happens if sample has no expression)
#' @export
Deconvolute <- function(sc.eset, bulk.eset, markers=NULL, cell.types="cellType",
                        subject.names="SubjectName", verbose=TRUE) {
  if ((! methods::is(sc.eset, "ExpressionSet")) || 
      (! methods::is(bulk.eset, "ExpressionSet"))) {
    base::stop("Expression data should be in ExpressionSet")
  }
  else if (! cell.types %in% Biobase::varLabels(sc.eset)) {
    base::stop(base::cat(base::sprintf("Cell type label \"%s\"", cell.types),
                         "not found in single-cell ExpressionSet varLabels.",
                         sep=" "))
  }
  else if (! subject.names %in% Biobase::varLabels(sc.eset)) {
    base::stop(base::cat(base::sprintf("Individual label \"%s\"",
                                       subject.names),
                         "not found in single-cell ExpressionSet varLabels.",
                         sep=" "))
  }
  n.cell.types <-
    base::length(base::levels(base::factor(sc.eset[[cell.types]])))
  if (n.cell.types == 1) {
    base::stop(base::cat("Single-cell pheno data indicates only one cell type",
                         "present. No need to deconvolute.", sep=" "))
  }
  if (verbose) {
    base::cat(base::sprintf("Deconvoluting into %i cell types.\n",
                              n.cell.types))
  }
  samples <- GetOverlappingSamples(sc.eset, bulk.eset, subject.names, verbose)
  if (base::is.null(markers)) {
    markers <- Biobase::featureNames(sc.eset)
  }
  else {
    markers <- base::unique(base::unlist(markers))
  }
  genes <- GetOverlappingGenes(sc.eset, bulk.eset, markers, verbose)
  sc.eset <-
    Biobase::ExpressionSet(assayData=Biobase::exprs(sc.eset)[genes,],
                           phenoData=sc.eset@phenoData)
  bulk.eset <-
    Biobase::ExpressionSet(assayData=Biobase::exprs(bulk.eset)[genes,],
                           phenoData=bulk.eset@phenoData)
  if (verbose) {
    base::cat("Converting single-cell counts to CPM and",
              "filtering zero variance genes.\n", sep=" ")
  }
  sc.eset <- CountsToCPM(sc.eset)
  sc.eset <- FilterZeroVarianceGenes(sc.eset, verbose)
  if (verbose) {
    base::cat("Converting bulk counts to CPM and filtering",
              "unexpressed genes.\n", sep=" ")
  }
  bulk.eset <- CountsToCPM(bulk.eset)
  bulk.eset <- FilterUnexpressedGenes(bulk.eset, verbose)
  genes <- base::intersect(Biobase::featureNames(sc.eset),
                           Biobase::featureNames(bulk.eset))
  if (base::length(genes) == 0) {
    base::stop(base::cat("Zero genes remaining after filtering and",
                         "intersecting bulk, single-cell, and marker genes.",
                         sep=" "))
  }
  if (verbose) {
    n.cells <- base::ncol(sc.eset)
    base::cat("Generating single-cell based reference from",
              sprintf("%i cells.\n", n.cells), sep=" ")
  }
  sc.ref <- GenerateSCReference(sc.eset, cell.types)[genes,,drop=F]
  sc.props <- CalculateSCCellProportions(sc.eset, subject.names, cell.types)
  sc.props <- sc.props[base::colnames(sc.ref),,drop=F]
  if (verbose) {
    base::cat("Learning bulk transformation.\n")
  }
  # Y.train is pseudo-bulk expression based on reference profile weighted by
  #   cell type proportions estimated for single-cell samples.
  Y.train <- sc.ref %*% sc.props[,samples$overlapping,drop=F]
  # X.train is the actual bulk for the single-cell samples.
  X.train <- Biobase::exprs(bulk.eset)[genes,samples$overlapping,drop=F]
  # X.pred is the bulk for the remaining samples to be deconvoluted.
  X.pred <- Biobase::exprs(bulk.eset)[genes,samples$remaining,drop=F]
  template <- base::numeric(base::length(samples$remaining))
  base::names(template) <- samples$remaining
  if (verbose) {
    base::cat("Applying transformation to bulk samples and deconvoluting.\n")
  }
  # Y.pred is the transformed bulk for samples to be deconvoluted.
  Y.pred <- base::matrix(base::vapply(X=genes, FUN=TransformBulk,
                                      FUN.VALUE=template,
                                      Y.train, X.train, X.pred,
                                      USE.NAMES=TRUE),
                         nrow=base::length(samples$remaining))
  # Columns in Y.pred with NaN indicate transformation could not be learned
  #   for that gene. 
  indices <- base::apply(Y.pred, MARGIN=2,
                         FUN=function(column) {base::anyNA(column)})
  if (base::any(indices)) {
    if (verbose) {
      n.dropped <- base::sum(indices)
      base::cat(base::sprintf("Dropped an additional %i genes", n.dropped),
                "for which a transformation could not be learned.\n", sep=" ")
    }
    if (sum(!indices) == 0) {
      base::stop("Zero genes left for deconvolution.")
    }
    Y.pred <- Y.pred[,!indices,drop=F]
    sc.ref <- sc.ref[!indices,,drop=F]
  }
  bulk.props <- base::as.matrix(base::apply(Y.pred, 1,
                                            function(b) {
                                              lsei::pnnls(sc.ref, b, sum=1)$x
                                            }))
  base::rownames(bulk.props) <- base::colnames(sc.ref)
  base::colnames(bulk.props) <- samples$remaining
  #props <- base::cbind(sc.props, bulk.props)
  #return(props)
  return(bulk.props)
}
