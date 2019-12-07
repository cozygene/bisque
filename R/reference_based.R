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
    base::message(base::sprintf("Found %i samples ", n.overlapping),
                                "with bulk and single-cell expression.")
    base::message(base::sprintf("Remaining %i ", n.remaining),
                                "bulk samples will be decomposed.")
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
    base::message(base::sprintf("Using %i genes in both", n.genes),
                                " bulk and single-cell expression.")
  }
  return(overlapping.genes)
}

#' Generate reference profile for cell types identified in single-cell data
#'
#' Averages expression within each cell type across all samples to use as 
#' reference profile.
#'
#' @param sc.eset Expression Set with single-cell data
#' @param cell.types A character string. Name of phenoData attribute in sc.eset
#'   that indicates cell type
#' @return sc.ref Matrix. Reference profile with number of gene rows by number
#'   of cell types columns. 
GenerateSCReference <- function(sc.eset, cell.types) {
  cell.labels <- base::factor(sc.eset[[cell.types]])
  all.cell.types <- base::levels(cell.labels)
  aggr.fn <- function(cell.type) {
    base::rowMeans(Biobase::exprs(sc.eset)[,cell.labels == cell.type, drop=F])
  }
  template <- base::numeric(base::nrow(sc.eset))
  sc.ref <- base::vapply(all.cell.types, aggr.fn, template)
  return(sc.ref)
}

#' Calculate cell proportions based on single-cell data
#'
#' Returns proportion of each cell type out of total cells for each individual
#' in the single-cell Expression Set
#'
#' @param sc.eset Expression Set with single-cell data
#' @param subject.names A character string. Name of phenoData attribute in 
#'   sc.eset that indicates individual ID.
#' @param cell.types A character string. Name of phenoData attribute in sc.eset
#'   that indicates cell type
#' @return sc.props Matrix. Cell proportions with number of cell types rows 
#'   by number of individuals columns
CalculateSCCellProportions <- function(sc.eset, subject.names, cell.types) {
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

#' Transforms bulk expression of a gene given overlapping data
#'
#' For a specific gene, this function uses linear regression to learn a
#' transformation of the bulk expression to match the values produced
#' by the single-cell based reference and observed single-cell based cell
#' proportions.
#'
#' If a linear transformation cannot be learned for a gene (zero variance in
#' observed bulk or single-cell based weighted sums), a vector of NaNs will 
#' be returned of the expected length (length of X.pred)
#'
#' @param gene Character string. Gene name that corresponds to row in Y.train
#' @param Y.train Numeric Matrix. Number of gene rows by number of overlapping
#'   individuals columns. Contains weighted sum of reference profile by 
#'   single-cell based cell proportion estimates for each individual
#' @param X.train Numeric Matrix. Number of gene rows by number of overlapping
#'   individuals columns. Contains observed bulk expression for each individual
#' @param X.pred Numeric Matrix. Number of gene rows by number of remaining
#'   individuals columns. Contains observed bulk expression for each individual
#'   to be transformed.
#' @return Y.pred Numeric Matrix. One row for given gene by number of remaining
#'   individuals columns. Contains transformed bulk expression for each
#'   individual.
SupervisedTransformBulk <- function(gene, Y.train, X.train, X.pred) {
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
  # gene for decomposition.
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

#' Transforms bulk expression of a gene using only single-cell data
#'
#' For a specific gene, this function learns a transformation of
#' the bulk expression to match the distribution produced
#' by the single-cell based reference and observed single-cell based cell
#' proportions.
#'
#' @param gene Character string. Gene name that corresponds to row in Y.train
#' @param Y.train Numeric Matrix. Number of gene rows by number of overlapping
#'   individuals columns. Contains weighted sum of reference profile by 
#'   single-cell based cell proportion estimates for each individual
#' @param X.pred Numeric Matrix. Number of gene rows by number of remaining
#'   individuals columns. Contains observed bulk expression for each individual
#'   to be transformed.
#' @return Y.pred Numeric Matrix. One row for given gene by number of remaining
#'   individuals columns. Contains transformed bulk expression for each
#'   individual.
SemisupervisedTransformBulk <- function(gene, Y.train, X.pred) {
  # Learns linear transformation of observed bulk to match distribution of
  # weighted sum of reference
  #
  # Used with vapply, processes one gene
  Y.train.scaled <- base::scale(Y.train[gene,,drop=T])
  Y.center <- base::attr(Y.train.scaled, "scaled:center")
  Y.scale <- base::attr(Y.train.scaled, "scaled:scale")
  n <- base::length(Y.train.scaled)
  # Shrinkage estimator that minimizes MSE for scaling factor
  shrink.scale <- base::sqrt(base::sum((Y.train[gene,,drop=T]-Y.center)^2)/n+1)
  X.pred.scaled <- base::scale(X.pred[gene,,drop=T])
  Y.pred <- base::matrix((X.pred.scaled * shrink.scale) + Y.center,
                         dimnames=base::list(base::colnames(X.pred), gene))
  return(Y.pred)
}

#' Performs reference-based decomposition of bulk expression using single-cell
#' data
#'
#' Generates a reference profile based on single-cell data. Learns a
#' transformation of bulk expression based on observed single-cell proportions
#' and performs  NNLS regression on these transformed values to estimate cell
#' proportions.
#'
#' Expects read counts for both datasets, as they will be converted to 
#' counts per million (CPM). Two options available: Use overlapping indivudals
#' found in both single-cell and bulk datasets to learn transformation or 
#' learn transformation from single-cell alone. The overlapping option is
#' expected to have better performance. 
#'
#' @param bulk.eset Expression Set containin bulk data. No PhenoData required
#'   but if overlapping option used, IDs returned by sampleNames(bulk.eset) 
#'   should match those found in sc.eset phenoData individual labels.
#' @param sc.eset Expression Set containing single-cell data. PhenoData of this
#'   Expression Set should contain cell type and individual labels for each
#'   cell. Names of these fields specified by arguments below.
#' @param markers Structure, such as character vector, containing marker genes
#'   to be used in decomposition. `base::unique(base::unlist(markers))` should
#'   return a simple vector containing each gene name. If no argument or NULL
#'   provided, the method will use all available genes for decomposition.
#' @param cell.types Character string. Name of phenoData attribute in sc.eset
#'   indicating cell type label for each cell
#' @param subject.names Character string. Name of phenoData attribute in sc.eset
#'   indicating individual label for each cell
#' @param use.overlap Boolean. Whether to use and expect overlapping samples 
#'   in decomposition.
#' @param verbose Boolean. Whether to print log info during decomposition.
#'   Errors will be printed regardless. 
#' @return A list. Slot \strong{bulk.props} contains a matrix of cell type
#'   proportion estimates with cell types as rows and individuals as columns.
#'   Slot \strong{sc.props} contains a matrix of cell type proportions 
#'   estimated directly from counting single-cell data. 
#'   Slot \strong{rnorm} contains Euclidean norm of the residuals for each
#'   individual's proportion estimates. Slot \strong{genes.used} contains
#'   vector of genes used in decomposition
#' @examples
#' library(Biobase)
#' sim.data <- SimulateData(n.ind=10, n.genes=100, n.cells=100,
#'                          cell.types=c("Neurons", "Astrocytes", "Microglia"),
#'                          avg.props=c(.5, .3, .2))
#' sim.data$sc.eset <- sim.data$sc.eset[,sim.data$sc.eset$SubjectName %in% as.character(6:10)]
#' res <- ReferenceBasedDecomposition(sim.data$bulk.eset, sim.data$sc.eset)
#' estimated.cell.proportions <- res$bulk.props
#' 
#' @export
ReferenceBasedDecomposition <- function(bulk.eset,
                                        sc.eset,
                                        markers=NULL,
                                        cell.types="cellType",
                                        subject.names="SubjectName",
                                        use.overlap=TRUE, 
                                        verbose=TRUE) {
  if ((! methods::is(sc.eset, "ExpressionSet")) || 
      (! methods::is(bulk.eset, "ExpressionSet"))) {
    base::stop("Expression data should be in ExpressionSet")
  }
  else if (! cell.types %in% Biobase::varLabels(sc.eset)) {
    base::stop(base::sprintf("Cell type label \"%s\" ", cell.types),
               "not found in single-cell ExpressionSet varLabels.")
  }
  else if (! subject.names %in% Biobase::varLabels(sc.eset)) {
    base::stop(base::sprintf("Individual label \"%s\"", subject.names),
               " not found in single-cell ExpressionSet varLabels.")
  }
  n.sc.individuals <-
    base::length(base::levels(base::factor(sc.eset[[subject.names]])))
  if (n.sc.individuals == 1) {
    base::stop("Only one individual detected in single-cell data. At least ",
               "two subjects are needed (three or more recommended).")
  }
  else if (n.sc.individuals == 2) {
    base::warning("Only two individuals detected in single-cell data. While ",
                  "Bisque will run, we recommend at least three subjects for",
                  " reliable performance.")
  }
  n.cell.types <-
    base::length(base::levels(base::factor(sc.eset[[cell.types]])))
  if (n.cell.types == 1) {
    base::stop("Single-cell pheno data indicates only one cell type",
               " present. No need for decomposition.")
  }
  if (verbose) {
    base::message(base::sprintf("Decomposing into %i cell types.",
                                n.cell.types))
  }
  if (use.overlap) {
    samples <- GetOverlappingSamples(sc.eset, bulk.eset, subject.names, verbose)
  }
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
    base::message("Converting single-cell counts to CPM and ",
              "filtering zero variance genes.")
  }
  sc.eset <- CountsToCPM(sc.eset)
  sc.eset <- FilterZeroVarianceGenes(sc.eset, verbose)
  if (verbose) {
    base::message("Converting bulk counts to CPM and filtering",
                  " unexpressed genes.")
  }
  bulk.eset <- CountsToCPM(bulk.eset)
  bulk.eset <- FilterUnexpressedGenes(bulk.eset, verbose)
  genes <- base::intersect(Biobase::featureNames(sc.eset),
                           Biobase::featureNames(bulk.eset))
  if (base::length(genes) == 0) {
    base::stop("Zero genes remaining after filtering and ",
               "intersecting bulk, single-cell, and marker genes.")
  }
  if (verbose) {
    n.cells <- base::ncol(sc.eset)
    base::message("Generating single-cell based reference from ",
                  sprintf("%i cells.\n", n.cells))
  }
  sc.ref <- GenerateSCReference(sc.eset, cell.types)[genes,,drop=F]
  sc.props <- CalculateSCCellProportions(sc.eset, subject.names, cell.types)
  sc.props <- sc.props[base::colnames(sc.ref),,drop=F]
  if (use.overlap) {
    if (verbose) {
      base::message("Learning bulk transformation from overlapping samples.")
    }
    # Y.train is pseudo-bulk expression based on reference profile weighted by
    #   cell type proportions estimated for single-cell samples.
    Y.train <- sc.ref %*% sc.props[,samples$overlapping,drop=F]
    # X.train is the actual bulk for the single-cell samples.
    X.train <- Biobase::exprs(bulk.eset)[genes,samples$overlapping,drop=F]
    # X.pred is the bulk for the remaining samples to be decomposed.
    X.pred <- Biobase::exprs(bulk.eset)[genes,samples$remaining,drop=F]
    template <- base::numeric(base::length(samples$remaining))
    base::names(template) <- samples$remaining
    if (verbose) {
      base::message("Applying transformation to bulk samples and decomposing.")
    }
    # Y.pred is the transformed bulk for samples to be decomposed.
    Y.pred <- base::matrix(base::vapply(X=genes, FUN=SupervisedTransformBulk,
                                        FUN.VALUE=template,
                                        Y.train, X.train, X.pred,
                                        USE.NAMES=TRUE),
                           nrow=base::length(samples$remaining))
    sample.names <- samples$remaining
  }
  else {
    if (verbose) {
      base::message("Inferring bulk transformation from single-cell alone.")
    }
    # Y.train is pseudo-bulk expression based on reference profile weighted by
    #   cell type proportions estimated for single-cell samples.
    Y.train <- sc.ref %*% sc.props
    # X.pred is the bulk for the remaining samples to be decomposed.
    X.pred <- Biobase::exprs(bulk.eset)[genes,,drop=F]
    sample.names <- base::colnames(Biobase::exprs(bulk.eset))
    template <- base::numeric(base::length(sample.names))
    base::names(template) <- sample.names
    if (verbose) {
      base::message("Applying transformation to bulk samples and decomposing.")
    }
    # Y.pred is the transformed bulk for samples to be decomposed.
    Y.pred <- base::matrix(base::vapply(X=genes,
                                        FUN=SemisupervisedTransformBulk,
                                        FUN.VALUE=template,
                                        Y.train, X.pred,
                                        USE.NAMES=TRUE),
                           nrow=base::length(sample.names))
  }
  # Columns in Y.pred with NaN indicate transformation could not be learned
  #   for that gene. 
  indices <- base::apply(Y.pred, MARGIN=2,
                         FUN=function(column) {base::anyNA(column)})
  if (base::any(indices)) {
    if (verbose) {
      n.dropped <- base::sum(indices)
      base::message(base::sprintf("Dropped an additional %i genes", n.dropped),
                    " for which a transformation could not be learned.")
    }
    if (sum(!indices) == 0) {
      base::stop("Zero genes left for decomposition.")
    }
    Y.pred <- Y.pred[,!indices,drop=F]
    sc.ref <- sc.ref[!indices,,drop=F]
  }
  results <- base::as.matrix(base::apply(Y.pred, 1,
                                         function(b) {
                                           sol <- lsei::pnnls(sc.ref,
                                                              b, sum=1)
                                           return(base::append(sol$x,
                                                               sol$rnorm))
                                         }))
  base::rownames(results) <- base::append(base::colnames(sc.ref), "rnorm")
  base::colnames(results) <- sample.names
  rnorm <- results["rnorm",,drop=T]
  names(rnorm) <- sample.names
  results <- base::list(bulk.props=results[base::colnames(sc.ref),,drop=F],
                        sc.props=sc.props,
                        rnorm=rnorm,
                        genes.used=base::rownames(sc.ref))
  return(results)
}
