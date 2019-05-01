#' Converts Seurat object to Expression Set
#'
#' `SeuratToExpressionSet()` returns an Expression Set with phenotype data
#' indicating cell type (cellType) and individual (SubjectName) for each cell
#' in a Seurat object. Raw counts data is used for assay data.
#'
#' Note that the \emph{Seurat} and \emph{Biobase} libraries should be attached
#' before running this function. The \emph{delimiter} and \emph{position} 
#' arguments are used to infer the individual ID from the cell ID. For example,
#' a delimiter of "-" and position of "2" indicates that the individual ID for
#' the cell ID \strong{ACTG-2} would be \strong{2}. 
#'
#' @param seurat.object Seurat object with attributes \emph{raw.data},
#'   \emph{ident}, and \emph{cell.names}
#' @param delimiter Character to split cell names with to find individual ID.
#' @param position Integer indicating 1-indexed position of individual ID after
#'   splitting cell name with \emph{delimiter}.
#' @param version Character string. Either "v2" or "v3. Seurat version used to
#'   create Seurat object. 
#' @return sc.eset Expression set containing relevant phenotype and individual
#'   data, \emph{cellType} and \emph{SubjectName}.
#'
#' @export
SeuratToExpressionSet <- function(seurat.object, delimiter, position,
                                  version = c("v2", "v3")) {
  if (! "Seurat" %in% base::.packages()) {
    base::stop("Seurat package is not attached")
  }
  version <- base::match.arg(version)
  if (version == "v2") {
    get.cell.names <- function(obj) obj@cell.names
    get.ident <- function(obj) obj@ident
    get.raw.data <- function(obj) obj@raw.data
  }
  else if (version == "v3") {
    get.cell.names <- function(obj) base::colnames(obj)
    get.ident <- function(obj) Seurat::Idents(object=obj)
    get.raw.data <- function(obj) Seurat::GetAssayData(object = obj,
                                                       slot = "counts")
  }
  individual.ids <- base::sapply(base::strsplit(get.cell.names(seurat.object),
                                                delimiter),
                                 `[[`, position)
  base::names(individual.ids) <- get.cell.names(seurat.object)
  individual.ids <- base::factor(individual.ids)
  n.individuals <- base::length(base::levels(individual.ids))
  base::cat(base::sprintf("Split sample names by \"%s\"", delimiter),
            base::sprintf("and checked position %i.", position),
            base::sprintf("Found %i individuals.\n", n.individuals),
            sep=" ")
  base::cat(base::sprintf("Example: \"%s\" corresponds to individual \"%s\".\n",
                          get.cell.names(seurat.object)[1], individual.ids[1]))
  sample.ids <- base::names(get.ident(seurat.object))
  sc.pheno <- base::data.frame(check.names=F, check.rows=F,
                               stringsAsFactors=F,
                               row.names=sample.ids,
                               SubjectName=individual.ids,
                               cellType=get.ident(seurat.object))
  sc.meta <- base::data.frame(labelDescription=base::c("SubjectName",
                                                       "cellType"),
                              row.names=base::c("SubjectName",
                                                "cellType"))
  sc.pdata <- methods::new("AnnotatedDataFrame",
                           data=sc.pheno,
                           varMetadata=sc.meta)
  sc.data <- base::as.matrix(get.raw.data(seurat.object)[,sample.ids])
  sc.eset <- Biobase::ExpressionSet(assayData=sc.data,
                                    phenoData=sc.pdata)
  return(sc.eset)
}

#' Convert counts data in Expression Set to counts per million (CPM)
#'
#' @param eset Expression Set containing counts assay data.
#' @return eset Expression Set containing CPM assay data
CountsToCPM <- function(eset) {
  Biobase::exprs(eset) <- base::sweep(Biobase::exprs(eset),
                                      2, base::colSums(Biobase::exprs(eset)),
                                      `/`) * 1000000
  indices <- base::apply(Biobase::exprs(eset), MARGIN=2,
                         FUN=function(column) {base::anyNA(column)})
  if (base::any(indices)) {
    n.cells <- base::sum(indices)
    base::stop(base::sprintf("Zero expression in selected genes for %i cells",
                             n.cells))
  }
  return(eset)
}

#' Remove genes in Expression Set with zero variance across samples
#'
#' @param eset Expression Set 
#' @param verbose Boolean. Print logging info
#' @return eset Expression Set with zero variance genes removed
FilterZeroVarianceGenes <- function(eset, verbose=TRUE) {
  indices <- (base::apply(Biobase::exprs(eset), 1, stats::var) != 0)
  indices <- indices & (! base::is.na(indices))
  if (base::sum(indices) > 0) {
    eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(eset)[indices,,drop=F],
                             phenoData=eset@phenoData)
  }
  if (verbose) {
    genes.filtered <- base::length(indices) - base::nrow(eset)
    base::cat(base::sprintf("Filtered %i zero variance genes.\n",
                            genes.filtered))
  }
  return(eset)
}

#' Remove genes in Expression Set with zero expression in all samples
#'
#' @param eset Expression Set
#' @param verbose Boolean. Print logging info
#' @return eset Expression Set with zero expression genes removed
FilterUnexpressedGenes <- function(eset, verbose=TRUE) {
  indices <- (base::apply(Biobase::exprs(eset), 1, base::sum) != 0)
  indices <- indices & (! base::is.na(indices))
  if (base::sum(indices) > 0) {
    eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(eset)[indices,,drop=F],
                             phenoData=eset@phenoData)
  }
  if (verbose) {
    genes.filtered <- base::length(indices) - base::nrow(eset)
    base::cat(base::sprintf("Filtered %i unexpressed genes.\n", genes.filtered))
  }
  return(eset)
}
