bisque.sc <- setClass(
  "bisque.sc",
  slots=c(
    data = "ANY",
    cell.labels = "factor",
    individual.ids = "factor"
  )
)


#' @export
`[.bisque.sc` <- function(x, i){
  index <- x@individual.ids %in% i
  bisque.sc.object <- CreateBisqueSCObject(x@data[,index],
                                           droplevels(x@cell.labels[index]),
                                           droplevels(x@individual.ids[index]))
  return(bisque.sc.object)
}

#' @export
CreateBisqueSCObject <- function(sc.data, cell.type.factor, individual.factor){
  bisque.sc.object <- new(Class = "bisque.sc",
                          data = sc.data,
                          cell.labels = cell.type.factor,
                          individual.ids = individual.factor)
  return(bisque.sc.object)
}

#' @export
LoadFromCSV <- function(sc.csv, cell.csv, markers.txt, delimiter="-",
                        position=2, sparse=TRUE) {
  # Load single-cell info from csv files for deconvolution
  #
  # Args:
  #   sc.csv: Path to csv with single-cell data (genes as rows)
  #   cell.csv: Path to csv with cell ID as first column and cell-type as
  #             second column
  #   markers.txt: Path to file with 1 marker gene per line, no header
  sc.data <- as.matrix(read.csv(sc.csv, header=TRUE,
                                check.names=FALSE, row.names=1))
  if (sparse){
    sc.col.names <- colnames(sc.data)
    sc.row.names <- rownames(sc.data)
    sc.data <- as(sc.data, "dgTMatrix")
    colnames(sc.data) <- sc.col.names
    rownames(sc.data) <- sc.row.names
  }
  cell.table <- read.csv(cell.csv, header=TRUE, check.names=FALSE, row.names=1)
  cell.labels <- factor(cell.table[,1])
  names(cell.labels) <- rownames(cell.table)
  individual.ids <- sapply(base::strsplit(names(cell.labels), delimiter),
                           `[[`, position)
  names(individual.ids) <- names(cell.labels)
  markers <- base::unique(base::scan(markers.txt, what = character()))
  genes <- base::intersect(rownames(sc.data), markers)
  bisque.sc.object <- new(Class="bisque.sc",
                          data = sc.data[genes,names(cell.labels)],
                          cell.labels = cell.labels,
                          individual.ids = individual.ids)
  return(bisque.sc.object)
}

#' @export
LoadFromSeurat <- function(seurat.object, markers, delimiter="-", position=2) {
  # Extracts data from seurat object for deconvolution
  #
  # Note that only cells with labels are used.
  #
  # Args:
  #   seurat.object: Seurat object with raw.data (dgTMatrix) and ident
  #        (from Seurat::FindClusters)
  #   markers: Object returned by Seurat::FindAllMarkers
  #   delimiter: Character to split cell.names with
  #   position: After splitting cell.name with delimiter, position of individual
  #             ID. 1-indexed.
  #
  # Returns:
  #   sc@counts: Raw count matrix
  #   sc@markers: Marker genes for all cell types
  #   sc@cells.labels: Cell type label for all cells
  #   sc@individual.ids: Sample label for all cells
  individual.ids <- sapply(base::strsplit(seurat.object@cell.names, delimiter),
                           `[[`, position)
  names(individual.ids) <- seurat.object@cell.names
  individual.ids <- factor(individual.ids)
  genes <- base::intersect(rownames(seurat.object@raw.data),
                           base::unique(unlist(markers)))
  bisque.sc.object <- new(Class="bisque.sc",
                          data = seurat.object@raw.data[genes,names(seurat.object@ident)],
                          cell.labels = seurat.object@ident,
                          individual.ids = individual.ids)
  return(bisque.sc.object)
}

GetOverlappingSamples <- function(sc, bulk.data) {
  # Returns vector of samples found in both sc and bulk data
  bulk.samples <- colnames(bulk.data)
  sc.samples <- levels(sc@individual.ids)
  overlapping.samples <- base::intersect(bulk.samples, sc.samples)
  if (length(overlapping.samples) == 0) {
    stop("No overlapping samples found in bulk and single-cell")
  }
  return(overlapping.samples)
}

GetOverlappingGenes <- function(sc, bulk.data) {
  #Returns vector of genes found in single-cell and bulk data
  bulk.genes <- rownames(bulk.data)
  sc.genes <- rownames(sc@data)
  overlapping.genes <- base::intersect(bulk.genes, sc.genes)
  if (length(overlapping.genes) == 0) {
    stop("No overlapping genes found between bulk and single-cell markers.")
  }
  return(overlapping.genes)
}

CountsToCPM <- function(counts) {
  # Returns input matrix (dgTMatrix or base::matrix) converted to CPM
  count.class <- class(counts)
  if (count.class == "dgTMatrix") {
    counts@x <- (counts@x / Matrix::colSums(counts)[counts@j + 1L]) * 1000000 
  }
  else if (count.class == "matrix") {
    counts <- sweep(counts, 2, base::colSums(counts), `/`) * 1000000
  }
  else {
    stop("Input matrix not base::matrix or dgTMatrix")
  }
  return(as(counts, count.class))
}

GenerateSCReference <- function(sc) {
  # Averages expression within each cell type cluster, returning a n.genes
  # by n.celltype dense matrix
  sc.class <- class(sc@data)
  if (sc.class == "dgTMatrix") {
    Z <- vapply(levels(sc@cell.labels),
                function(x) { Matrix::rowMeans(sc@data[,sc@cell.labels == x,
                                                       drop=F]) },
                numeric(nrow(sc@data)))
  }
  else if (sc.class == "matrix") {
    Z <- vapply(levels(sc@cell.labels),
                function(x) { base::rowMeans(sc@data[,sc@cell.labels == x,
                                                     drop=F]) },
                numeric(nrow(sc@data)))
  }
  else {
    stop("Input matrix not base::matrix or dgTMatrix")
  }
  return(Z)
}

#' @export
CalculateCellProportions <- function(sc) {
  # Counts cells belonging to each individual, returning a n.celltype by
  # n.individuals matrix of proportions
  sapply(levels(sc@individual.ids),
         function(x){
           base::table(sc@cell.labels[sc@individual.ids == x]) /
             length(sc@cell.labels[sc@individual.ids == x])
         })
}

TransformBulk <- function(gene, Y.train, X.train, X.pred) {
  # Learns linear transformation between observed bulk expression and linear
  # combination of sc-based reference and known proportions. Applies this
  # transformation to bulk samples to be deconvoluted. 
  #
  # Used with vapply, processes one gene
  Y.train.scaled <- base::scale(Y.train[gene,])
  Y.center <- attr(Y.train.scaled, "scaled:center")
  Y.scale <- attr(Y.train.scaled, "scaled:scale")
  X.train.scaled <- base::scale(X.train[gene,])
  X.center <- attr(X.train.scaled, "scaled:center")
  X.scale <- attr(X.train.scaled, "scaled:scale")
  X.pred.scaled <- base::scale(X.pred[gene,],
                               center=X.center,
                               scale=X.scale)
  coeff <- as.numeric(coefficients(stats::lm(Y.train.scaled ~ X.train.scaled +0)))
  Y.pred.scaled <- X.pred.scaled * coeff
  Y.pred <- matrix((Y.pred.scaled * Y.scale) + Y.center,
                   dimnames=list(colnames(X.pred), gene))
  return(Y.pred)
}

#' @export
Deconvolute <- function(sc, bulk.data) {
  # Deconvolute bulk.data using the input sc class.
  if (class(sc) != "bisque.sc") {
    stop("Single-cell argument should be bisque.sc object")
  }
  overlapping.samples <- GetOverlappingSamples(sc, bulk.data)
  remaining.samples <- base::setdiff(colnames(bulk.data), overlapping.samples)
  # Remove 0 variance genes from bulk
  # Ideally do this to single-cell too, but marker genes chosen from single cell
  # so these genes should have nonzero variance. 
  bulk.data <- bulk.data[apply(bulk.data, 1, stats::var) != 0,]
  overlapping.genes <- GetOverlappingGenes(sc, bulk.data)
  sc@data <- CountsToCPM(sc@data)
  bulk.data <- CountsToCPM(bulk.data)
  ref.data <- GenerateSCReference(sc)[overlapping.genes,,drop=F]
  sc.props <- CalculateCellProportions(sc)
  Y.train <- ref.data %*% sc.props[,overlapping.samples,drop=F]
  X.train <- bulk.data[overlapping.genes,overlapping.samples,drop=F]
  X.pred <- bulk.data[overlapping.genes,remaining.samples,drop=F]
  template <- numeric(length(remaining.samples))
  names(template) <- remaining.samples
  Y.pred <- matrix(vapply(X=overlapping.genes, FUN=TransformBulk,
                          FUN.VALUE=template, Y.train, X.train, X.pred,
                          USE.NAMES=TRUE),
                   nrow=length(remaining.samples))
  bulk.props <- as.matrix(apply(Y.pred, 1, function(b) {
                                             lsei::pnnls(ref.data, b, sum=1)$x
                                           }))
  rownames(bulk.props) <- colnames(ref.data)
  colnames(bulk.props) <- remaining.samples
  return(bulk.props)
}
