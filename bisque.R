LoadSeuratObject <- function(obj, markers, delimiter="-", position=2) {
  # Extracts data from seurat object for deconvolution
  #
  # Note that only cells with labels are used.
  #
  # Args:
  #   obj: Seurat object with raw.data (dgTMatrix) and ident
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
  individual.ids <- sapply(strsplit(obj@cell.names, delimiter), `[[`, position)
  names(individual.ids) <- obj@cell.names
  individual.ids <- factor(individual.ids)
  sc <- list("data" = obj@raw.data[,names(obj@ident)],
             "markers" = unique(unlist(markers)),
             "cell.labels" = obj@ident,
             "individual.ids" = individual.ids)
  return(sc)
}

GetOverlappingSamples <- function(sc, bulk.data) {
  # Returns vector of samples found in both sc and bulk data
  bulk.samples <- colnames(bulk.data)
  sc.samples <- levels(sc@individual.ids)
  overlapping.samples <- intersect(bulk.samples, sc.samples)
  if (length(overlapping.samples) == 0) {
    stop("No overlapping samples found in bulk and single-cell")
  }
  return(overlapping.samples)
}

GetOverlappingGenes <- function(sc, bulk.data) {
  #Returns vector of genes found in markers, single-cell data, and bulk data
  bulk.genes <- row.names(bulk.data)
  sc.genes <- row.names(sc@data)
  overlapping.genes <- intersect(bulk.genes, intersect(sc@markers, sc.genes))
  if (length(overlapping.genes) == 0) {
    stop("No overlapping genes found between bulk and single-cell markers.")
  }
  return(overlapping.genes)
}

CountsToCPM <- function(counts, sparse=FALSE) {
  # Returns input matrix (sparse indicated by boolean argument) converted to CPM
  if (sparse) {
    counts@x <- (counts@x / colSums(counts)[counts@j + 1L]) * 1000000 
  }
  else {
    counts <- sweep(counts, 2, colSums(counts), `/`) * 1000000
  }
  return(counts)
}

GenerateSCReference <- function(sc) {
  # Averages expression within each cell type cluster, returning a n.genes
  # by n.celltype dense matrix
  Z <- vapply(levels(sc@cell.labels),
              function(x) { rowMeans(sc@data[,sc@cell.labels == x, drop=F]) },
              numeric(nrow(sc@data)))
  return(Z)
}

CalculateCellProportions <- function(sc) {
  # Counts cells belonging to each individual, returning a n.celltype by
  # n.individuals matrix of proportions
  sapply(levels(sc@individual.ids),
         function(x){
           table(sc@cell.labels[sc@individual.ids == x]) /
             length(sc@cell.labels[sc@individual.ids == x])
         })
}

TransformBulk <- function(gene, Y.train, X.train, X.pred) {
  # Learns linear transformation between observed bulk expression and linear
  # combination of sc-based reference and known proportions. Applies this
  # transformation to bulk samples to be deconvoluted. 
  #
  # Used with vapply, processes one gene
  Y.train.scaled <- scale(Y.train[gene,])
  Y.center <- attr(Y.train.scaled, "scaled:center")
  Y.scale <- attr(Y.train.scaled, "scaled:scale")
  X.train.scaled <- scale(X.train[gene,])
  X.center <- attr(X.train.scaled, "scaled:center")
  X.scale <- attr(X.train.scaled, "scaled:scale")
  X.pred.scaled <- scale(X.pred[gene,],
                         center=X.center,
                         scale=X.scale)
  coeff <- as.numeric(coefficients(lm(Y.train.scaled ~ X.train.scaled +0)))
  Y.pred.scaled <- X.scaled * coeff
  Y.pred <- matrix((Y.pred.scaled * Y.scale) + Y.center,
                   dimnames=list(colnames(X.pred), gene))
  return(Y.pred)
}

Deconvolute <- function(sc, bulk.data) {
  # Deconvolute bulk.data using the input sc class.
  overlapping.samples <- GetOverlappingSamples(sc, bulk.data)
  remaining.samples <- setdiff(colnames(bulk.data), overlapping.samples)
  overlapping.genes <- GetOverlappingGenes(sc, bulk.data)
  sc@data <- CountsToCPM(sc@data, sparse=TRUE)
  bulk.data <- CountsToCPM(bulk.data)
  ref.data <- GenerateSCReference(sc)[overlapping.genes,]
  sc.props <- CalculateCellProportions(sc)
  Y.train <- (ref.data %*% sc.props)[,overlapping.samples]
  X.train <- bulk.data[overlapping.genes,overlapping.samples]
  X.pred <- bulk.data[overlapping.genes,remaining.samples]
  template <- numeric(length(remaining.samples))
  names(template) <- remaining.samples
  Y.pred <- vapply(X=overlapping.genes, FUN=TransformBulk,
                   FUN.VALUE=template, Y.train, X.train, X.pred, USE.NAMES=TRUE)
  bulk.props <- apply(Y.pred, 1, function(b) {lsei::pnnls(ref.data, b, sum=1)}) 
  rownames(bulk.props) <- colnames(ref.data)
  colnames(bulk.props) <- remaining.samples
  return(bulk.props)
}
