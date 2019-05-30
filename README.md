# Bisque

An R toolkit for accurate and efficient estimation of cell composition ('decomposition') from bulk expression data with single-cell information.

Bisque provides two modes of operation:

### Reference-based decomposition
This method utilizes single-cell data to decompose bulk expression.
We assume that both single-cell and bulk counts are measured from the same tissue.
Specifically, the cell composition of the labeled single-cell data should match the expected physiological composition.
While we don't explicitly require matched samples, we expect having samples with both single-cell and bulk expression measured will provide more accurate results.

### Marker-based decomposition
This method utilizes marker genes alone to decompose bulk expression when a reference profile is not available.
Single-cell data is not explicitly required but can be used to identify these marker genes.
This method captures relative abundances of a cell type across individuals. Note that these abundances are not proportions, so they cannot be compared between different cell types. 

## Installation

`devtools::install_github("cozygene/bisque")`

## Getting Started
You can load Bisque as follows:

`library(Bisque)`

The two modes of operation described above are called as follows:

`res <- Bisque::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers)`

`res <- Bisque::MarkerBasedDecomposition(bulk.eset, markers)`

Each method returns a list of results with estimated cell proportions/abundances stored in `res$bulk.props`.

To see examples of these methods on simulated data, check out the vignette:

`browseVignettes("Bisque")`
