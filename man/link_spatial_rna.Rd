\name{link_spatial_rna}
\alias{link_spatial_rna}
\title{Link spatial RNA spots to multiome RNA cells}
\usage{
link_spatial_rna(
  spatial, multiome, spatial_assay = "Spatial", rna_assay = "RNA",
  spatial_coord_cols = c("x", "y"), celltype_col = NULL,
  nfeatures = 10000L, integration_dims = 1:30, harmony_dims = 1:20,
  k_neighbors = 50L, seed = 123L, output_dir = NULL,
  overwrite = FALSE, return_integrated = FALSE, verbose = TRUE
)
}
\arguments{
\item{spatial}{A prepared spatial RNA Seurat object.}
\item{multiome}{A prepared multiome Seurat object containing an RNA assay.}
\item{spatial_assay}{Assay name in \code{spatial}.}
\item{rna_assay}{RNA assay name in \code{multiome}.}
\item{spatial_coord_cols}{Two spatial metadata columns holding x and y.}
\item{celltype_col}{Optional multiome metadata column copied to links.}
\item{nfeatures}{Number of integration features.}
\item{integration_dims}{Dimensions used for integration.}
\item{harmony_dims}{Dimensions used for Harmony and neighbour finding.}
\item{k_neighbors}{Joint nearest-neighbour list size.}
\item{seed}{Random seed.}
\item{output_dir}{Optional standard-output directory.}
\item{overwrite}{Whether existing standard outputs may be replaced.}
\item{return_integrated}{Whether to return and save the integrated object.}
\item{verbose}{Whether Seurat and Harmony print progress.}
}
\value{A list with \code{links}, \code{spots}, \code{qc}, and
\code{unmatched_spots}; optionally also \code{integrated}.}
\description{Uses SCT/CCA integration and Harmony followed by deterministic,
constrained-first, one-to-one greedy assignment. Sample-specific annotation and
QC are intentionally outside this function.}

