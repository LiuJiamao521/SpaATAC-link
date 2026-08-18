suppressPackageStartupMessages(library(Seurat))
set.seed(2026)
output_dir <- "pseudodata/data"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
genes <- sprintf("GENE_%03d", 1:120)
spot_ids <- sprintf("spot_%03d", 1:150)
cell_ids <- sprintf("cell_%03d", 1:180)
gene_baseline <- rgamma(length(genes), 2, 1)
spot_group <- rep(1:4, length.out = length(spot_ids))
cell_group <- rep(1:4, length.out = length(cell_ids))
effects <- matrix(rgamma(length(genes) * 4, 1.5, 2), nrow = length(genes))
make_counts <- function(ids, groups) {
  x <- vapply(groups, function(g) {
    rpois(length(genes), gene_baseline + 2 * effects[, g])
  }, numeric(length(genes)))
  dimnames(x) <- list(genes, ids)
  methods::as(x, "dgCMatrix")
}
spatial_obj <- CreateSeuratObject(
  make_counts(spot_ids, spot_group), assay = "Spatial"
)
spatial_obj$x <- rep(seq(0, 90, 10), length.out = 150)
spatial_obj$y <- rep(seq(0, 50, 10), each = 10, length.out = 150)
spatial_obj$region <- paste0("region_", spot_group)
multiome_obj <- CreateSeuratObject(
  make_counts(cell_ids, cell_group), assay = "RNA"
)
multiome_obj$celltype <- paste0("type_", cell_group)
saveRDS(spatial_obj, file.path(output_dir, "spatial.rds"))
saveRDS(multiome_obj, file.path(output_dir, "multiome.rds"))
write.csv(data.frame(id_st = colnames(spatial_obj), spatial_obj[[]]),
          file.path(output_dir, "spatial_metadata.csv"), row.names = FALSE)
write.csv(data.frame(id_rna = colnames(multiome_obj), multiome_obj[[]]),
          file.path(output_dir, "multiome_metadata.csv"), row.names = FALSE)
message("Wrote pseudodata Seurat inputs to ", output_dir)
