#' Link spatial RNA spots to multiome RNA cells
#'
#' Integrates a prepared spatial RNA Seurat object with the RNA assay of a
#' prepared single-cell multiome Seurat object, constructs a joint Harmony
#' nearest-neighbour space, and assigns at most one RNA cell to each spatial
#' spot. Spots with fewer RNA candidates are assigned first, which preserves
#' the constrained-first greedy strategy of the original workflow.
#'
#' This function deliberately performs no sample-specific QC, annotation
#' remapping, or cell-type downsampling. Cell names are preserved as the public
#' identifiers even though collision-safe prefixes are used internally.
#'
#' @param spatial A prepared spatial RNA Seurat object.
#' @param multiome A prepared multiome Seurat object containing an RNA assay.
#' @param spatial_assay Assay name in `spatial`.
#' @param rna_assay RNA assay name in `multiome`.
#' @param spatial_coord_cols Two columns in `spatial[[]]` holding x and y.
#' @param celltype_col Optional metadata column in `multiome` to copy to
#'   `cluster_rna` in the link table.
#' @param nfeatures Number of integration features.
#' @param integration_dims Dimensions used for anchor finding/integration.
#' @param harmony_dims Dimensions used for Harmony and neighbour finding.
#' @param k_neighbors Size of the joint nearest-neighbour list.
#' @param seed Random seed for Seurat and Harmony operations.
#' @param output_dir Optional directory for standard CSV/RDS outputs.
#' @param overwrite Whether existing outputs may be replaced.
#' @param return_integrated Include the large integrated Seurat object in the
#'   returned list and, when `output_dir` is set, save it as an RDS.
#' @param verbose Forward progress output to Seurat/Harmony where supported.
#'
#' @return A list containing `links`, `spots`, `qc`, `unmatched_spots`, and,
#'   optionally, `integrated`.
#' @export
link_spatial_rna <- function(
    spatial,
    multiome,
    spatial_assay = "Spatial",
    rna_assay = "RNA",
    spatial_coord_cols = c("x", "y"),
    celltype_col = NULL,
    nfeatures = 10000L,
    integration_dims = 1:30,
    harmony_dims = 1:20,
    k_neighbors = 50L,
    seed = 123L,
    output_dir = NULL,
    overwrite = FALSE,
    return_integrated = FALSE,
    verbose = TRUE) {
  .validate_link_inputs(
    spatial, multiome, spatial_assay, rna_assay, spatial_coord_cols,
    celltype_col, nfeatures, integration_dims, harmony_dims, k_neighbors
  )

  spatial_ids <- colnames(spatial)
  rna_ids <- colnames(multiome)

  # Seurat 4 expects a common RNA assay after merge/SplitObject. Build that
  # internal assay from the user-selected inputs and discard unrelated assays;
  # metadata and original identifiers remain unchanged.
  spatial_counts <- Seurat::GetAssayData(
    spatial, assay = spatial_assay, slot = "counts"
  )
  multiome_counts <- Seurat::GetAssayData(
    multiome, assay = rna_assay, slot = "counts"
  )
  spatial@assays <- list(
    RNA = Seurat::CreateAssayObject(counts = spatial_counts)
  )
  multiome@assays <- list(
    RNA = Seurat::CreateAssayObject(counts = multiome_counts)
  )
  Seurat::DefaultAssay(spatial) <- "RNA"
  Seurat::DefaultAssay(multiome) <- "RNA"

  spatial[[".spatialatac_source"]] <- "spatial"
  spatial[[".spatialatac_original_id"]] <- spatial_ids
  multiome[[".spatialatac_source"]] <- "rna"
  multiome[[".spatialatac_original_id"]] <- rna_ids

  spatial <- Seurat::RenameCells(spatial, add.cell.id = "st")
  multiome <- Seurat::RenameCells(multiome, add.cell.id = "rna")

  combined <- merge(spatial, multiome)
  object_list <- Seurat::SplitObject(combined, split.by = ".spatialatac_source")

  integrated <- withr::with_seed(seed, {
    object_list <- lapply(
      object_list,
      function(x) Seurat::SCTransform(x, verbose = verbose)
    )
    features <- Seurat::SelectIntegrationFeatures(
      object.list = object_list,
      nfeatures = nfeatures
    )
    object_list <- Seurat::PrepSCTIntegration(
      object.list = object_list,
      anchor.features = features,
      verbose = verbose
    )
    object_list <- lapply(
      object_list,
      function(x) Seurat::RunPCA(x, features = features, verbose = verbose)
    )
    anchors <- Seurat::FindIntegrationAnchors(
      object.list = object_list,
      normalization.method = "SCT",
      anchor.features = features,
      dims = integration_dims,
      reduction = "cca",
      k.anchor = 20,
      verbose = verbose
    )
    x <- Seurat::IntegrateData(
      anchorset = anchors,
      normalization.method = "SCT",
      dims = integration_dims,
      verbose = verbose
    )
    x <- Seurat::RunPCA(x, verbose = verbose)
    x <- harmony::RunHarmony(
      x,
      group.by.vars = ".spatialatac_source",
      reduction.use = "pca",
      dims.use = harmony_dims,
      project.dim = FALSE,
      reduction.save = "harmony",
      verbose = verbose
    )
    Seurat::FindNeighbors(
      x,
      k.param = k_neighbors,
      reduction = "harmony",
      dims = harmony_dims,
      return.neighbor = TRUE,
      assay = "integrated",
      verbose = verbose
    )
  })

  neighbour <- .get_joint_neighbour(integrated)
  assignment <- .greedy_one_to_one(integrated, neighbour)
  links <- assignment$links

  rna_meta <- multiome[[]]
  rna_lookup <- stats::setNames(
    rna_meta$.spatialatac_original_id,
    rownames(rna_meta)
  )
  st_meta_internal <- spatial[[]]
  st_lookup <- stats::setNames(
    st_meta_internal$.spatialatac_original_id,
    rownames(st_meta_internal)
  )

  if (nrow(links) > 0L) {
    links$id_st <- unname(st_lookup[links$internal_st])
    links$id_rna <- unname(rna_lookup[links$internal_rna])
    links <- links[, c("id_st", "id_rna", "dist", "neighbor_rank"), drop = FALSE]
    if (!is.null(celltype_col)) {
      type_lookup <- stats::setNames(rna_meta[[celltype_col]], rna_meta$.spatialatac_original_id)
      links$cluster_rna <- unname(type_lookup[links$id_rna])
    }
  } else {
    links <- data.frame(
      id_st = character(), id_rna = character(), dist = numeric(),
      neighbor_rank = integer(), stringsAsFactors = FALSE
    )
    if (!is.null(celltype_col)) links$cluster_rna <- character()
  }

  spatial_meta <- spatial[[]]
  spots <- data.frame(
    id_st = spatial_meta$.spatialatac_original_id,
    x = spatial_meta[[spatial_coord_cols[[1L]]]],
    y = spatial_meta[[spatial_coord_cols[[2L]]]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  public_meta <- spatial_meta[, setdiff(
    colnames(spatial_meta),
    c(".spatialatac_source", ".spatialatac_original_id", spatial_coord_cols)
  ), drop = FALSE]
  rownames(public_meta) <- NULL
  spots <- cbind(spots, public_meta)

  matched_ids <- links$id_st
  unmatched_spots <- spots[!spots$id_st %in% matched_ids, c("id_st", "x", "y"), drop = FALSE]
  qc <- .link_qc(spatial_ids, rna_ids, links, unmatched_spots, k_neighbors, seed)
  result <- list(
    links = links,
    spots = spots,
    qc = qc,
    unmatched_spots = unmatched_spots
  )
  if (isTRUE(return_integrated)) result$integrated <- integrated

  if (!is.null(output_dir)) {
    .write_link_outputs(result, output_dir, overwrite, return_integrated)
  }
  result
}

.validate_link_inputs <- function(
    spatial, multiome, spatial_assay, rna_assay, spatial_coord_cols,
    celltype_col, nfeatures, integration_dims, harmony_dims, k_neighbors) {
  if (!inherits(spatial, "Seurat") || !inherits(multiome, "Seurat")) {
    stop("`spatial` and `multiome` must both be Seurat objects.", call. = FALSE)
  }
  if (!spatial_assay %in% names(spatial@assays)) {
    stop("Spatial assay not found: ", spatial_assay, call. = FALSE)
  }
  if (!rna_assay %in% names(multiome@assays)) {
    stop("Multiome RNA assay not found: ", rna_assay, call. = FALSE)
  }
  if (length(spatial_coord_cols) != 2L ||
      !all(spatial_coord_cols %in% colnames(spatial[[]]))) {
    stop("`spatial_coord_cols` must name two spatial metadata columns.", call. = FALSE)
  }
  if (!is.null(celltype_col) && !celltype_col %in% colnames(multiome[[]])) {
    stop("Cell-type column not found: ", celltype_col, call. = FALSE)
  }
  if (anyDuplicated(colnames(spatial)) || anyDuplicated(colnames(multiome))) {
    stop("Cell names must be unique within each input object.", call. = FALSE)
  }
  if (ncol(spatial) < 2L || ncol(multiome) < 2L) {
    stop("Each input must contain at least two cells.", call. = FALSE)
  }
  if (nfeatures < 1L || k_neighbors < 2L ||
      !length(integration_dims) || !length(harmony_dims)) {
    stop("Integration features/dimensions and neighbour parameters are invalid.", call. = FALSE)
  }
  invisible(TRUE)
}

.get_joint_neighbour <- function(integrated) {
  neighbours <- integrated@neighbors
  if (!length(neighbours)) {
    stop("Seurat did not return a neighbour object.", call. = FALSE)
  }
  valid <- vapply(neighbours, function(x) {
    inherits(x, "Neighbor") && length(x@cell.names) == ncol(integrated)
  }, logical(1))
  if (!any(valid)) {
    stop("No joint neighbour object matches the integrated cells.", call. = FALSE)
  }
  neighbours[[tail(which(valid), 1L)]]
}

.greedy_one_to_one <- function(integrated, neighbour) {
  cell_names <- neighbour@cell.names
  nn_idx <- neighbour@nn.idx
  nn_dist <- neighbour@nn.dist
  source <- integrated[[]][cell_names, ".spatialatac_source", drop = TRUE]
  spatial_rows <- which(source == "spatial")
  rna_cells <- cell_names[source == "rna"]

  candidates <- lapply(spatial_rows, function(i) {
    idx <- nn_idx[i, ]
    keep <- idx > 0L & idx <= length(cell_names)
    ids <- cell_names[idx[keep]]
    dists <- nn_dist[i, keep]
    is_rna <- ids %in% rna_cells
    frame <- data.frame(
      internal_rna = ids[is_rna],
      dist = as.numeric(dists[is_rna]),
      stringsAsFactors = FALSE
    )
    if (nrow(frame)) {
      frame <- frame[order(frame$dist, frame$internal_rna), , drop = FALSE]
      frame$neighbor_rank <- seq_len(nrow(frame))
    } else {
      frame$neighbor_rank <- integer()
    }
    frame
  })
  names(candidates) <- cell_names[spatial_rows]
  candidate_count <- lengths(candidates)
  process_order <- order(candidate_count, names(candidates))
  used_rna <- character()
  rows <- vector("list", length(candidates))
  out_i <- 0L

  for (candidate_i in process_order) {
    available <- candidates[[candidate_i]]
    available <- available[!available$internal_rna %in% used_rna, , drop = FALSE]
    if (!nrow(available)) next
    chosen <- available[1L, , drop = FALSE]
    out_i <- out_i + 1L
    rows[[out_i]] <- data.frame(
      internal_st = names(candidates)[candidate_i],
      internal_rna = chosen$internal_rna,
      dist = chosen$dist,
      neighbor_rank = chosen$neighbor_rank,
      stringsAsFactors = FALSE
    )
    used_rna <- c(used_rna, chosen$internal_rna)
  }
  rows <- rows[seq_len(out_i)]
  links <- if (length(rows)) do.call(rbind, rows) else data.frame()
  list(links = links, candidate_count = candidate_count)
}

.link_qc <- function(spatial_ids, rna_ids, links, unmatched, k, seed) {
  distances <- links$dist
  data.frame(
    n_spatial_spots = length(spatial_ids),
    n_rna_cells = length(rna_ids),
    n_links = nrow(links),
    n_unmatched_spots = nrow(unmatched),
    link_fraction = if (length(spatial_ids)) nrow(links) / length(spatial_ids) else NA_real_,
    median_distance = if (length(distances)) stats::median(distances) else NA_real_,
    max_distance = if (length(distances)) max(distances) else NA_real_,
    k_neighbors = k,
    seed = seed,
    stringsAsFactors = FALSE
  )
}

.write_link_outputs <- function(result, output_dir, overwrite, return_integrated) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    links = file.path(output_dir, "spatial_rna_links.csv"),
    spots = file.path(output_dir, "spatial_spots.csv"),
    qc = file.path(output_dir, "link_qc.csv"),
    unmatched = file.path(output_dir, "unmatched_spots.csv")
  )
  if (isTRUE(return_integrated)) {
    paths <- c(paths, integrated = file.path(output_dir, "integrated_seurat.rds"))
  }
  existing <- paths[file.exists(paths)]
  if (length(existing) && !isTRUE(overwrite)) {
    stop("Refusing to overwrite existing output(s): ",
         paste(existing, collapse = ", "), call. = FALSE)
  }
  utils::write.csv(result$links, paths[["links"]], row.names = FALSE)
  utils::write.csv(result$spots, paths[["spots"]], row.names = FALSE)
  utils::write.csv(result$qc, paths[["qc"]], row.names = FALSE)
  utils::write.csv(result$unmatched_spots, paths[["unmatched"]], row.names = FALSE)
  if (isTRUE(return_integrated)) saveRDS(result$integrated, paths[["integrated"]])
  invisible(paths)
}

