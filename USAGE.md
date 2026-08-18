# SpatialATAC Usage Guide

SpatialATAC links spatial transcriptomics data to paired single-cell multiome
RNA+ATAC data and produces an ATAC peak matrix on the complete spatial grid.

The workflow has two public functions:

1. R: `link_spatial_rna()` links Spatial RNA spots one-to-one to multiome RNA
   cells.
2. Python: `place_atac_on_spots()` retrieves the same cells' ATAC profiles,
   places them at spatial coordinates, and applies MAGIC to estimate values for
   spots without a direct match.

## 1. Required input data

Users must prepare three input objects.

### 1.1 Spatial RNA: Seurat object

Requirements:

- Matrix orientation: genes x spatial spots.
- Each column is a spatial spot, and `colnames(spatial_obj)` contains unique
  spot identifiers.
- One assay contains raw RNA counts. Select it with `spatial_assay`; the
  default name is `Spatial`.
- `spatial_obj[[]]` must contain two finite numeric coordinate columns.
  Their names are supplied through `spatial_coord_cols`, which defaults to
  `c("x", "y")`.
- Gene names should have sufficient overlap with the multiome RNA genes.
- Platform-specific QC, low-quality spot filtering, and coordinate preparation
  must be completed before calling this module.

Example validation:

```r
spatial_obj
DefaultAssay(spatial_obj)
head(colnames(spatial_obj))
head(spatial_obj[[]][, c("x", "y")])

stopifnot(
  !anyDuplicated(colnames(spatial_obj)),
  all(c("x", "y") %in% colnames(spatial_obj[[]])),
  all(is.finite(spatial_obj$x)),
  all(is.finite(spatial_obj$y))
)
```

### 1.2 Multiome RNA: Seurat object

This object contains the RNA measurements from the paired multiome cells.

Requirements:

- Matrix orientation: genes x cells.
- Each column is one cell, and `colnames(multiome_obj)` contains unique cell
  identifiers.
- One assay contains raw RNA counts. Select it with `rna_assay`; the default
  name is `RNA`.
- An optional metadata column such as `celltype` can be copied into the link
  table by setting `celltype_col`.
- These cells must represent the same paired cells as the ATAC AnnData object
  described below.

```r
multiome_obj
head(colnames(multiome_obj))
head(multiome_obj[[]])

stopifnot(!anyDuplicated(colnames(multiome_obj)))
```

### 1.3 Multiome ATAC: H5AD / AnnData

This object contains the ATAC measurements from the same multiome cells.

Requirements:

- Matrix orientation: cells x peaks.
- `adata.obs_names` contains unique cell identifiers.
- `adata.var_names` contains unique peak identifiers.
- The recommended peak naming convention is `chr:start-end`, for example
  `chr1:1000-1500`.
- `adata.X` contains the selected ATAC accessibility counts.
- Sparse CSR or CSC storage is strongly recommended because peak matrices can
  be very large.
- Additional `adata.obs` fields such as cell type, QC metrics, and batch are
  allowed. They are copied to the output spot metadata with an `atac_` prefix.

```python
import anndata as ad

atac = ad.read_h5ad("atac.h5ad")
print(atac)                 # observations x variables: cells x peaks
print(atac.obs_names[:3])
print(atac.var_names[:3])

assert atac.obs_names.is_unique
assert atac.var_names.is_unique
```

## 2. RNA and ATAC cell identifiers must match exactly

The workflow uses the known within-cell pairing of multiome RNA and ATAC. It
does not infer a new RNA-to-ATAC correspondence.

The following identifier sets must match:

```text
colnames(multiome_obj)  <==>  atac.obs_names
```

The order may differ, but every identifier must match exactly, including sample
prefixes, barcode suffixes, hyphens, and letter case. The module does not remove
`-1`, add sample names, or perform fuzzy barcode matching.

A recommended check is:

```python
import anndata as ad
import pandas as pd

atac = ad.read_h5ad("atac.h5ad")
rna_ids = pd.read_csv("rna_cell_ids.txt", header=None)[0].astype(str)
assert set(rna_ids) == set(atac.obs_names.astype(str))
```

A filtered ATAC subset can still be used. Linked RNA IDs missing from ATAC are
reported in `missing_atac_ids.txt` and skipped; replacement cells are never
guessed.

## 3. Stage 1: link Spatial RNA to multiome RNA

From the repository root, install the R package if it is not already installed:

```bash
R CMD INSTALL .
```

Run the R stage:

```r
library(spatialatac)

spatial_obj <- readRDS("spatial.rds")
multiome_obj <- readRDS("multiome.rds")

mapping <- link_spatial_rna(
  spatial = spatial_obj,
  multiome = multiome_obj,
  spatial_assay = "Spatial",
  rna_assay = "RNA",
  spatial_coord_cols = c("x", "y"),
  celltype_col = "celltype",
  nfeatures = 10000,
  integration_dims = 1:30,
  harmony_dims = 1:20,
  k_neighbors = 50,
  seed = 123,
  output_dir = "results",
  overwrite = FALSE
)
```

For small datasets, reduce `nfeatures`, the integration and Harmony dimensions,
and `k_neighbors`. Input QC, annotation harmonization, and optional cell-type
subsampling are intentionally left to the user.

Standard outputs:

- `spatial_rna_links.csv`
- `spatial_spots.csv`
- `unmatched_spots.csv`
- `link_qc.csv`

### spatial_rna_links.csv

| Column | Meaning |
|---|---|
| `id_st` | Unique spatial spot identifier |
| `id_rna` | Unique linked multiome RNA cell identifier |
| `dist` | Distance in the Harmony-integrated neighbour space |
| `neighbor_rank` | Rank of the selected RNA neighbour |
| `cluster_rna` | Optional annotation copied from `celltype_col` |

The links are strictly one-to-one for both `id_st` and `id_rna`.

### spatial_spots.csv

Required columns:

| Column | Meaning |
|---|---|
| `id_st` | Unique spatial spot identifier |
| `x` | Numeric x coordinate |
| `y` | Numeric y coordinate |

Other Spatial Seurat metadata fields are also written to this table. Every row
becomes an observation in the final spatial ATAC object, including spots without
a direct RNA-cell link.

## 4. Stage 2: place ATAC profiles on spatial spots

If the Python package is not installed, install it from the repository root:

```bash
python -m pip install -e python
```

Run the Python stage:

```python
from spatial_atac import place_atac_on_spots

spatial_atac = place_atac_on_spots(
    atac="atac.h5ad",
    links="results/spatial_rna_links.csv",
    spots="results/spatial_spots.csv",
    output_dir="results",
    overwrite=False,
    magic_kwargs={"solver": "approximate"},
)

print(spatial_atac)
print(spatial_atac.uns["spatial_atac_qc"])
```

Standard outputs:

- `spatial_atac.h5ad`
- `projection_qc.json`
- `missing_atac_ids.txt`

The final H5AD contains:

- Observations: the complete spatial spot grid.
- Variables: ATAC peaks that are nonzero in at least one successfully mapped
  cell.
- `layers["counts"]`: sparse raw ATAC counts after spatial placement.
- `X`: normalized, log-transformed, MAGIC-imputed values.
- `obsm["spatial"]`: x/y coordinates.
- `obs["has_observed_atac"]`: whether a spot has a directly mapped ATAC cell.
- `obs["source_cell_id"]`: the directly mapped multiome cell identifier.

MAGIC may create a dense matrix. Before running a large dataset, check
`estimated_dense_gb_float32` and ensure that sufficient memory is available.

## 5. Software environment

The project does not create a new environment. It can run directly in an
existing environment containing both R and Python dependencies.

The environment used for validation is named:

```text
R-4.3.0
```

The environment is named `R-4.3.0`, but its currently detected R version is
4.3.2.

### Validated R versions

| Software | Validated version |
|---|---:|
| R | 4.3.2 |
| Seurat | 4.4.0 |
| SeuratObject | 5.0.0 |
| harmony | 1.2.0 |
| withr | 3.0.2 |
| Matrix | 1.6.1 |

The formal R dependencies are declared in `DESCRIPTION`: `Seurat`,
`harmony`, and `withr`. `Matrix` provides sparse matrix support.
`testthat` is needed only for tests.

### Validated Python versions

| Software | Validated version |
|---|---:|
| Python | 3.12.0 |
| anndata | 0.10.8 |
| numpy | 1.26.4 |
| pandas | 2.0.3 |
| scipy | 1.13.0 |
| scanpy | 1.10.1 |
| graphtools | 1.5.3 |
| magic-impute | 3.0.0 |

The formal Python dependencies are declared in `python/pyproject.toml`.
`magic-impute` must be available through `scanpy.external.pp.magic`.

To run explicitly with the validated environment:

```bash
conda activate R-4.3.0
R_BIN=Rscript
PY_BIN=python

"$R_BIN" your_link_script.R
PYTHONPATH=python/src "$PY_BIN" your_projection_script.py
```

## 6. Pseudodata

The root-level `pseudodata/` directory contains fully artificial format
examples and no real sample information.

Generate the inputs:

```bash
Rscript \
  pseudodata/generate_seurat_inputs.R

python \
  pseudodata/generate_atac_input.py
```

The generated data contain:

- Spatial RNA: 120 genes x 150 spots.
- Multiome RNA: 120 genes x 180 cells.
- Multiome ATAC: 180 cells x 150 peaks.
- Identical RNA and ATAC cell IDs: `cell_001` through `cell_180`.

See `pseudodata/README.md` for the complete pseudodata example.
