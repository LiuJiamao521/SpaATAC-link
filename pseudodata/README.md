# Pseudodata example

从仓库根目录、现有 R-4.3.0 环境运行：

```bash
Rscript pseudodata/generate_seurat_inputs.R
python pseudodata/generate_atac_input.py
```

生成：

- `data/spatial.rds`：120 genes × 150 spots，assay 为 `Spatial`
- `data/multiome.rds`：120 genes × 180 cells，assay 为 `RNA`
- `data/atac.h5ad`：180 cells × 150 peaks
- 四个 CSV：便于直接查看 metadata、cell IDs 和 peak IDs

RNA 和 ATAC 的细胞名都严格为 `cell_001` 至 `cell_180`。

R 阶段：

```r
library(spatialatac)
mapping <- link_spatial_rna(
  readRDS("pseudodata/data/spatial.rds"),
  readRDS("pseudodata/data/multiome.rds"),
  spatial_assay = "Spatial", rna_assay = "RNA",
  spatial_coord_cols = c("x", "y"), celltype_col = "celltype",
  nfeatures = 100, integration_dims = 1:10, harmony_dims = 1:8,
  k_neighbors = 10,
  output_dir = "pseudodata/results", overwrite = TRUE
)
```

Python 阶段：

```python
from spatial_atac import place_atac_on_spots
result = place_atac_on_spots(
    "pseudodata/data/atac.h5ad",
    "pseudodata/results/spatial_rna_links.csv",
    "pseudodata/results/spatial_spots.csv",
    output_dir="pseudodata/results", overwrite=True,
)
```
