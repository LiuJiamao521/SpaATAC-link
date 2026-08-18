# 用户输入数据与格式

本流程需要三项生物学输入。R 阶段使用 Spatial RNA 和 multiome RNA；
Python 阶段使用配对的 multiome ATAC，以及 R 阶段输出的两个 CSV。

## 1. Spatial RNA：Seurat 对象

- 矩阵方向为 genes × spots，列名是唯一的空间 spot ID。
- 一个 assay 保存原始 counts；名称通过 `spatial_assay` 指定，默认 `Spatial`。
- `spatial[[]]` 中必须有两个无缺失的数值坐标列；通过
  `spatial_coord_cols` 指定，默认 `c("x", "y")`。
- 基因名称应与 multiome RNA 有足够交集。
- 平台特异的 QC、过滤和注释整理应在调用本模块之前完成。

## 2. Multiome RNA：Seurat 对象

- 矩阵方向为 genes × cells，列名是唯一的细胞 ID。
- 一个 assay 保存原始 RNA counts；名称通过 `rna_assay` 指定，默认 `RNA`。
- 可选 metadata 列（例如 `celltype`）可通过 `celltype_col` 写入 links。
- 最重要：其细胞名必须与同一批细胞的 ATAC `obs_names` 完全一致。

顺序可以不同，但字符、前后缀和大小写必须一致。本模块不会猜测或修改 barcode：

```r
stopifnot(!anyDuplicated(colnames(multiome_obj)))
# 应满足 setequal(colnames(multiome_obj), atac_cell_ids)
```

## 3. Multiome ATAC：H5AD / AnnData

- 矩阵方向为 cells × peaks。
- `adata.obs_names` 是唯一细胞 ID，与 multiome RNA 列名完全一致。
- `adata.var_names` 是唯一 peak ID，推荐 `chr:start-end`，如
  `chr1:1000-1500`。
- `adata.X` 保存 ATAC accessibility counts；强烈建议稀疏 CSR/CSC 矩阵。
- 额外 `obs` 列会以 `atac_` 前缀复制到输出 spot metadata。

```python
import anndata as ad
atac = ad.read_h5ad("atac.h5ad")
assert atac.obs_names.is_unique
assert atac.var_names.is_unique
print(atac.shape)  # cells x peaks
```

## 两阶段之间的 CSV

R 函数自动生成：

### spatial_rna_links.csv

| 列 | 要求 | 含义 |
|---|---|---|
| `id_st` | 必需、唯一 | spatial spot ID |
| `id_rna` | 必需、唯一 | multiome cell ID，精确匹配 ATAC obs_names |
| `dist` | R 自动生成 | 整合空间中的距离 |
| `neighbor_rank` | R 自动生成 | 被选邻居的排名 |
| `cluster_rna` | 可选 | 从 celltype_col 复制的注释 |

`id_st` 和 `id_rna` 都必须一对一。

### spatial_spots.csv

| 列 | 要求 | 含义 |
|---|---|---|
| `id_st` | 必需、唯一 | spatial spot ID |
| `x` | 必需、数值 | x 坐标 |
| `y` | 必需、数值 | y 坐标 |
| 其他列 | 可选 | 传递到最终输出的 spatial metadata |

该表的每一行都会进入最终空间 ATAC，包括未直接匹配到细胞的 spots。

## 最小检查

```r
stopifnot(
  !anyDuplicated(colnames(spatial_obj)),
  !anyDuplicated(colnames(multiome_obj)),
  all(c("x", "y") %in% colnames(spatial_obj[[]]))
)
```

```python
assert atac.obs_names.is_unique
assert atac.var_names.is_unique
assert set(rna_cell_ids) == set(atac.obs_names)
```

可运行的人工示例见 `pseudodata/`。这些数据只用于展示格式和
smoke test，没有任何生物学含义。
