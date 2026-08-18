"""Generate an artificial cells-by-peaks ATAC AnnData input."""
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

rng = np.random.default_rng(2026)
out = Path("pseudodata/data")
out.mkdir(parents=True, exist_ok=True)
cells = pd.Index([f"cell_{i:03d}" for i in range(1, 181)], name="id_rna")
starts = np.arange(1000, 1000 + 150 * 500, 500)
peaks = pd.Index([f"chr1:{x}-{x + 499}" for x in starts], name="peak_id")
counts = rng.negative_binomial(1, 0.92, size=(len(cells), len(peaks)))
counts[rng.random(counts.shape) < 0.85] = 0
obs = pd.DataFrame({
    "celltype": [f"type_{i % 4 + 1}" for i in range(len(cells))],
    "batch": "pseudo_batch",
}, index=cells)
adata = ad.AnnData(
    X=sparse.csr_matrix(counts, dtype=np.float32),
    obs=obs, var=pd.DataFrame(index=peaks),
)
adata.write_h5ad(out / "atac.h5ad", compression="gzip")
obs.reset_index().to_csv(out / "atac_obs.csv", index=False)
adata.var.reset_index().to_csv(out / "atac_var.csv", index=False)
print(adata)
