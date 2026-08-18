import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from spatial_atac import place_atac_on_spots


@pytest.fixture
def inputs():
    atac = ad.AnnData(
        X=sparse.csr_matrix([[1, 0, 2], [0, 3, 0]], dtype=np.float32),
        obs=pd.DataFrame({"celltype": ["A", "B"]}, index=["c1", "c2"]),
        var=pd.DataFrame(index=["p1", "p2", "p3"]),
    )
    links = pd.DataFrame(
        {
            "id_st": ["s1", "s2", "s3"],
            "id_rna": ["c1", "missing", "c2"],
            "dist": [0.1, 0.2, 0.3],
        }
    )
    spots = pd.DataFrame(
        {"id_st": ["s1", "s2", "s3", "s4"], "x": [1, 2, 3, 4], "y": [5, 6, 7, 8]}
    )
    return atac, links, spots


def test_full_grid_counts_and_missing_ids(inputs, monkeypatch):
    atac, links, spots = inputs

    def fake_magic(adata, **kwargs):
        adata.X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X

    monkeypatch.setattr("scanpy.external.pp.magic", fake_magic)
    with pytest.warns(UserWarning, match="absent"):
        result = place_atac_on_spots(atac, links, spots)

    assert result.shape == (4, 3)
    assert result.obs_names.tolist() == ["s1", "s2", "s3", "s4"]
    assert result.obs["has_observed_atac"].tolist() == [True, False, True, False]
    np.testing.assert_array_equal(
        result.layers["counts"].toarray(),
        np.array([[1, 0, 2], [0, 0, 0], [0, 3, 0], [0, 0, 0]], dtype=np.float32),
    )
    np.testing.assert_array_equal(result.obsm["spatial"], spots[["x", "y"]].to_numpy())
    assert result.uns["spatial_atac_qc"]["missing_atac_ids"] == ["missing"]


def test_rejects_non_one_to_one_links(inputs):
    atac, links, spots = inputs
    links.loc[1, "id_st"] = "s1"
    with pytest.raises(ValueError, match="one-to-one"):
        place_atac_on_spots(atac, links, spots)

