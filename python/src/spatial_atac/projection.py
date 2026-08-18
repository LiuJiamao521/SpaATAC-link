"""Projection of paired single-cell ATAC profiles onto spatial RNA spots."""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Any, Mapping

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def place_atac_on_spots(
    atac: ad.AnnData | str | Path,
    links: pd.DataFrame | str | Path,
    spots: pd.DataFrame | str | Path,
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    normalize_target_sum: float | None = None,
    magic_kwargs: Mapping[str, Any] | None = None,
    dense_warning_gb: float = 8.0,
) -> ad.AnnData:
    """Place paired multiome ATAC profiles on a complete spatial spot grid.

    RNA and ATAC cell names are matched exactly through ``links.id_rna``. A
    linked RNA cell absent from the filtered ATAC object is reported and
    skipped. Unlinked spots and spots whose linked ATAC cell is absent receive
    sparse zero-count rows before normalization and MAGIC imputation.

    Parameters
    ----------
    atac
        AnnData peak matrix, or a path to an H5AD file. Observation names must
        be the paired RNA/ATAC cell identifiers.
    links
        DataFrame or CSV containing unique ``id_st`` and ``id_rna`` columns.
        Optional ``dist``, ``neighbor_rank``, and ``cluster_rna`` columns are
        propagated into spot metadata.
    spots
        DataFrame or CSV containing unique ``id_st``, ``x``, and ``y`` columns.
        Every row becomes one output observation, in the supplied order.
    output_dir
        If supplied, write ``spatial_atac.h5ad``, ``projection_qc.json``, and
        ``missing_atac_ids.txt``. Nothing is written otherwise.
    overwrite
        Permit replacement of standard output files.
    normalize_target_sum
        Passed to :func:`scanpy.pp.normalize_total`. The default ``None``
        preserves Scanpy's median-total normalization used by the notebook.
    magic_kwargs
        Extra arguments for :func:`scanpy.external.pp.magic`. The solver
        defaults to ``"approximate"``.
    dense_warning_gb
        Warn when one float32 dense matrix is estimated to exceed this size.

    Returns
    -------
    anndata.AnnData
        Full spatial grid. ``layers['counts']`` holds sparse raw counts and
        ``X`` holds normalized, log-transformed, MAGIC-imputed values.
    """
    atac_obj = _read_atac(atac)
    links_df = _read_table(links, "links")
    spots_df = _read_table(spots, "spots")
    _validate_tables(atac_obj, links_df, spots_df)

    # Exact matching is intentional: no prefixes, suffixes, or fuzzy barcode
    # transformations are inferred by this package.
    atac_ids = pd.Index(atac_obj.obs_names.astype(str))
    link_ids = links_df["id_rna"].astype(str)
    present = link_ids.isin(atac_ids)
    matched_links = links_df.loc[present].copy()
    missing_ids = links_df.loc[~present, "id_rna"].astype(str).tolist()
    if missing_ids:
        warnings.warn(
            f"{len(missing_ids)} of {len(links_df)} linked RNA cells are absent "
            "from the ATAC object and will be skipped.",
            UserWarning,
            stacklevel=2,
        )

    matched_links["id_st"] = matched_links["id_st"].astype(str)
    matched_links["id_rna"] = matched_links["id_rna"].astype(str)
    spots_df = spots_df.copy()
    spots_df["id_st"] = spots_df["id_st"].astype(str)
    spot_index = pd.Index(spots_df["id_st"], name="id_st")
    matched_links = matched_links[matched_links["id_st"].isin(spot_index)].copy()

    linked_atac = atac_obj[matched_links["id_rna"].tolist(), :].copy()
    linked_x = linked_atac.X
    if not sparse.issparse(linked_x):
        linked_x = sparse.csr_matrix(np.asarray(linked_x))
    else:
        linked_x = linked_x.tocsr()

    output_rows = spot_index.get_indexer(matched_links["id_st"])
    placement = sparse.csr_matrix(
        (
            np.ones(len(matched_links), dtype=np.float32),
            (output_rows, np.arange(len(matched_links))),
        ),
        shape=(len(spots_df), len(matched_links)),
    )
    full_counts = (placement @ linked_x).tocsr()

    # "All non-zero peaks" means peaks observed in at least one successfully
    # mapped ATAC cell. Zero-filled spatial rows do not affect this criterion.
    nonzero_features = np.asarray(full_counts.getnnz(axis=0)).ravel() > 0
    full_counts = full_counts[:, nonzero_features].astype(np.float32)
    var = linked_atac.var.iloc[np.flatnonzero(nonzero_features)].copy()

    obs = spots_df.set_index("id_st", drop=True).copy()
    obs.index = obs.index.astype(str)
    source_cell = pd.Series(
        matched_links["id_rna"].to_numpy(),
        index=matched_links["id_st"],
        dtype="object",
    )
    obs["source_cell_id"] = source_cell.reindex(obs.index)
    obs["has_observed_atac"] = obs["source_cell_id"].notna()
    for column in ("dist", "neighbor_rank", "cluster_rna"):
        if column in matched_links.columns:
            values = pd.Series(
                matched_links[column].to_numpy(), index=matched_links["id_st"]
            )
            obs[f"link_{column}"] = values.reindex(obs.index)

    atac_metadata = linked_atac.obs.copy()
    atac_metadata.index = matched_links["id_st"].to_numpy()
    atac_metadata = atac_metadata.add_prefix("atac_").reindex(obs.index)
    collisions = obs.columns.intersection(atac_metadata.columns)
    if len(collisions):
        raise ValueError(
            "Spatial metadata collides with prefixed ATAC metadata: "
            + ", ".join(collisions)
        )
    obs = pd.concat([obs, atac_metadata], axis=1)

    result = ad.AnnData(X=full_counts.copy(), obs=obs, var=var)
    result.layers["counts"] = full_counts.copy()
    result.obsm["spatial"] = spots_df[["x", "y"]].to_numpy(dtype=float)

    sc.pp.normalize_total(result, target_sum=normalize_target_sum)
    sc.pp.log1p(result)

    estimated_dense_gb = (
        result.n_obs * result.n_vars * np.dtype(np.float32).itemsize / 1024**3
    )
    if estimated_dense_gb > dense_warning_gb:
        warnings.warn(
            "MAGIC may materialize dense matrices; one float32 matrix is "
            f"approximately {estimated_dense_gb:.1f} GiB for this dataset.",
            ResourceWarning,
            stacklevel=2,
        )
    kwargs: dict[str, Any] = {"solver": "approximate"}
    if magic_kwargs is not None:
        kwargs.update(dict(magic_kwargs))
    try:
        sc.external.pp.magic(result, **kwargs)
    except ImportError as exc:
        raise RuntimeError(
            "MAGIC imputation requires the existing environment to provide "
            "the `magic-impute` package. No environment is created automatically."
        ) from exc

    qc = {
        "n_spots": int(len(spots_df)),
        "n_input_links": int(len(links_df)),
        "n_atac_cells": int(atac_obj.n_obs),
        "n_observed_spots": int(len(matched_links)),
        "n_missing_atac_ids": int(len(missing_ids)),
        "observed_spot_fraction": float(len(matched_links) / len(spots_df)),
        "n_nonzero_features": int(result.n_vars),
        "estimated_dense_gb_float32": float(estimated_dense_gb),
        "normalization_target_sum": normalize_target_sum,
        "magic_kwargs": kwargs,
        "missing_atac_ids": missing_ids,
    }
    result.uns["spatial_atac_qc"] = qc

    if output_dir is not None:
        _write_outputs(result, qc, missing_ids, Path(output_dir), overwrite)
    return result


def _read_atac(value: ad.AnnData | str | Path) -> ad.AnnData:
    if isinstance(value, ad.AnnData):
        return value
    path = Path(value)
    if not path.is_file():
        raise FileNotFoundError(f"ATAC H5AD not found: {path}")
    return sc.read_h5ad(path)


def _read_table(value: pd.DataFrame | str | Path, label: str) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        return value.copy()
    path = Path(value)
    if not path.is_file():
        raise FileNotFoundError(f"{label} table not found: {path}")
    return pd.read_csv(path)


def _validate_tables(
    atac: ad.AnnData, links: pd.DataFrame, spots: pd.DataFrame
) -> None:
    link_required = {"id_st", "id_rna"}
    spot_required = {"id_st", "x", "y"}
    if missing := link_required.difference(links.columns):
        raise ValueError(f"links is missing required columns: {sorted(missing)}")
    if missing := spot_required.difference(spots.columns):
        raise ValueError(f"spots is missing required columns: {sorted(missing)}")
    if links[["id_st", "id_rna"]].isna().any().any():
        raise ValueError("links contains missing spot or RNA cell IDs")
    if spots[["id_st", "x", "y"]].isna().any().any():
        raise ValueError("spots contains missing IDs or coordinates")
    if links["id_st"].duplicated().any() or links["id_rna"].duplicated().any():
        raise ValueError("links must be one-to-one in both id_st and id_rna")
    if spots["id_st"].duplicated().any():
        raise ValueError("spots.id_st must be unique")
    if not atac.obs_names.is_unique or not atac.var_names.is_unique:
        raise ValueError("ATAC observation and feature names must be unique")
    unknown_spots = ~links["id_st"].astype(str).isin(spots["id_st"].astype(str))
    if unknown_spots.any():
        bad = links.loc[unknown_spots, "id_st"].astype(str).tolist()[:5]
        raise ValueError(f"links contains spot IDs absent from spots: {bad}")


def _write_outputs(
    result: ad.AnnData,
    qc: Mapping[str, Any],
    missing_ids: list[str],
    output_dir: Path,
    overwrite: bool,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    h5ad_path = output_dir / "spatial_atac.h5ad"
    qc_path = output_dir / "projection_qc.json"
    missing_path = output_dir / "missing_atac_ids.txt"
    paths = (h5ad_path, qc_path, missing_path)
    existing = [str(path) for path in paths if path.exists()]
    if existing and not overwrite:
        raise FileExistsError("Refusing to overwrite: " + ", ".join(existing))
    result.write_h5ad(h5ad_path, compression="gzip")
    qc_path.write_text(json.dumps(qc, indent=2, ensure_ascii=False) + "\n")
    missing_path.write_text("".join(f"{cell_id}\n" for cell_id in missing_ids))

