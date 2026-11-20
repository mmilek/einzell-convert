from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from sklearn.neighbors import NearestNeighbors

MAX_CELLS = 1024
MIN_CELL_TYPES = 8
MAX_CELL_TYPES = 12

REQUIRED_COLUMNS = ("x", "y", "cell_type")
OPTIONAL_COLUMNS = ("pseudotime", "local_density", "global_density", "dist_to_peak")


class ConversionError(ValueError):
    """Raised when the input doesn't qualify or is malformed."""


def validate_input(df: pd.DataFrame) -> None:
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ConversionError(f"missing required column(s): {', '.join(missing)}")
    if df[list(REQUIRED_COLUMNS)].isna().any().any():
        raise ConversionError("required columns must not contain empty/NaN values")

    n_types = df["cell_type"].nunique()
    if n_types < MIN_CELL_TYPES:
        raise ConversionError(
            f"input has only {n_types} distinct cell types; Einzell requires at least "
            f"{MIN_CELL_TYPES} (matching census_api.py's MIN_CELL_TYPES qualification)."
        )
    if len(df) < MAX_CELLS:
        raise ConversionError(
            f"input has only {len(df)} rows; Einzell requires at least {MAX_CELLS} cells "
            f"to draw its proportional {MAX_CELLS}-cell sample from (matching census_api.py's MAX_CELLS)."
        )


def sample_cells(df: pd.DataFrame, max_cells: int = MAX_CELLS, max_cell_types: int = MAX_CELL_TYPES,
                  random_state: int = 42) -> pd.DataFrame:
    """
    Proportional multinomial sample down to `max_cells` rows: keep only the
    `max_cell_types` most frequent cell types, compute their proportions,
    allocate `max_cells` slots across them via a multinomial draw, then
    randomly sample that many rows from within each type.
    """
    top_types = df["cell_type"].value_counts().head(max_cell_types).index
    df = df[df["cell_type"].isin(top_types)]

    proportions = df["cell_type"].value_counts(normalize=True)
    target_counts = pd.Series(
        np.random.RandomState(random_state).multinomial(max_cells, proportions.values),
        index=proportions.index,
    )

    sampled = []
    for cell_type, n in target_counts.items():
        subset = df[df["cell_type"] == cell_type]
        if n > 0 and len(subset) > 0:
            sampled.append(subset.sample(n=min(n, len(subset)), random_state=random_state))
    if not sampled:
        raise ConversionError("no cells sampled -- check cell_type proportions in the input")
    return pd.concat(sampled).reset_index(drop=True)


def compute_local_density(xy: np.ndarray, n_neighbors: int = 15) -> np.ndarray:
    k = min(n_neighbors, len(xy) - 1)
    nbrs = NearestNeighbors(n_neighbors=k).fit(xy)
    distances, _ = nbrs.kneighbors()
    return 1.0 / np.maximum(distances.mean(axis=1), 1e-9)


def compute_global_density(xy: np.ndarray, bw_method: float = 0.1) -> np.ndarray:
    kde = gaussian_kde(xy.T, bw_method=bw_method)
    return kde(xy.T)


def compute_dist_to_peak(xy: np.ndarray, cell_types: np.ndarray, local_density: np.ndarray) -> np.ndarray:
    """For each cell type, find the
    cell with the highest local density ("the peak") and measure every
    same-type cell's distance to it."""
    peak_dist = np.zeros(len(xy))
    for cell_type in np.unique(cell_types):
        mask = cell_types == cell_type
        if mask.any():
            peak_idx = np.flatnonzero(mask)[np.argmax(local_density[mask])]
            peak_coord = xy[peak_idx]
            peak_dist[mask] = np.linalg.norm(xy[mask] - peak_coord, axis=1)
    return peak_dist


def compute_pseudotime(xy: np.ndarray) -> np.ndarray:
    import anndata
    import scanpy as sc

    adata = anndata.AnnData(X=xy.astype(np.float64))
    sc.pp.neighbors(adata, use_rep="X", n_neighbors=15)
    sc.tl.diffmap(adata)
    adata.uns["iroot"] = 0
    sc.tl.dpt(adata)
    pseudotime = adata.obs["dpt_pseudotime"].to_numpy()

    # Cells in a different connected component than the root (common when
    # (x, y) clusters are tightly separated, e.g. synthetic/toy data, or a
    # cell type with no neighbors near the root) get NaN from dpt. 
    # Treat them as "maximally distant from the root" rather than letting
    # NaN slip into the output JSON.
    finite = pseudotime[np.isfinite(pseudotime)]
    fallback = finite.max() if finite.size else 1.0
    pseudotime = np.where(np.isfinite(pseudotime), pseudotime, fallback)
    return pseudotime


def compute_dataset_output(df: pd.DataFrame, max_cells: int = MAX_CELLS) -> dict:
    """
    Full pipeline: validate -> proportional sample -> derive any missing
    optional columns -> output Einzell's dataset JSON schema: a dict keyed by a
    fresh 0-based sequential string index, each value
    [cell_type, x, y, pseudotime, local_density, global_density, dist_to_peak].
    """
    validate_input(df)
    sampled = sample_cells(df, max_cells=max_cells)

    xy = sampled[["x", "y"]].to_numpy(dtype=np.float64)
    cell_types = sampled["cell_type"].to_numpy()

    local_density = (
        sampled["local_density"].to_numpy(dtype=np.float64)
        if "local_density" in sampled.columns
        else compute_local_density(xy)
    )
    global_density = (
        sampled["global_density"].to_numpy(dtype=np.float64)
        if "global_density" in sampled.columns
        else compute_global_density(xy)
    )
    dist_to_peak = (
        sampled["dist_to_peak"].to_numpy(dtype=np.float64)
        if "dist_to_peak" in sampled.columns
        else compute_dist_to_peak(xy, cell_types, local_density)
    )
    pseudotime = (
        sampled["pseudotime"].to_numpy(dtype=np.float64)
        if "pseudotime" in sampled.columns
        else compute_pseudotime(xy)
    )

    output_map = {}
    for i in range(len(sampled)):
        output_map[str(i)] = [
            str(cell_types[i]),
            float(xy[i, 0]),
            float(xy[i, 1]),
            float(pseudotime[i]),
            float(local_density[i]),
            float(global_density[i]),
            float(dist_to_peak[i]),
        ]
    return output_map
