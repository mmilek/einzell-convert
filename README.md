# Convert your scRNAseq data for Einzell

Converts single-cell RNAseq datasets into the JSON format needed by Einzell - single cell synthesizer.

## Required input: a CSV with `x-coordinate`, `y-coordinate`, `cell_type`

```csv
x,y,cell_type
1.38,-6.83,Astrocyte
1.70,11.61,Microglia
0.60,10.91,Astrocyte
...
```

That's the minimum: a 2D embedding coordinate (UMAP, t-SNE, or any other 2D
layout you already have) and a cell-type/cluster label, one row per cell.
Everything else Einzell needs is derived automatically.

This tool deliberately does **not** read `.h5ad` or `.rds` files directly.
If that's what you have, export `x, y, cell_type` to a CSV yourself first —
it's a one-liner either way:

```python
# Python / anndata
adata.obs[["cell_type"]].assign(x=adata.obsm["X_umap"][:, 0], y=adata.obsm["X_umap"][:, 1]) \
    .to_csv("cells.csv", index=False)
```

```r
# R / Seurat
umap <- Embeddings(seurat_obj, "umap")
write.csv(data.frame(x = umap[,1], y = umap[,2], cell_type = seurat_obj$cell_type), "cells.csv", row.names = FALSE)
```

### Optional columns

If you already have them, `pseudotime`, `local_density`, `global_density`,
and `dist_to_peak` columns are used as-is. Otherwise they're computed from
`x`/`y`/`cell_type` directly (see "How conversion works" below).

### Cell-count handling

Einzell always works with exactly **1024 cells** per dataset. Cells are
selected with a *proportional* sample, not a plain random sample: cell-type
proportions are computed first, a multinomial draw allocates how many of
the 1024 slots each cell type gets, and that many cells are randomly sampled
from within each type. See `sample_cells()` in `einzell_convert/core.py` for 
details.

Your CSV needs **at least 8 distinct cell types** and **at least 1024 rows total**, or
conversion will fail. Only the 12 most frequent cell types are kept if you
have more than that.

## Usage

```bash
pip install -r requirements.txt
python convert.py cells.csv output.json
```

## How conversion works

1. **Sample** proportionally down to 1024 cells (see above).
2. **Local density**: inverse mean distance to each cell's 15 nearest
   neighbors in (x, y) space.
3. **Global density**: a Gaussian KDE evaluated over the (x, y) embedding.
4. **Distance to peak**: for each cell type, the cell with the highest local
   density is its "peak"; every other cell of that type gets its distance to
   that peak.
5. **Pseudotime**: This tool runs a diffusion-pseudotime algorithm
   (`scanpy`'s `sc.tl.diffmap` + `sc.tl.dpt`) directly on the 2D (x, y)
   embedding you provide. Since this is a lower-dimensional input
   it might not reflect pseudotime values computed from raw expression
   data.

Einzell-compatible output is a
dict keyed by a fresh 0-based sequential string index, each value
`[cell_type, x, y, pseudotime, local_density, global_density, dist_to_peak]`.

