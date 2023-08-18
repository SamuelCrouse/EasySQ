# Sam Crouse
# scrouse2@uwyo.edu
# follows the squidpy analyze vizgen tutorial

# python imports
import matplotlib.pyplot as plt
import os

# my imports
import EasySQ as eSQ
import EasySQTools as tools


# sq demo imports
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq


def esq_demo_1():
    path = 'F:/sunlabmerfishdata/QSFL01222023/'
    path = os.getcwd().replace('\\', '/') + '/tutorial_data/'
    esq = eSQ.Analysis(data_path=path)
    esq.print()

    perUn = esq.qcMetrics()
    print("percent unassigned: {}".format(perUn))

    esq.plotTranscripts(show=False)

    esq.filterCells(minCounts=10)

    print("layers")
    esq.layers()
    print("highly variable genes")
    esq.highlyVariableGenes()
    print("normalize total")
    esq.normalizeTotal()
    print("log1p")
    esq.log1p()
    print("pca")
    esq.pp_pca()
    print("neighbors")
    esq.neighbors()
    print("umap")
    esq.tl_umap()
    print("leiden")
    print("\n\ncolors before leiden: {}".format(tools.getLeidenColors(esq.getAdata())))
    esq.leiden(resolution=10)
    print("\n\ncolors after leiden: {}".format(tools.getLeidenColors(esq.getAdata())))

    print("pl_UMAP")

    esq.pl_umap(graphs=["total_counts", "n_genes_by_counts", "leiden"])
    esq.spatialScatter(graphs=["leiden"])
    esq.showPlots()

    print("spatial neighbors")
    esq.spatialNeighbors(delaunay=True)

    print(tools.getLeidenColors(esq.getAdata()))


def sq_demo_1():
    print("import")
    vizgen_dir = Path().resolve() / "tutorial_data" / "vizgen_data"

    adata = sq.read.vizgen(
        path=vizgen_dir,
        counts_file="datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_by_gene_S1R1.csv",
        meta_file="datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_metadata_S1R1.csv",
        transformation_file="datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_images_micron_to_mosaic_pixel_transform.csv",
    )

    sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100, 200, 300), inplace=True)

    adata.obsm["blank_genes"].to_numpy().sum() / adata.var["total_counts"].sum() * 100

    fig, axs = plt.subplots(1, 4, figsize=(15, 4))

    axs[0].set_title("Total transcripts per cell")
    sns.histplot(
        adata.obs["total_counts"],
        kde=False,
        ax=axs[0],
    )

    axs[1].set_title("Unique transcripts per cell")
    sns.histplot(
        adata.obs["n_genes_by_counts"],
        kde=False,
        ax=axs[1],
    )

    axs[2].set_title("Transcripts per FOV")
    sns.histplot(
        adata.obs.groupby("fov").sum()["total_counts"],
        kde=False,
        ax=axs[2],
    )

    axs[3].set_title("Volume of segmented cells")
    sns.histplot(
        adata.obs["volume"],
        kde=False,
        ax=axs[3],
    )

    sc.pp.filter_cells(adata, min_counts=10)

    adata.layers["counts"] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    sc.pl.umap(
        adata,
        color=[
            "total_counts",
            "n_genes_by_counts",
            "leiden",
        ],
        wspace=0.4,
    )

    sq.pl.spatial_scatter(
        adata,
        shape=None,
        color=[
            "leiden",
        ],
        wspace=0.4,
    )

    plt.show()


if __name__ == "__main__":
    esq_demo_1()
    # sq_demo_1()