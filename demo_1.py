# Sam Crouse
# scrouse2@uwyo.edu
# follows the squidpy analyze vizgen tutorial

# python imports
import matplotlib.pyplot as plt
import os
import time

# my imports
import EasySQ as esq
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
    t0 = time.time()

    # path = 'F:/sunlabmerfishdata/QSFL01222023/'
    path = os.getcwd() + '\\tutorial_data_1\\'
    esqAn = esq.Analysis(data_path=path)
    esqAn.print()

    esqAn.setLeidenColors(color_file="leiden_generated_random_1.txt")
    # esqAn.setLeidenColors(color_file="leiden_color_set_1_gradient.csv")

    perUn = esqAn.qcMetrics(percentTop=(50, 100, 125))  # (1, 50, 100, 150, 483)
    print("percent unassigned: {}%".format(perUn))

    esqAn.plotTranscripts(show=False)

    esqAn.filterCells(minCounts=1)

    print("layers")
    esqAn.layers()
    print("highly variable genes")
    esqAn.highlyVariableGenes()
    print("normalize total")
    esqAn.normalizeTotal()
    print("log1p")
    esqAn.log1p()
    print("pca")
    esqAn.pp_pca()
    print("neighbors")
    esqAn.pp_neighbors()
    print("umap")
    esqAn.tl_umap()
    print("leiden")
    esqAn.leiden()

    print("pl_UMAP")
    esqAn.pl_umap(graphs=["total_counts", "n_genes_by_counts", "leiden"])
    esqAn.spatialScatter(graphs=["leiden"])

    print("spatial neighbors")
    esqAn.gr_spatialNeighbors(delaunay=True)

    print("compute centrality scores")
    esqAn.gr_centralityScores()
    esqAn.pl_centralityScores(figsize=(16, 5))

    print("co-occurrence probability")
    adata_subsample = esqAn.adataSubsample(fraction=0.5, copy=True)
    esqAn.gr_co_occurrence(adata=adata_subsample)
    esqAn.pl_co_occurrence(adata=adata_subsample)

    esqAn.spatialScatter(adata=adata_subsample, graphs=["leiden"], size=1)

    print("neighbors enrichment analysis")
    esqAn.gr_nhoodEnrichment()
    esqAn.pl_nhoodEnrichment()

    esqAn.spatialScatter(adata=adata_subsample, graphs="leiden", size=0.1)

    print("Ripley's statistics")
    esqAn.gr_ripley()
    esqAn.pl_ripley()
    esqAn.spatialScatter(adata=adata_subsample, graphs="leiden", groups=["0", "1", "3"])

    print("Moran's I score")
    esqAn.gr_spatialNeighbors(adata=adata_subsample, delaunay=True)
    esqAn.gr_spatialAutocorr(adata=adata_subsample, nPerms=100, nJobs=1)
    print(adata_subsample.uns["moranI"].head(10))

    esqAn.spatialScatter(adata=adata_subsample, graphs=["Slc17a7", "Chrm1"])

    t1 = time.time()
    totalTime = t1 - t0
    print("time elapsed: {}".format(totalTime))

    esqAn.showPlots()


# """
def sq_demo_1():
    t0 = time.time()

    print("\nImporting data...")
    vizgen_path = os.getcwd() + '\\'
    meta_data_path = 'F:\\sunlabmerfishdata\\QSFL01222023\\cell_metadata.csv'
    cell_by_gene_path = 'F:\\sunlabmerfishdata\\QSFL01222023\\cell_by_gene.csv'
    transformation_file_path = 'micron_to_mosaic_pixel_transform.csv'

    adata = sq.read.vizgen(path=vizgen_path,
                           counts_file=cell_by_gene_path,
                           meta_file=meta_data_path,
                           transformation_file=transformation_file_path)

    print(adata)

    """

    print("\nCalculating QC metrics...")
    sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100, 125), inplace=True)
    per_unassigned = adata.obsm["blank_genes"].to_numpy().sum() / adata.var["total_counts"].sum() * 100
    print("Percent Unassigned {}%".format(per_unassigned))

    print("\nPlotting transcript and volume data...")
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

    print("\nRunning a bunch of data calculations...")
    sc.pp.filter_cells(adata, min_counts=10)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    # set the leiden colors so the leiden pallets aren't grey
    color_path = os.getcwd().replace('\\', '/') + '/' + 'colors/'
    color_file = 'leiden_generated_random_3.txt'
    color_path += color_file
    with open(color_path, "r") as color_file:
        colors = color_file.read().split('\n')

    while len(colors) > 26:
        colors.pop()

    # you can choose either of these to run this with, you will get different colors
    # if you change the path above you can set your own .csv of hexes
    # adata.uns['leiden_colors'] = matplotlib.colors.CSS4_COLORS.values()
    adata.uns['leiden_colors'] = colors

    sc.pl.umap(
        adata,
        color=[
            "total_counts",
            "n_genes_by_counts",
            "leiden",
        ],
        wspace=0.4,
        show=False,
    )

    print("Visualize annotation on UMAP and spatial coordinates.")
    sq.pl.spatial_scatter(
        adata,
        shape=None,
        color=[
            "leiden",
        ],
        wspace=0.4,
    )

    sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

    sq.gr.centrality_scores(adata, cluster_key="leiden")

    sq.pl.centrality_scores(adata, cluster_key="leiden", figsize=(16, 5))

    adata_subsample = sc.pp.subsample(adata, fraction=0.5, copy=True)

    sq.gr.co_occurrence(
        adata_subsample,
        cluster_key="leiden",
    )
    sq.pl.co_occurrence(
        adata_subsample,
        cluster_key="leiden",
        clusters="12",
        figsize=(10, 10),
    )
    sq.pl.spatial_scatter(
        adata_subsample,
        color="leiden",
        shape=None,
        size=2,
    )

    sq.gr.nhood_enrichment(adata, cluster_key="leiden")

    fig, ax = plt.subplots(1, 2, figsize=(13, 7))
    sq.pl.nhood_enrichment(
        adata,
        cluster_key="leiden",
        figsize=(8, 8),
        title="Neighborhood enrichment adata",
        ax=ax[0],
    )
    sq.pl.spatial_scatter(adata_subsample, color="leiden", shape=None, size=2, ax=ax[1])

    # note: calculate ripley's statistics

    fig, ax = plt.subplots(1, 2, figsize=(15, 7))
    mode = "L"

    sq.gr.ripley(adata, cluster_key="leiden", mode=mode)
    sq.pl.ripley(adata, cluster_key="leiden", mode=mode, ax=ax[0])

    sq.pl.spatial_scatter(
        adata_subsample,
        color="leiden",
        groups=["0", "1", "3"],
        shape=None,
        size=2,
        ax=ax[1],
    )

    # sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)
    # sq.gr.spatial_autocorr(
    #     adata_subsample,
    #     mode="moran",
    #     n_perms=100,
    #     n_jobs=1,
    # )
    # adata_subsample.uns["moranI"].head(10)

    sq.pl.spatial_scatter(
        adata_subsample,
        color=[
            "Slc17a7",
            "Npy2r",
            "leiden"
        ],
        shape=None,
        size=2,
        img=False,
    )

    t1 = time.time()
    totalTime = t1 - t0
    print("time elapsed: {}".format(totalTime))

    plt.show()
# """


if __name__ == "__main__":
    esq_demo_1()
    # sq_demo_1()
