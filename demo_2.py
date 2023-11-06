# Sam Crouse
# scrouse2@uwyo.edu
# follows the squidpy vizgen liver vignette tutorial

# python imports
from copy import deepcopy
import pandas as pd
from scipy.cluster import hierarchy as sch
from matplotlib import pyplot as plt
import scanpy as sc
import squidpy as sq
import os
import time

# my imports
import EasySQ as esq
import EasySQTools as tools


def esq_demo_2():
    t0 = time.time()  # start timer

    # path = 'F:/sunlabmerfishdata/QSFL01222023/'
    path = os.getcwd() + '/tutorial_data_2/'
    esqAn = esq.Analysis(data_path=path)
    esqAn.print()

    esqAn.setLeidenColors("leiden_generated_random_5.txt")

    esqAn.qcMetrics(percentTop=(50, 100, 200, 300))  # check

    esqAn.filterCells(minCounts=50)
    esqAn.filterGenes(minCells=10)

    print("normalize total")
    esqAn.normalizeTotal()
    print("log transform")
    esqAn.log1p()
    print("scale")
    esqAn.scale(maxValue=10)

    resolution = 1.5  # 1.5
    print("PCA")
    esqAn.tl_pca()
    print("neighbors")
    esqAn.pp_neighbors(nNeighbors=10, nPcs=20)
    print("UMAP")
    esqAn.tl_umap()
    print("leiden")
    esqAn.leiden(resolution=resolution)

    esqAn.pl_umap(graphs=["leiden"], size=1)
    esqAn.availableGraphs()
    esqAn.spatialScatter(graphs=["leiden"], libraryID="spatial", figsize=(10, 10), size=0.5)

    print("Get and Assign Cell Types")
    esqAn.assignReferenceCells()

    print("Calculate Leiden Cluster Average Expression Signatures")
    esqAn.calculateLeidenExpressionSignatures()

    print("Assign Cell Type Based on Top Expressed Marker Genes")
    esqAn.assignCellTypesOnExpression()

    print("Neighborhood Enrichment")
    esqAn.gr_spatialNeighbors()
    esqAn.gr_nhoodEnrichment()
    esqAn.pl_nhoodEnrichment(vmin=-50, vmax=100)

    esqAn.showPlots()

    print("Neighborhood Enrichment Clusters")
    esqAn.plotNHoodEnrichmentClusters()

    print("Network Centrality Scores")
    esqAn.gr_centralityScores()

    esqAn.setDfCentral(deepcopy(esqAn.getAdata().uns["leiden_centrality_scores"]))
    esqAn.getDfCentral().index = esqAn.getMetaLeiden().index.tolist()

    # sort clusters based on centrality scores
    # closeness centrality - measure of how close the group is to other nodes.
    esqAn.calcSerCloseness()

    # degree centrality - fraction of non-group members connected to group members.
    esqAn.calcSerDegree()

    # clustering coefficient - measure of the degree to which nodes cluster together.
    esqAn.calcSerCluster()

    print("High Closeness Score")
    esqAn.highClosenessScore()

    print("Low Closeness Score")
    esqAn.lowClosenessScore()

    print("High Degree Centrality")
    esqAn.highDegreeCentrality()

    print("Low Degree Centrality")
    esqAn.lowDegreeCentrality()

    print("High Clustering Coefficient")
    esqAn.highClusteringCoefficient()

    print("Low Clustering Coefficient")
    esqAn.lowClusteringCoefficient()

    print("Autocorrelation: Moran's I Score")
    esqAn.calcMoransIScore(numView=12)
    esqAn.plotMoransIScore(size=0.5, figsize=(3, 3))

    t1 = time.time()  # end timer
    totalTime = t1 - t0
    print("time elapsed: {}".format(totalTime))

    esqAn.showPlots()
    return


def sq_demo_2():
    t0 = time.time()  # start timer

    # In [3]:
    print("\nImporting data...")
    vizgen_path = vizgen_path = os.getcwd() + '\\'
    meta_data_path = os.getcwd() + "/tutorial_data_2/" + "Liver1Slice1_cell_metadata.csv"
    cell_by_gene_path = os.getcwd() + "/tutorial_data_2/" + "Liver1Slice1_cell_by_gene.csv"
    transformation_file_path = 'micron_to_mosaic_pixel_transform.csv'

    adata = sq.read.vizgen(path=vizgen_path,
                           counts_file=cell_by_gene_path,
                           meta_file=meta_data_path,
                           transformation_file=transformation_file_path
                           )

    # path = os.getcwd() + '/tutorial_data_2/'
    # esqAn = esq.Analysis(data_path=path)
    # esqAn.print()
    # adata = esqAn.getAdata()

    print(adata)

    """
    # note: set the leiden colors so the leiden pallets aren't grey
    color_path = vizgen_path + "colors\\"
    color_file1 = 'leiden_color_set_1_gradient.csv'
    color_file2 = 'leiden_color_set_1_random.csv'
    color_file3 = 'leiden_color_set_2_random.csv'
    color_file4 = 'leiden_color_set_3_random.csv'
    color_file5 = 'leiden_generated_random_3.txt'

    color_files = [color_file1, color_file2, color_file3, color_file4, color_file5]
    color_data = []

    for file in color_files:
        current_color_path = color_path + file
        with open(current_color_path, "r") as color_file:
            colors = color_file.read().split('\n')

        while len(colors) > 32:
            colors.pop()

        color_data.append(colors)

    # you can choose either of these to run this with, you will get different colors
    # if you change the path above you can set your own .csv of hexes
    # adata.uns['leiden_colors'] = matplotlib.colors.CSS4_COLORS.values()
    adata.uns['leiden_colors'] = color_data[4]
    # adata.uns['random_colors_1'] = color_data[1]
    # adata.uns['random_colors_2'] = color_data[2]
    # adata.uns['random_colors_3'] = color_data[3]
    # """

    # In [4]:
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=(50, 100, 200, 300), inplace=True
    )

    # In [5]:
    sc.pp.filter_cells(adata, min_counts=50)
    sc.pp.filter_genes(adata, min_cells=10)

    # In [6]:
    print("normalize total")
    sc.pp.normalize_total(adata)
    print("log transform")
    sc.pp.log1p(adata)
    print("scale")
    sc.pp.scale(adata, max_value=10)

    # In [7]:
    resolution = 1.5
    print("PCA")
    sc.tl.pca(adata, svd_solver="arpack")
    print("neighbors")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    print("UMAP")
    sc.tl.umap(adata)
    print("Leiden")
    sc.tl.leiden(adata, resolution=resolution)

    # In [9]:
    # note: display graphs up to this point
    sc.set_figure_params(figsize=(10, 10))
    sc.pl.umap(adata, color=["leiden"], size=5)

    # In [10]:
    sq.pl.spatial_scatter(
        adata, shape=None, color="leiden", size=0.5, library_id="spatial", figsize=(10, 10)
    )

    # In [12]:
    gene_panel = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41421-021-00266-1/MediaObjects/41421_2021_266_MOESM1_ESM.xlsx"
    df_ref_panel_ini = pd.read_excel(gene_panel, index_col=0)
    df_ref_panel = df_ref_panel_ini.iloc[1:, :1]
    df_ref_panel.index.name = None
    df_ref_panel.columns = ["Function"]

    # Assign marker gene metadata using reference dataset
    marker_genes = df_ref_panel[
        df_ref_panel["Function"].str.contains("marker")
    ].index.tolist()

    meta_gene = deepcopy(adata.var)
    common_marker_genes = list(set(meta_gene.index.tolist()).intersection(marker_genes))
    meta_gene.loc[common_marker_genes, "Markers"] = df_ref_panel.loc[
        common_marker_genes, "Function"
    ]
    meta_gene["Markers"] = meta_gene["Markers"].apply(
        lambda x: "N.A." if "marker" not in str(x) else x
    )
    counts = meta_gene["Markers"].value_counts()
    print(counts)

    # In [13]:
    ser_counts = adata.obs["leiden"].value_counts()
    ser_counts.name = "cell counts"
    meta_leiden = pd.DataFrame(ser_counts)

    cat_name = "leiden"
    sig_leiden = pd.DataFrame(
        columns=adata.var_names, index=adata.obs[cat_name].cat.categories
    )
    for clust in adata.obs[cat_name].cat.categories:
        sig_leiden.loc[clust] = adata[adata.obs[cat_name].isin([clust]), :].X.mean(0)
    sig_leiden = sig_leiden.transpose()
    leiden_clusters = ["Leiden-" + str(x) for x in sig_leiden.columns.tolist()]
    sig_leiden.columns = leiden_clusters
    meta_leiden.index = sig_leiden.columns.tolist()
    meta_leiden["leiden"] = pd.Series(
        meta_leiden.index.tolist(), index=meta_leiden.index.tolist()
    )

    # In [14]:
    meta_gene = pd.DataFrame(index=sig_leiden.index.tolist())
    meta_gene["info"] = pd.Series("", index=meta_gene.index.tolist())
    meta_gene["Markers"] = pd.Series("N.A.", index=sig_leiden.index.tolist())
    meta_gene.loc[common_marker_genes, "Markers"] = df_ref_panel.loc[
        common_marker_genes, "Function"
    ]

    meta_leiden["Cell_Type"] = pd.Series("N.A.", index=meta_leiden.index.tolist())
    num_top_genes = 30
    for inst_cluster in sig_leiden.columns.tolist():
        top_genes = (
            sig_leiden[inst_cluster]
            .sort_values(ascending=False)
            .index.tolist()[:num_top_genes]
        )

        inst_ser = meta_gene.loc[top_genes, "Markers"]
        inst_ser = inst_ser[inst_ser != "N.A."]
        ser_counts = inst_ser.value_counts()

        max_count = ser_counts.max()

        max_cat = "_".join(sorted(ser_counts[ser_counts == max_count].index.tolist()))
        max_cat = max_cat.replace(" marker", "").replace(" ", "-")

        print(inst_cluster, max_cat)
        meta_leiden.loc[inst_cluster, "Cell_Type"] = max_cat

    # rename clusters
    meta_leiden["name"] = meta_leiden.apply(
        lambda x: x["Cell_Type"] + "_" + x["leiden"], axis=1
    )
    leiden_names = meta_leiden["name"].values.tolist()
    meta_leiden.index = leiden_names

    # transfer cell type labels to single cells
    leiden_to_cell_type = deepcopy(meta_leiden)
    leiden_to_cell_type.set_index("leiden", inplace=True)
    leiden_to_cell_type.index.name = None

    adata.obs["Cell_Type"] = adata.obs["leiden"].apply(
        lambda x: leiden_to_cell_type.loc["Leiden-" + str(x), "Cell_Type"]
    )
    adata.obs["Cluster"] = adata.obs["leiden"].apply(
        lambda x: leiden_to_cell_type.loc["Leiden-" + str(x), "name"]
    )

    sq.pl.spatial_scatter(
        adata,
        shape=None,
        # color=["Vwf", "Axin2"],
        color=["leiden"],
        size=15,
        # cmap="Reds",
        img=False,
        figsize=(12, 8),
    )

    # In [19]:
    sq.gr.spatial_neighbors(adata, coord_type="generic", spatial_key="spatial")
    sq.gr.nhood_enrichment(adata, cluster_key="leiden")
    sq.pl.nhood_enrichment(
        adata,
        cluster_key="leiden",
        method="average",
        cmap="inferno",
        vmin=-50,
        vmax=100,
        figsize=(5, 5),
    )

    # In [20]:
    n_clusters = [4]
    df_nhood_enr = pd.DataFrame(
        adata.uns["leiden_nhood_enrichment"]["zscore"],
        columns=leiden_clusters,
        index=leiden_clusters,
    )
    nhood_cluster_levels = ["Level-" + str(x) for x in n_clusters]
    linkage = sch.linkage(df_nhood_enr, method="average")
    mat_nhood_clusters = sch.cut_tree(linkage, n_clusters=n_clusters)
    df_cluster = pd.DataFrame(
        mat_nhood_clusters, columns=nhood_cluster_levels, index=meta_leiden.index.tolist()
    )

    inst_level = "Level-" + str(n_clusters[0])
    all_clusters = list(df_cluster[inst_level].unique())
    # sc.set_figure_params(figsize=(10,10))
    """
    for inst_cluster in all_clusters:
        inst_clusters = df_cluster[df_cluster[inst_level] == inst_cluster].index.tolist()
        sq.pl.spatial_scatter(
            adata,
            groups=inst_clusters,
            shape=None,
            color=["Th"],
            size=15,
            img=False,
            figsize=(12, 8),
        )

        # plt.show()
    # """

    sq.pl.spatial_scatter(
        adata,
        shape=None,
        color=["Cluster"],
        size=15,
        figsize=(12, 12),
        img=False,
    )

    sq.pl.spatial_scatter(
        adata,
        shape=None,
        color=["Cell_Type"],
        size=15,
        img=False,
    )

    sq.pl.spatial_scatter(
        adata,
        shape=None,
        color=["n_counts"],
        size=15,
        img=False,
    )

    # note: Network Centrality Scores

    # In [21]:
    sq.gr.centrality_scores(adata, "leiden")
    sc.set_figure_params(figsize=(20, 8))

    # copy centrality data to new DataFrame
    df_central = deepcopy(adata.uns["leiden_centrality_scores"])
    df_central.index = meta_leiden.index.tolist()

    # sort clusters based on centrality scores
    # closeness centrality - measure of how close the group is to other nodes.
    ser_closeness = df_central["closeness_centrality"].sort_values(ascending=False)

    # The degree centrality for a node v is the fraction of nodes it is connected to.
    ser_degree = df_central["degree_centrality"].sort_values(ascending=False)

    # clustering coefficient - measure of the degree to which nodes cluster together.
    ser_cluster = df_central["average_clustering"].sort_values(ascending=False)

    # In [22]:
    inst_clusters = ser_closeness.index.tolist()[:5]
    print(inst_clusters)

    sq.pl.spatial_scatter(
        adata,
        shape=None,
        color=["Cluster"],
        size=15,
        figsize=(12, 24),
        img=False
    )

    sq.pl.spatial_scatter(
        adata, shape=None, groups=inst_clusters, color=["Cluster"], size=15, img=False, figsize=(10, 10)
    )

    # In [23]:
    inst_clusters = ser_closeness.index.tolist()[-5:]
    print(inst_clusters)
    sq.pl.spatial_scatter(
        adata, shape=None, groups=inst_clusters, color="Cluster", size=15, img=False, figsize=(10, 10)
    )

    # In [24]:
    inst_clusters = ser_degree.index.tolist()[:5]
    print(inst_clusters)
    sq.pl.spatial_scatter(
        adata, shape=None, groups=inst_clusters, color="Cluster", size=15, img=False, figsize=(10, 10)
    )

    # In [25]:
    inst_clusters = ser_degree.index.tolist()[-5:]
    print(inst_clusters)
    sq.pl.spatial_scatter(
        adata, shape=None, groups=inst_clusters, color="Cluster", size=15, img=False, figsize=(10, 10)
    )

    # In [26]:
    inst_clusters = ser_cluster.index.tolist()[:5]
    print(inst_clusters)
    sq.pl.spatial_scatter(
        adata, shape=None, groups=inst_clusters, color="Cluster", size=15, img=False, figsize=(10, 10)
    )

    # In [27]:
    inst_clusters = ser_cluster.index.tolist()[-5:]
    print(inst_clusters)
    sq.pl.spatial_scatter(
        adata, shape=None, groups=inst_clusters, color="Cluster", size=15, img=False, figsize=(15, 15)
    )

    # In [28]:
    sq.gr.spatial_autocorr(adata, mode="moran")
    num_view = 12
    top_autocorr = (
        adata.uns["moranI"]["I"].sort_values(ascending=False).head(num_view).index.tolist()
    )
    bot_autocorr = (
        adata.uns["moranI"]["I"].sort_values(ascending=True).head(num_view).index.tolist()
    )

    # In [34]:
    sq.pl.spatial_scatter(
        adata, shape=None, color=top_autocorr, size=20, cmap="Reds", img=False, figsize=(5, 5)
    )

    # In [30]:
    # top cell types based on average expression of top_autocorr genes
    top_cell_types = sig_leiden.loc[top_autocorr].mean(axis=0).sort_values(ascending=False).index.tolist()[:5]
    print(top_cell_types)


    # In [35]:
    sq.pl.spatial_scatter(
        adata, shape=None, color=bot_autocorr, size=20, cmap="Reds", img=False, figsize=(5, 5)
    )

    # top cell types based on average expression of bot_autocorr genes
    # In [32]:
    top_cell_types = sig_leiden.loc[bot_autocorr].mean(axis=0).sort_values(ascending=False).index.tolist()[:5]
    print(top_cell_types)

    t1 = time.time()  # end timer
    totalTime = t1 - t0
    print("time elapsed: {}".format(totalTime))

    plt.show()
    return


if __name__ == "__main__":
    esq_demo_2()
    # sq_demo_2()
