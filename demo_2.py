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
import EasySQ as eSQ


def esq_demo_2():
    t0 = time.time()  # start timer

    path = 'F:/sunlabmerfishdata/QSFL01222023/'
    esq = eSQ.Analysis(data_path=path)
    esq.print()

    esq.qcMetrics(percentTop=(50, 100, 140))

    esq.filterCells()
    esq.filterGenes()

    print("normalize total")
    esq.normalizeTotal()
    print("log transform")
    esq.log1p()
    print("scale")
    esq.scale()

    resolution = 1.5
    print("PCA")
    esq.tl_pca()
    print("neighbors")
    esq.pp_neighbors()
    print("UMAP")
    esq.tl_umap()
    print("leiden")
    esq.leiden(resolution=resolution)

    esq.pl_umap(graphs=["leiden"])

    esq.spatialScatter(graphs="leiden", libraryID="spatial")

    print("assign cell types")
    esq.assignReferenceCells()

    print("Calculate Leiden Cluster Average Expression Signatures")
    esq.calculateLeidenExpressionSignatures()

    print("Assign Cell Type Based on Top Expressed Marker Genes")
    esq.assignCellTypesOnExpression()

    print("Neighborhood Enrichment")
    esq.gr_spatialNeighbors()
    esq.gr_nhoodEnrichment()
    esq.pl_nhoodEnrichment()

    print("Neighborhood Enrichment Clusters")
    n_clusters = [4]
    df_nhood_enr = pd.DataFrame(
        esq.getAdata().uns["leiden_nhood_enrichment"]["zscore"],
        columns=esq.getLeidenClusters(),
        index=esq.getLeidenClusters(),
    )

    nhood_cluster_levels = ["Level-" + str(x) for x in n_clusters]
    linkage = sch.linkage(df_nhood_enr, method="average")
    mat_nhood_clusters = sch.cut_tree(linkage, n_clusters=n_clusters)
    df_cluster = pd.DataFrame(
        mat_nhood_clusters, columns=nhood_cluster_levels, index=esq.getMetaLeiden().index.tolist()
    )

    inst_level = "Level-" + str(n_clusters[0])
    all_clusters = list(df_cluster[inst_level].unique())
    # sc.set_figure_params(figsize=(10,10))
    for inst_cluster in all_clusters:
        inst_clusters = df_cluster[df_cluster[inst_level] == inst_cluster].index.tolist()

        esq.spatialScatter(graphs="Cluster", groups=inst_clusters, size=15)

    print("Network Centrality Scores")
    esq.gr_centrality_scores()

    # copy centrality data to new DataFrame
    # df_central = deepcopy(esq.getAdata().uns["leiden_centrality_scores"])
    # df_central.index = esq.getMetaLeiden().index.tolist()

    esq.setDfCentral(deepcopy(esq.getAdata().uns["leiden_centrality_scores"]))
    esq.getDfCentral().index = esq.getMetaLeiden().index.tolist()

    # sort clusters based on centrality scores
    ################################################
    # closeness centrality - measure of how close the group is to other nodes.
    esq.calcSerCloseness()

    # degree centrality - fraction of non-group members connected to group members.
    esq.calcSerDegree()

    # clustering coefficient - measure of the degree to which nodes cluster together.
    esq.calcSerCluster()

    print("High Closeness Score")
    esq.highClosenessScore()

    print("Low Closeness Score")
    esq.lowClosenessScore()

    print("High Degree Centrality")
    esq.highDegreeCentrality()

    print("Low Degree Centrality")
    esq.lowDegreeCentrality()

    print("High Clustering Coefficient")
    esq.highClusteringCoefficient()

    print("Low Clustering Coefficient")
    esq.lowClusteringCoefficient()

    print("Autocorrelation: Moran's I Score")
    sq.gr.spatial_autocorr(esq.getAdata(), mode="moran")
    num_view = 12
    top_autocorr = (
        esq.getAdata().uns["moranI"]["I"].sort_values(ascending=False).head(num_view).index.tolist()
    )
    bot_autocorr = (
        esq.getAdata().uns["moranI"]["I"].sort_values(ascending=True).head(num_view).index.tolist()
    )
    esq.spatialScatter(graphs=top_autocorr, cmap="Reds", size=10, figsize=(3, 3))

    # top cell types based on average expression of top_autocorr genes
    print(esq.getSigLeiden().loc[top_autocorr].mean(axis=0).sort_values(ascending=False).index.tolist()[
    :5
    ])

    print("Genes with low autocorrelation")
    esq.spatialScatter(graphs=bot_autocorr, cmap="Reds",  size=10, figsize=(3, 3))

    # top cell types based on average expression of bot_autocorr genes
    print(esq.getSigLeiden().loc[bot_autocorr].mean(axis=0).sort_values(ascending=False).index.tolist()[
    :5
    ])

    t1 = time.time()  # end timer
    totalTime = t1 - t0
    print("time elapsed: {}".format(totalTime))

    esq.showPlots()
    return


def sq_demo_2():
    ####################################################################################################################
    ####################################################################################################################
    # note: Single-Cell Clustering of Vizgen MERFISH Mouse Liver Data
    #  Obtain Data from Vizgen

    # note: We will use the Liver1Slice1 dataset from Vizgen's MERFISH Mouse Liver Map:
    #  https://info.vizgen.com/mouse-liver-access. In order to run this tutorial we will download the cell_by_gene.csv
    #  and meta_cell.csv. Please follow the instructions to obtain access to the showcase data and download the
    #  data - here we save the data to a directory called tutorial_data/ in the same directory as this notebook.
    # note: create an adata object with the Vizgen liver data
    # In [3]:
    print("\nImporting data...")
    vizgen_path = vizgen_path = os.getcwd() + '\\'
    meta_data_path = 'F:\\sunlabmerfishdata\\QSFL01222023\\cell_metadata.csv'
    cell_by_gene_path = 'F:\\sunlabmerfishdata\\QSFL01222023\\cell_by_gene.csv'
    transformation_file_path = 'micron_to_mosaic_pixel_transform.csv'

    adata = sq.read.vizgen(path=vizgen_path,
                           counts_file=cell_by_gene_path,
                           meta_file=meta_data_path,
                           transformation_file=transformation_file_path
                           )
    print(adata)

    ####################################################################################################################
    ####################################################################################################################
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

        while len(colors) > 33:
            colors.pop()

        color_data.append(colors)

    # you can choose either of these to run this with, you will get different colors
    # if you change the path above you can set your own .csv of hexes
    # adata.uns['leiden_colors'] = matplotlib.colors.CSS4_COLORS.values()
    adata.uns['leiden_colors'] = color_data[4]
    adata.uns['random_colors_1'] = color_data[1]
    adata.uns['random_colors_2'] = color_data[2]
    adata.uns['random_colors_3'] = color_data[3]

    ####################################################################################################################
    ####################################################################################################################
    # note: Make gene names unique and calculate QC metrics

    # In [4]:
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=(50, 100), inplace=True
    )

    ####################################################################################################################
    ####################################################################################################################
    # note: Filter cells with low expression and genes that are expressed in too few cells.

    # In [5]:
    sc.pp.filter_cells(adata, min_counts=50)
    sc.pp.filter_genes(adata, min_cells=10)

    ####################################################################################################################
    ####################################################################################################################
    # note: Data Pre-processing
    # note: Here we use Scanpy total-count normalize, logarithmize, and scale gene expression to unit
    #  variance (clipping values that exceed 10 standard deviations).

    # In [6]:
    print("normalize total")
    sc.pp.normalize_total(adata)
    print("log transform")
    sc.pp.log1p(adata)
    print("scale")
    sc.pp.scale(adata, max_value=10)

    ####################################################################################################################
    ####################################################################################################################
    # note: Dimensionality Reduction, Neighbor Calculation, and Clustering
    #  Here we use Scanpy to: reduce the dimensionality of our data by running principal component analysis
    #  (PCA{, calculate the neighborhood graph of cells in PCA space, display cells in a two-dimensional UMAP embedding,
    #  and finally identify clusters of cells using Leiden clustering.

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

    ####################################################################################################################
    ####################################################################################################################
    # note: UMAP with Leiden Clustering Labels
    #  Here we visualize the distributions of the Leiden clusters in the UMAP plot. We see Leiden clusters tend to
    #  segregate into distinct regions within the UMAP plot.

    # In [9]:
    # note: display graphs up to this point
    # sc.set_figure_params(figsize=(10, 10))
    # sc.pl.umap(adata, color=["leiden"], size=5)

    ####################################################################################################################
    ####################################################################################################################
    # note: Spatial Distributions of Cells
    #  Here we visualize the spatial locations of cells in the mouse liver colored by Leiden cluster. We observe
    #  distinct spatial localizations of Leiden clusters throughout the tissue. We see some clusters line blood vessels
    #  and while others form concentric patterns reflecting hepatic zonation (Cunningham et. al. 2021). Our next step is
    #  to assign tentative cell type to our Leiden clusters and assess their spatial localizations in the liver.

    # In [10]:
    sq.pl.spatial_scatter(
        adata, shape=None, color="leiden", size=0.5, library_id="spatial", figsize=(10, 10)
    )

    ####################################################################################################################
    ####################################################################################################################
    # note: Assign Cell Types
    #  Reference Cell Type Marker Gene Sets
    #  In order to tentatively assign liver cell types we utilize a gene-level cell type marker reference from the
    #  publication Spatial transcriptome profiling by MERFISH reveals fetal liver hematopoietic stem cell niche
    #  architecture. These marker genes will be used to assess cell type composition of the Leiden clusters. See
    #  MERFISH gene panel metadata: https://www.nature.com/articles/s41421-021-00266-1.

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

    ####################################################################################################################
    ####################################################################################################################
    # note: Calculate Leiden Cluster Average Expression Signatures
    #  Here we calculate the average gene expression signatures of the Leiden clusters, which will be used to assign
    #  cell type composition.

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

    ####################################################################################################################
    ####################################################################################################################
    # note: Assign Cell Type Based on Top Expressed Marker Genes
    #  Here we assign cell type composition to the Leiden clusters by counting the frequency of cell type marker genes
    #  in the top 30 most up-regulated genes for each cluster such that cell type is assigned based on the most
    #  frequently occuring marker genes. If there is a tie in the number of marker genes, we assign the cluster to more
    #  than one cell type. Using this approach eleven Leiden clusters are assigned to be Hepatocyte containing clusters.

    print("\n")

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

    ####################################################################################################################
    ####################################################################################################################
    # note: Hepatocyte Zonation
    #  Hepatocytes are the most abundant cell in the liver and have multiple roles in metabolism, endocrine production,
    #  protein synthesis, and detoxification. Hepatocytes form complex, radial structures called lobules that contain a
    #  central vein (with low blood oxygen level) surrounded by peripheral portal veins (with high blood oxygen level).
    #  Hepatocytes can also be broadly classified as peri-central or peri-portal based on their proximity to central and
    #  portal veins, respectively.

    # note Central and Portal Blood Vessels
    #  We can use the gene Vwf to identify endothelial cells (Horvath et. al. 2004) that line liver blood vessels.
    #  Plotting Vwf expression level in single cells (below) shows clear enrichment at blood vessel borders.
    # Note that Vwf expression can also be used to distinguish artifactual holes/tears in the liver from blood vessels.
    #  Next, we use the expression of Axin2 to mark peri-central regions (Sun et. al. 2020). Plotting Axin2 expression
    #  level in single cells shows Axin2 lining a subset of all Vwf positive blood vessels, which allows us to
    #  distinguish peri-portal (Axin2 negative) and peri-central (Axin2 positive) blood vessels. Similarly, we can use
    #  Axin2 expression to identify peri-central hepatocyte Leiden clusters (see next section).

    # note: display graphs up to this point
    # sc.set_figure_params(figsize=(10, 10))
    # sc.pl.umap(adata, color=["leiden"], size=5)

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

    ####################################################################################################################
    ####################################################################################################################
    # note: Distinguishing Peri-Portal and Peri-Central Hepatocytes
    #  As described above, we use the expression of Axin2 as a marker for peri-central hepatocytes

    # In [16]:
    # all_hepatocyte_clusters = [x for x in meta_leiden.index.tolist() if "Hepatocyte" in x]
    # sig_leiden.columns = meta_leiden.index.tolist()
    # ser_axin2 = sig_leiden[all_hepatocyte_clusters].loc["Axin2"]
    # peri_central = ser_axin2[ser_axin2 > 0].index.tolist()
    # peri_portal = ser_axin2[ser_axin2 <= 0].index.tolist()

    ####################################################################################################################
    ####################################################################################################################
    # note: Peri-Central Hepatocytes
    #  Plotting peri-central and peri-portal hepatocytes separately displays their distinct interlocking morphologies
    #  with respect to blood vessels.

    # In [33]:
    # sq.pl.spatial_scatter(
    #     adata, groups=peri_central, color="Cluster", size=15, img=False, figsize=(15, 15)
    # )

    ####################################################################################################################
    ####################################################################################################################
    # note: Peri-Portal Hepatocytes

    # In [18]:
    # sq.pl.spatial_scatter(
    #     adata, groups=peri_portal, color="Cluster", size=15, img=False, figsize=(15, 15)
    # )

    ####################################################################################################################
    ####################################################################################################################
    # note: Neighborhood Enrichment
    #  In this section we will use Squidpy to identify clusters that are spatially enriched for one another using a
    #  neighborhood enrichment test {func}squidpy.gr.nhood_enrichment. This test determines if cells belonging to two
    #  different clusters are close to each other more often than expected.

    # note: In order to run this test we first have to
    #  calculate a connectivity graph using the {func}squidpy.gr.spatial_neighbors method. This graph consists of cells
    #  (nodes) and cell-cell interactions (edges).

    # note: We also visualize the neighborhood enrichment using a hierarchically
    #  clustered heatmap which shows clusters of enriched neighborhoods in our tissue.

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

    ####################################################################################################################
    ####################################################################################################################
    # note: Neighborhood Enrichment Clusters
    #  Here we visualize clusters from our neighborhood enrichment data obtained by hierarchically clustering the
    #  Z-scored neighborhood enrichment scores. We observe a cluster containing three peri-portal hepatocyte leiden
    #  clusters (Hepatocyte_Leiden-0, Hepatocyte_Leiden-4, and Hepatocyte_Leiden-6 ) and a large cluster containing
    #  several peri-central hepatocytes (Hepatocyte_Leiden-1, Hepatocyte_Leiden-2, Hepatocyte_Leiden-25). This
    #  demonstrates that we can recapitulate known spatial enrichment of peri-portal and peri-central hepatocytes using
    #  neighborhood enrichment.

    print("In [20]")

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
    # """
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

    var_list = list(adata.var_names)
    obs_columns_list = list(adata.obs.columns)

    print(var_list)
    print(obs_columns_list)

    # stand-outs: Cluster, Cell_Type, leiden, n_counts, pct_counts_in_top_50_genes, log1p_total_counts,
    # log1p_n_genes_by_counts, max_y, min_y, max_x, min_x, volume, fov,

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

    ####################################################################################################################
    ####################################################################################################################
    # note: Network Centrality Scores

    # In [21]:
    sq.gr.centrality_scores(adata, "leiden")
    sc.set_figure_params(figsize=(20, 8))

    # copy centrality data to new DataFrame
    df_central = deepcopy(adata.uns["leiden_centrality_scores"])
    df_central.index = meta_leiden.index.tolist()

    # sort clusters based on centrality scores
    ################################################
    # closeness centrality - measure of how close the group is to other nodes.
    ser_closeness = df_central["closeness_centrality"].sort_values(ascending=False)

    # degree centrality - fraction of non-group members connected to group members.
    # [Networkx](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.degree_centrality.html#networkx.algorithms.centrality.degree_centrality)
    # The degree centrality for a node v is the fraction of nodes it is connected to.
    ser_degree = df_central["degree_centrality"].sort_values(ascending=False)

    # clustering coefficient - measure of the degree to which nodes cluster together.
    ser_cluster = df_central["average_clustering"].sort_values(ascending=False)

    ####################################################################################################################
    ####################################################################################################################
    # note: High Closeness Score
    #  Groups/clusters with high closeness are close to other groups and will tend to display a dispersed distribution
    #  throughout the tissue. Three of the top five clusters based on the closeness score are epithelial cells
    #  (SEC: sinusoidal epithelial cells, AEC: arterial epithelial cells) and one cluster consists of macrophages.
    #  We see that these groups are indeed evenly distributed across the tissue which results in them being relatively
    #  close to many other groups.

    # works to here

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

    ####################################################################################################################
    ####################################################################################################################
    # note: Low Closeness Score
    #  Groups with low closeness are not close to other groups and tend to display an uneven and isolated distribution
    #  throughout the tissue. We see that clusters with low closeness scores tend to be located near blood vessels and
    #  consist of megakaryocyte, neutrophil, and erythroid cells. Their distinct localization and proximity to blood
    #  vessel (e.g., only neighboring cells on one side while the other side faces the blood vessel) contributes to
    #  their low level of interactions with other clusters.

    # In [23]:
    inst_clusters = ser_closeness.index.tolist()[-5:]
    print(inst_clusters)
    sq.pl.spatial_scatter(
        adata, shape=None, groups=inst_clusters, color="Cluster", size=15, img=False, figsize=(10, 10)
    )

    ####################################################################################################################
    ####################################################################################################################
    # note: High Degree Centrality
    #  Similarly to the results from the closeness scores we observed above, we see that the SEC and Macrophage
    #  clusters (Leiden-3, Leiden-5, and Leiden-8) have high degree centrality scores indicating that they have a high
    #  fraction of non-group member connections. We also n*te that these clusters tend to be more evenly distributed
    #  throughout the tissue.

    # In [24]:
    inst_clusters = ser_degree.index.tolist()[:5]
    print(inst_clusters)
    sq.pl.spatial_scatter(
        adata, shape=None, groups=inst_clusters, color="Cluster", size=15, img=False, figsize=(10, 10)
    )

    ####################################################################################################################
    ####################################################################################################################
    # note: Low Degree Centrality
    #  These cluster have particularly low non-group member connections and similarly to the results from the closeness
    #  score we see that these clusters tend to be near blood vessels. We also n*te that more of the clusters tend to
    #  be lower abundance clusters (Leiden cluster numbers are ranked by abundance).

    # In [25]:
    inst_clusters = ser_degree.index.tolist()[-5:]
    print(inst_clusters)
    sq.pl.spatial_scatter(
        adata, shape=None, groups=inst_clusters, color="Cluster", size=15, img=False, figsize=(10, 10)
    )

    ####################################################################################################################
    ####################################################################################################################
    # note: High Clustering Coefficient
    #  The clustering coefficient indicates the degree to which nodes cluster together. We see that the top scoring
    #  clusters tend to be located near blood vessels. This distribution is somewhat counter-intuitive since one might
    #  expect that well segregated cells would form dense blobs rather than thin lines. However, their position along
    #  the edge of of blood vessels likely reduces the total number of neighbors each cell has since they are positioned
    #  at the border of a hole in the tissue.

    # In [26]:
    inst_clusters = ser_cluster.index.tolist()[:5]
    print(inst_clusters)
    sq.pl.spatial_scatter(
        adata, shape=None, groups=inst_clusters, color="Cluster", size=15, img=False, figsize=(10, 10)
    )

    ####################################################################################################################
    ####################################################################################################################
    # note: Low Clustering Coefficient
    #  Again, we see that non-isolated groups tend to be more evenly distributed throughout the tissue. Interestingly,
    #  these clusters mostly consist of hepatocytes.

    # In [27]:
    inst_clusters = ser_cluster.index.tolist()[-5:]
    print(inst_clusters)
    sq.pl.spatial_scatter(
        adata, shape=None, groups=inst_clusters, color="Cluster", size=15, img=False, figsize=(15, 15)
    )

    ####################################################################################################################
    ####################################################################################################################
    # note: Auto-correlation: Moran's I Score
    #  Our previous focus has been mainly on the distribution of cell clusters throughout the tissue. However we can
    #  also use Squidpy to investigate the spatial distributions of genes expressed in the tissue.

    # note: Here we use Squidpy to calculate the Moran's I global spatial auto-correlation statistic, which can be used
    #  to identify genes that are non-randomly distributed in the tissue. We will visualize the top and bottom 20
    #  scoring genes to highlight specific examples of genes with high and low auto-correlation.

    # In [28]:
    sq.gr.spatial_autocorr(adata, mode="moran")
    num_view = 12
    top_autocorr = (
        adata.uns["moranI"]["I"].sort_values(ascending=False).head(num_view).index.tolist()
    )
    bot_autocorr = (
        adata.uns["moranI"]["I"].sort_values(ascending=True).head(num_view).index.tolist()
    )

    ####################################################################################################################
    ####################################################################################################################
    # note: Genes with high spatial auto-correlation
    #  We n*te that many of the top scoring genes show expression patterns following the spatial pattern of hepatocyte
    #   Leiden clusters (e.g. Aldh1b1, Aldh3a20) and blood vessels (e.g. Dpt, Sfrp1). We can also check which cell types
    #   are most associated with these highly auto-correlated genes by ranking cell clusters based on their expression
    #   of these genes (e.g. mean expression). Below we see that three of the top five cell clusters are hepatocytes,
    #   which agrees with the gene's localization patterns. These results indicate the primary spatial patterns in gene
    #   expression auto-correlation are being driven by their expression in hepatocytes and blood vessel associated
    #   cell types.

    # In [34]:
    sq.pl.spatial_scatter(
        adata, shape=None, color=top_autocorr, size=20, cmap="Reds", img=False, figsize=(5, 5)
    )

    # In [30]:
    # top cell types based on average expression of top_autocorr genes
    top_cell_types = sig_leiden.loc[top_autocorr].mean(axis=0).sort_values(ascending=False).index.tolist()[:5]
    print(top_cell_types)

    ####################################################################################################################
    ####################################################################################################################
    # note: Genes with low auto-correlation
    #  Genes with low auto-correlation show more evenly distributed expression patterns that do not follow
    #  hepatic zonation.

    # In [35]:
    sq.pl.spatial_scatter(
        adata, shape=None, color=bot_autocorr, size=20, cmap="Reds", img=False, figsize=(5, 5)
    )

    # top cell types based on average expression of bot_autocorr genes
    # In [32]:
    top_cell_types = sig_leiden.loc[bot_autocorr].mean(axis=0).sort_values(ascending=False).index.tolist()[:5]
    print(top_cell_types)

    ####################################################################################################################
    ####################################################################################################################
    # note: Conclusion
    #  This tutorial shows how we can us Sqidpy and Scanpy to explore the spatial distribution of single-cells and gene
    #  expression in mouse liver tissue from Vizgen's MERFISH Mouse Liver Map. One of the main highlights is the
    #  intricate hepatic zonation patterns and their relationship to portal/central veins in the mouse liver.

    plt.show()


if __name__ == "__main__":
    esq_demo_2()
    # sq_demo_2()
