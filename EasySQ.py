# Sam Crouse
# scrouse2@uwyo.edu

# Class which implements the tools built in EasySQTools.py
# Allows creation of a squidpy analysis object with built in tools

# python imports
from copy import deepcopy
import pandas as pd
from scipy.cluster import hierarchy as sch

# my imports
import EasySQTools as tools
import time


class Analysis:
    """
    # note: some rules for data
    # 1. if a var isn't set, it should be None
    # 2. all values initialized in constructor. If nothing is passed, set to None.
    # 3. object.variable should not be used. Variables should only be accessed through getter and setter functions.
    #      Except in constructor.
    # """

    ####################################################################################################################
    ####################################################################################################################
    # region constructor
    # note: constructor
    def __init__(self, data_path):
        self.ID = tools.createCode(9)

        # passed value initialization
        self.data_path = data_path
        self.adata = self.readVizgen(self.data_path)

        # default value init
        self.metaGene = None
        self.commonMarkerGenes = None
        self.dfRefPanel = None
        self.sigLeiden = None
        self.metaLeiden = None
        self.leidenClusters = None

        self.dfCentral = None
        self.serCloseness = None
        self.serDegree = None
        self.serCluster = None

    # endregion constructor

    ####################################################################################################################
    ####################################################################################################################
    # region setter functions
    # note: primary setter funcs

    # set the ID of this object
    def setID(self, ID):
        self.ID = str(ID)

    def setDataPath(self, data_path):
        self.data_path = data_path

    def setAdata(self, adata):
        self.adata = adata

    def setMetaGene(self, value):
        self.metaGene = value

    def setCommonMarkerGenes(self, value):
        self.commonMarkerGenes = value

    def setDfRefPanel(self, value):
        self.dfRefPanel = value

    def setSigLeiden(self, value):
        self.sigLeiden = value

    def setMetaLeiden(self, value):
        self.metaLeiden = value

    def setLeidenClusters(self, value):
        self.leidenClusters = value

    def setDfCentral(self, value):
        self.dfCentral = value

    def setSerCloseness(self, value):
        self.serCloseness = value

    def setSerDegree(self, value):
        self.serDegree = value

    def setSerCluster(self, value):
        self.serCluster = value

    # endregion

    ####################################################################################################################
    ####################################################################################################################
    # region getter functions
    # note: primary getter funcs
    def getID(self):
        return self.ID

    def getDataPath(self):
        return self.data_path

    def getAdata(self):
        return self.adata

    def getMetaGene(self):
        return self.metaGene

    def getCommonMarkerGenes(self):
        return self.commonMarkerGenes

    def getDfRefPanel(self):
        return self.dfRefPanel

    def getSigLeiden(self):
        return self.sigLeiden

    def getMetaLeiden(self):
        return self.metaLeiden

    def getLeidenClusters(self):
        return self.leidenClusters

    def getDfCentral(self):
        return self.dfCentral

    def getSerCloseness(self):
        return self.serCloseness

    def getSerDegree(self):
        return self.serDegree

    def getSerCluster(self):
        return self.serCluster

    # endregion

    ####################################################################################################################
    ####################################################################################################################
    # region functions
    # note: given a filepath, reads and creates an anndata object
    def readVizgen(self, data_path):
        return tools.readVizgen(data_path)

    def availableGraphs(self, log=True):
        return tools.availableGraphs(self.getAdata(), log=log)

    # plot a spatial scatter. use colorInit to generate a custom color palette by default
    def spatialScatter(self, graphs, adata=None, show=False, libraryID=None, wspace=0.4, size=None,
                       shape=None, groups=None, cmap=None, figsize=None):
        if adata is None:
            return tools.spatialScatter(adata=self.getAdata(), graphs=graphs, show=show, libraryID=libraryID,
                                        wspace=wspace, size=size, shape=shape, groups=groups, cmap=cmap,
                                        figsize=figsize)

        else:
            return tools.spatialScatter(adata=adata, graphs=graphs, show=show, libraryID=libraryID, wspace=wspace,
                                        size=size, shape=shape, groups=groups, cmap=cmap, figsize=figsize)

    # runs the adata setup and basic analysis tools for the user
    def adataSetup(self):
        return tools.adataSetup(self.getAdata())

    # adata functions
    # calculate qc metrics
    def qcMetrics(self, percentTop=(50, 100)):
        """
            Documentation

            qcMetrics() class function docs
            =============================
            ----

            * qcMetrics class function for EasySQ class.
            * Acts as a wrapper for EasySQTools qcMetrics function.
            * Calculates a number of qc metrics for the self AnnData object.

            ----

            Parameters:
            ----------
            * percentTop=(50, 100): Calculate the top % of genes. This cannot exceed size of adata columns.

            Returns:
            -------
            * float(perUnassigned) percent unassigned genes/cells.

            Examples:
            --------
            1. using percentTop:
                * Given an AnnData object with n_obs × n_vars = 78329 × 483.
                    percentTop could be equal to (1, X0, X1, 483)
                * For example, say percentTop=(1, 50, 483)
                * If percentTop is greater than 483, it will throw an IndexError.
            ..
            2. Calling qcMetrics:
                * Create EasySQ class object and use:
                * perUn = esq.qcMetrics()
                    Where 'perUn' is the perUnassigned returned and default args are used.
        """

        return tools.qcMetrics(self.getAdata(), percentTop)

    def layers(self):
        return tools.layers(self.getAdata())

    def highlyVariableGenes(self, nTopGenes=4000):
        return tools.highlyVariableGenes(self.getAdata(), nTopGenes)

    def normalizeTotal(self, inplace=True):
        return tools.normalizeTotal(self.getAdata(), inplace=inplace)

    def log1p(self):
        return tools.log1p(self.getAdata())

    def pp_pca(self):
        return tools.pp_pca(self.getAdata())

    def tl_pca(self, svdSolver="arpack"):
        return tools.tl_pca(self.getAdata(), svdSolver=svdSolver)

    def pp_neighbors(self):
        return tools.pp_neighbors(self.getAdata())

    # uses the sc.tl.umap function to calculate umap data
    def tl_umap(self):
        return tools.tl_umap(self.getAdata())

    def clustering(self):
        pass

    def leiden(self, resolution=1.0):
        return tools.leiden(self.getAdata(), resolution=resolution)

    # can be used to filter cells with low expression
    def filterCells(self, minCounts=50):
        return tools.filterCells(self.getAdata(), minCounts=minCounts)

    # can be used to filter genes that are expressed in too few cells
    def filterGenes(self, minCells=10):
        return tools.filterGenes(self.getAdata(), minCells=minCells)

    # can be used to scale gene expression. IE Clip values that exceed 10 ('max value') standard deviations
    def scale(self, maxValue=10):
        return tools.scale(self.getAdata(), maxValue=maxValue)

    def plotTranscripts(self, show=False, figsize=(15, 4)):
        return tools.plotTranscripts(self.getAdata(), show=show, figsize=figsize)

    # uses the sc.pl.umap function plot a umap. Use colorInit to generate a custom color palette by default
    def pl_umap(self, graphs=["leiden"], show=False, colorInit=False, size=None, wspace=0.4):
        return tools.pl_umap(self.getAdata(), graphs=graphs, show=show, size=size, wspace=wspace)

    def gr_spatialNeighbors(self, adata=None, coordType="generic", spatialKey="spatial", delaunay=False):
        if adata is None:
            return tools.gr_spatialNeighbors(adata=self.getAdata(), coordType=coordType, spatialKey=spatialKey,
                                             delaunay=delaunay)

        else:
            return tools.gr_spatialNeighbors(adata=adata, coordType=coordType, spatialKey=spatialKey, delaunay=delaunay)

    # calculate nhoodEnrichment
    def gr_nhoodEnrichment(self, clusterKey="leiden"):
        return tools.gr_nhoodEnrichment(adata=self.getAdata(), clusterKey=clusterKey)

    # plot nhoodEnrichment data
    def pl_nhoodEnrichment(self, show=False, clusterKey="leiden", method="average", cmap="inferno", vmin=-50, vmax=100,
                           figsize=(5, 5)):
        tools.pl_nhoodEnrichment(adata=self.getAdata(), show=show, clusterKey=clusterKey, method=method,
                                 cmap=cmap, vmin=vmin, vmax=vmax, figsize=figsize)

    # calculate neighborhood enrichment clusters
    def plotNHoodEnrichmentClusters(self, nClusters=[4], method="average", size=10, figsize=(5, 5)):
        n_clusters = nClusters
        df_nhood_enr = pd.DataFrame(
            self.getAdata().uns["leiden_nhood_enrichment"]["zscore"],
            columns=self.getLeidenClusters(),
            index=self.getLeidenClusters(),
        )

        nhood_cluster_levels = ["Level-" + str(x) for x in n_clusters]
        linkage = sch.linkage(df_nhood_enr, method=method)
        mat_nhood_clusters = sch.cut_tree(linkage, n_clusters=n_clusters)
        df_cluster = pd.DataFrame(
            mat_nhood_clusters, columns=nhood_cluster_levels, index=self.getMetaLeiden().index.tolist()
        )

        inst_level = "Level-" + str(n_clusters[0])
        all_clusters = list(df_cluster[inst_level].unique())
        for inst_cluster in all_clusters:
            inst_clusters = df_cluster[df_cluster[inst_level] == inst_cluster].index.tolist()

            self.spatialScatter(graphs="Cluster", groups=inst_clusters, size=size, figsize=figsize)

    # compute the centrality scores
    def gr_centrality_scores(self, cluster_key="leiden"):
        return tools.gr_centrality_scores(adata=self.getAdata(), cluster_key=cluster_key)

    # graph the centrality scores
    def pl_centrality_scores(self, cluster_key="leiden", figsize=None):
        return tools.pl_centrality_scores(adata=self.getAdata(), cluster_key=cluster_key, figsize=figsize)

    # compute the co-occurrence probability
    def gr_co_occurrence(self, adata=None, cluster_key="leiden"):
        if adata is None:
            return tools.gr_co_occurrence(adata=self.getAdata(), cluster_key=cluster_key)

        else:
            return tools.gr_co_occurrence(adata=adata, cluster_key=cluster_key)

    # graph the co-occurrence probability
    def pl_co_occurrence(self, adata=None, cluster_key="leiden", clusters="12", figsize=None):
        if adata is None:
            return tools.pl_co_occurrence(adata=self.getAdata(), cluster_key=cluster_key, clusters=clusters,
                                          figsize=figsize)

        else:
            return tools.pl_co_occurrence(adata=adata, cluster_key=cluster_key, clusters=clusters, figsize=figsize)

    # compute Ripley's statistics
    def gr_ripley(self, cluster_key="leiden", mode="L"):
        return tools.gr_ripley(adata=self.getAdata(), cluster_key=cluster_key, mode=mode)

    # graph Ripley's statistics
    def pl_ripley(self, cluster_key="leiden", mode="L"):
        return tools.pl_ripley(adata=self.getAdata(), cluster_key=cluster_key, mode=mode)

    def gr_spatialAutocorr(self, adata=None, mode="moran", nPerms=100, nJobs=1):
        if adata is None:
            tools.gr_spatialAutocorr(adata=self.getAdata(), mode=mode, nPerms=nPerms, nJobs=nJobs)

        else:
            tools.gr_spatialAutocorr(adata=adata, mode=mode, nPerms=nPerms, nJobs=nJobs)

    # Reference Gene Functions (metaGene)
    def assignReferenceCells(self):
        # Assign marker gene metadata using reference dataset

        # get the reference panel
        gene_panel = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41421-021-00266-1/MediaObjects/41421_2021_266_MOESM1_ESM.xlsx"
        df_ref_panel_ini = pd.read_excel(gene_panel, index_col=0)
        self.setDfRefPanel(df_ref_panel_ini.iloc[1:, :1])
        self.getDfRefPanel().index.name = None
        self.getDfRefPanel().columns = ["Function"]

        # assign the values
        marker_genes = self.getDfRefPanel()[
            self.getDfRefPanel()["Function"].str.contains("marker")
        ].index.tolist()

        self.setMetaGene(deepcopy(self.getAdata().var))
        self.setCommonMarkerGenes(list(set(self.getMetaGene().index.tolist()).intersection(marker_genes)))
        self.getMetaGene().loc[self.getCommonMarkerGenes(), "Markers"] = self.getDfRefPanel().loc[
            self.getCommonMarkerGenes(), "Function"
        ]
        self.getMetaGene()["Markers"] = self.getMetaGene()["Markers"].apply(
            lambda x: "N.A." if "marker" not in str(x) else x
        )
        print(self.getMetaGene()["Markers"].value_counts())

    def calculateLeidenExpressionSignatures(self):
        ser_counts = self.getAdata().obs["leiden"].value_counts()
        ser_counts.name = "cell counts"
        self.setMetaLeiden(pd.DataFrame(ser_counts))

        cat_name = "leiden"
        self.setSigLeiden(pd.DataFrame(
            columns=self.getAdata().var_names, index=self.getAdata().obs[cat_name].cat.categories
        ))
        for clust in self.getAdata().obs[cat_name].cat.categories:
            self.getSigLeiden().loc[clust] = self.getAdata()[self.getAdata().obs[cat_name].isin([clust]), :].X.mean(0)
        self.setSigLeiden(self.getSigLeiden().transpose())
        self.setLeidenClusters(["Leiden-" + str(x) for x in self.getSigLeiden().columns.tolist()])
        self.getSigLeiden().columns = self.getLeidenClusters()
        self.getMetaLeiden().index = self.getSigLeiden().columns.tolist()
        self.getMetaLeiden()["leiden"] = pd.Series(
            self.getMetaLeiden().index.tolist(), index=self.getMetaLeiden().index.tolist()
        )

    def assignCellTypesOnExpression(self):
        meta_gene = pd.DataFrame(index=self.getSigLeiden().index.tolist())
        meta_gene["info"] = pd.Series("", index=meta_gene.index.tolist())
        meta_gene["Markers"] = pd.Series("N.A.", index=self.getSigLeiden().index.tolist())
        meta_gene.loc[self.getCommonMarkerGenes(), "Markers"] = self.getDfRefPanel().loc[
            self.getCommonMarkerGenes(), "Function"
        ]

        print("GOT HERE 3")

        self.getMetaLeiden()["Cell_Type"] = pd.Series("N.A.", index=self.getMetaLeiden().index.tolist())
        num_top_genes = 30
        for inst_cluster in self.getSigLeiden().columns.tolist():
            top_genes = (
                self.getSigLeiden()[inst_cluster]
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
            self.getMetaLeiden().loc[inst_cluster, "Cell_Type"] = max_cat

        # rename clusters
        self.getMetaLeiden()["name"] = self.getMetaLeiden().apply(
            lambda x: x["Cell_Type"] + "_" + x["leiden"], axis=1
        )
        leiden_names = self.getMetaLeiden()["name"].values.tolist()
        self.getMetaLeiden().index = leiden_names

        # transfer cell type labels to single cells
        leiden_to_cell_type = deepcopy(self.getMetaLeiden())
        leiden_to_cell_type.set_index("leiden", inplace=True)
        leiden_to_cell_type.index.name = None

        self.getAdata().obs["Cell_Type"] = self.getAdata().obs["leiden"].apply(
            lambda x: leiden_to_cell_type.loc["Leiden-" + str(x), "Cell_Type"]
        )
        self.getAdata().obs["Cluster"] = self.getAdata().obs["leiden"].apply(
            lambda x: leiden_to_cell_type.loc["Leiden-" + str(x), "name"]
        )

    def calcSerCloseness(self):
        # closeness centrality - measure of how close the group is to other nodes.
        self.setSerCloseness(self.getDfCentral()["closeness_centrality"].sort_values(ascending=False))

    def calcSerDegree(self):
        # The degree centrality for a node v is the fraction of nodes it is connected to.
        self.setSerDegree(self.getDfCentral()["degree_centrality"].sort_values(ascending=False))

    def calcSerCluster(self):
        # clustering coefficient - measure of the degree to which nodes cluster together.
        self.setSerCluster(self.getDfCentral()["average_clustering"].sort_values(ascending=False))

    def highClosenessScore(self):
        inst_clusters = self.getSerCloseness().index.tolist()[:5]
        self.spatialScatter(groups=inst_clusters, graphs="Cluster")

    def lowClosenessScore(self):
        inst_clusters = self.getSerCloseness().index.tolist()[-5:]
        self.spatialScatter(graphs="Cluster", groups=inst_clusters)

    def highDegreeCentrality(self):
        inst_clusters = self.getSerDegree().index.tolist()[:5]
        self.spatialScatter(graphs="Cluster", groups=inst_clusters)

    def lowDegreeCentrality(self):
        inst_clusters = self.getSerDegree().index.tolist()[-5:]
        self.spatialScatter(graphs="Cluster", groups=inst_clusters)

    def highClusteringCoefficient(self):
        inst_clusters = self.getSerCluster().index.tolist()[:5]
        self.spatialScatter(graphs="Cluster", groups=inst_clusters)

    def lowClusteringCoefficient(self):
        inst_clusters = self.getSerCluster().index.tolist()[-5:]
        self.spatialScatter(graphs="Cluster", groups=inst_clusters)

    # endregion

    # searches the color directory for the given file. If it is found, it returns the colors and sets the leiden colors.
    # If it isn't found, an error will be thrown.
    def setLeidenColors(self, color_file):
        leidenColors = tools.getColors(color_file=color_file)
        tools.setLeidenColors(adata=self.getAdata(), colors=leidenColors)
        return leidenColors

    def showPlots(self):
        tools.showPlots()

    ####################################################################################################################
    ####################################################################################################################
    # region print functions
    def print(self):
        print("\nOBJECT INFO")
        print("object: {}".format(self))
        print("object ID: {}".format(self.getID()))

        print("\nDATA")
        print(" data path: {}".format(self.getDataPath()))
        print(" adata:\n {}".format(self.getAdata()))

        print("")
    # endregion


if __name__ == "__main__":
    path = 'F:/sunlabmerfishdata/QSFL01222023/'
    esqAnalysis = Analysis(data_path=path)
    esqAnalysis.print()

    esqAnalysis.availableGraphs()

    print(esqAnalysis.qcMetrics(percentTop=(10, 50, 100)))

    esqAnalysis.filterCells(minCounts=50)
    esqAnalysis.filterGenes(minCells=10)

    esqAnalysis.setLeidenColors(color_file="leiden_generated_random_2.txt")
    esqAnalysis.leiden(resolution=1)

    esqAnalysis.spatialScatter(graphs=["leiden"])

    time.sleep(10000)

    esqAnalysis.availableGraphs()

    print("normalize")
    esqAnalysis.normalizeTotal()
    esqAnalysis.log1p()
    esqAnalysis.scale(maxValue=10)

    print("tl pca")
    esqAnalysis.tl_pca()
    print("neighbors")
    esqAnalysis.pp_neighbors()
    print("tl umap")
    esqAnalysis.tl_umap()
    print("leiden")
    esqAnalysis.leiden(resolution=1.5)

    print("pl umap")

    esqAnalysis.availableGraphs()

    esqAnalysis.pl_umap()

    esqAnalysis.spatialScatter(graphs=['leiden'])

    # neighborhood enrichment
    esqAnalysis.gr_spatialNeighbors()
    esqAnalysis.gr_nhoodEnrichment()
    esqAnalysis.pl_nhoodEnrichment()
