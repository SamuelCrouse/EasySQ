# Sam Crouse
# scrouse2@uwyo.edu

# Class which implements the tools built in EasySQTools.py
# Allows creation of a squidpy analysis object with built in tools

# python imports

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
    def spatialScatter(self, graphs, colorInit=False, show=False, colors=None, libraryID=None, wspace=0.4, size=None,
                       shape=None):
        return tools.spatialScatter(adata=self.getAdata(), graphs=graphs, show=show,
                                    colors=colors, libraryID=libraryID, wspace=wspace, size=size, shape=None)

    # runs the adata setup and basic analysis tools for the user
    def adataSetup(self):
        return tools.adataSetup(self.getAdata())

    # adata functions
    # calculate qc metrics
    def qcMetrics(self, percentTop=(50, 100), inplace=True):
        return tools.qcMetrics(self.getAdata(), percentTop, inplace=inplace)

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

    def neighbors(self):
        return tools.neighbors(self.getAdata())

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

    def plotTranscripts(self, show=True, figsize=(15, 4)):
        return tools.plotTranscripts(self.getAdata(), show=show, figsize=figsize)

    # uses the sc.pl.umap function plot a umap. Use colorInit to generate a custom color palette by default
    def pl_umap(self, graphs=["leiden"], show=False, colorInit=False, size=None, wspace=0.4):
        return tools.pl_umap(self.getAdata(), graphs=graphs, show=show, size=size, wspace=wspace)

    def assignCellTypes(self):
        return tools.assignCellTypes(self.getAdata())

    def spatialNeighbors(self, coordType="generic", spatialKey="spatial", delaunay=False):
        return tools.spatialNeighbors(adata=self.getAdata(), coordType=coordType, spatialKey=spatialKey,
                                      delaunay=delaunay)

    # calculate nhoodEnrichment
    def gr_nhoodEnrichment(self, clusterKey="leiden"):
        return tools.gr_nhoodEnrichment(adata=self.getAdata(), clusterKey=clusterKey)

    # plot nhoodEnrichment data
    def pl_nhoodEnrichment(self, show=True, clusterKey="leiden", method="average", cmap="inferno", vmin=-50, vmax=100,
                           figsize=(5, 5)):
        tools.pl_nhoodEnrichment(adata=self.getAdata(), show=show, clusterKey=clusterKey, method=method,
                                 cmap=cmap, vmin=vmin, vmax=vmax, figsize=figsize)

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
    esqAnalysis.neighbors()
    print("tl umap")
    esqAnalysis.tl_umap()
    print("leiden")
    esqAnalysis.leiden(resolution=1.5)

    print("pl umap")

    esqAnalysis.availableGraphs()

    esqAnalysis.pl_umap()

    esqAnalysis.spatialScatter(graphs=['leiden'])

    print(esqAnalysis.assignCellTypes())

    # neighborhood enrichment
    esqAnalysis.spatialNeighbors()
    esqAnalysis.gr_nhoodEnrichment()
    esqAnalysis.pl_nhoodEnrichment()
