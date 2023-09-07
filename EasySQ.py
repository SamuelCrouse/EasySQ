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
        Documentation

        Analysis class object documentation in EasySQ.py
        =============================
        ----

        * Acts as a Squidpy interface. Many of the functions used in vizgen analysis of spatial data are represented
          here. For example, reading in data is now simplified. You only need to provide a path, have your meta files
          contain these names 'cell_metadata' and 'cell_by_gene' and EasySQ will automatically import your data into an
          AnnData object for use with Squidpy.
        * AnnData is kept track of in this class. Any Squidpy (supported) Squidpy function that uses AnnData can be
          called using this class object. Operations will then be done to the stored AnnData.
        * Follow these demos to see how to work with EasySQ.
        Todo, make demos and link them here.

        ----

        Examples:
        --------
        1. Creating and working with an Analysis class object:
            * import EasySQ as esq
            * esqAn = esq.Analysis(data_path=pathToData)
        ..
        2. Print example:
            * esqAn.print()
            * This calls the class function print. Which prints the AnnData object among other information.
    """

    """
    # Some rules for data in this class.
    # 1. If a variable isn't set, it should be None.
    # 2. All values initialized in constructor. If nothing is passed, set to None.
    # 3. Object.variable should not be used. Variables should only be accessed through getter and setter functions.
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

        self.topAutoCorr = None
        self.botAutoCorr = None

    # endregion constructor

    ####################################################################################################################
    ####################################################################################################################
    # region setter functions
    """
        Documentation

        EasySQ Analysis Class Setter Functions
        =============================
        
        These are used mostly internally to set class variables.
    """

    # set the ID of this object
    def setID(self, ID):
        """
            Documentation

            setID() EasySQ Class function.
            =============================
            ----

            * Used to set the ID class variable.
            * Which is used to keep track of this object.

            ----

            Parameters:
            ----------
            * ID: The value for the ID of this object.
            * type: any

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setID():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setID("myID")
        """

        self.ID = str(ID)
        return None

    def setDataPath(self, data_path):
        """
            Documentation

            setDataPath() EasySQ Class function.
            =============================
            ----

            * Used to set the data_path class variable.
            * Which is used to read user provided meta data.

            ----

            Parameters:
            ----------
            * data_path: The file path to the directory containing the data.
            * type: string

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setDataPath():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setDataPath("path/to/my/dir/")
        """

        self.data_path = data_path
        return None

    def setAdata(self, adata):
        """
            Documentation

            setAdata() EasySQ Class function.
            =============================
            ----

            * Used to set the adata class variable.
            * This is what all analyses are performed on.
            * Note, overwriting adata with this function could interfere with saved Analysis class data.
            * If you need to supply another adata object, I would recommend creating a new Analysis class object.

            ----

            Parameters:
            ----------
            * adata: The AnnData object you wish to set to the class object.
            * type: AnnData object.

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setAdata():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setAdata(myAnnDataObject)
            ..
            2. AnnData Note:
                * An AnnData object is created on class object instantiation, so calling this manually is not recommended.
        """

        self.adata = adata
        return None

    def setMetaGene(self, value):
        """
            Documentation

            setMetaGene() EasySQ Class function.
            =============================
            ----

            * Used to set the metaGene class variable.
            * Used in assigning reference cells.

            ----

            Parameters:
            ----------
            * value: The value you wish to set self.metaGene to.
            * type: any

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setMetaGene():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setMetaGene(value)
        """

        self.metaGene = value
        return None

    def setCommonMarkerGenes(self, value):
        """
            Documentation

            setCommonMarkerGenes() EasySQ Class function.
            =============================
            ----

            * Used to set the commonMarkerGenes class variable.
            * Used in assigning reference cells.

            ----

            Parameters:
            ----------
            * value: The value you wish to set self.commonMarkerGenes to.
            * type: any

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setCommonMarkerGenes():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setCommonMarkerGenes(value)
        """

        self.commonMarkerGenes = value
        return None

    def setDfRefPanel(self, value):
        """
            Documentation

            setDfRefPanel() EasySQ Class function.
            =============================
            ----

            * Used to set the dfRefPanel class variable.
            * Used in assigning reference cells.

            ----

            Parameters:
            ----------
            * value: The value you wish to set self.dfRefPanel to.
            * type: any

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setDfRefPanel():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setDfRefPanel(value)
        """

        self.dfRefPanel = value
        return None

    def setSigLeiden(self, value):
        """
            Documentation

            setSigLeiden() EasySQ Class function.
            =============================
            ----

            * Used to set the sigLeiden class variable.
            * Used in calculating leiden expression signatures.

            ----

            Parameters:
            ----------
            * value: The value you wish to set self.sigLeiden to.
            * type: any

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setSigLeiden():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setSigLeiden(value)
        """

        self.sigLeiden = value
        return None

    def setMetaLeiden(self, value):
        """
            Documentation

            setMetaLeiden() EasySQ Class function.
            =============================
            ----

            * Used to set the metaLeiden class variable.
            * Used in calculating leiden expression signatures.

            ----

            Parameters:
            ----------
            * value: The value you wish to set self.metaLeiden to.
            * type: any

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setMetaLeiden():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setMetaLeiden(value)
        """

        self.metaLeiden = value
        return None

    def setLeidenClusters(self, value):
        """
            Documentation

            setLeidenClusters() EasySQ Class function.
            =============================
            ----

            * Used to set the leidenClusters class variable.
            * Used in calculating leiden expression signatures.

            ----

            Parameters:
            ----------
            * value: The value you wish to set self.leidenClusters to.
            * type: any

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setLeidenClusters():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setLeidenClusters(value)
        """

        self.leidenClusters = value
        return None

    def setDfCentral(self, value):
        """
            Documentation

            setDfCentral() EasySQ Class function.
            =============================
            ----

            * Used to set the dfCentral class variable.
            * Used in demos to store various data.

            ----

            Parameters:
            ----------
            * value: The value you wish to set self.dfCentral to.
            * type: any

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setDfCentral():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setDfCentral(value)
        """

        self.dfCentral = value
        return None

    def setSerCloseness(self, value):
        """
            Documentation

            setSerCloseness() EasySQ Class function.
            =============================
            ----

            * Used to set the serCloseness class variable.
            * Used in calcSerCloseness() in tandem with getDfCentral().

            ----

            Parameters:
            ----------
            * value: The value you wish to set self.serCloseness to.
            * type: any

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setSerCloseness():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setSerCloseness(value)
        """

        self.serCloseness = value
        return None

    def setSerDegree(self, value):
        """
            Documentation

            setSerDegree() EasySQ Class function.
            =============================
            ----

            * Used to set the serDegree class variable.
            * Used in calcSerDegree() in tandem with getDfCentral().

            ----

            Parameters:
            ----------
            * value: The value you wish to set self.serDegree to.
            * type: any

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setSerDegree():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setSerDegree(value)
        """

        self.serDegree = value
        return None

    def setSerCluster(self, value):
        """
            Documentation

            setSerCluster() EasySQ Class function.
            =============================
            ----

            * Used to set the serCluster class variable.
            * Used in calcSerCluster() in tandem with getDfCentral().

            ----

            Parameters:
            ----------
            * value: The value you wish to set self.serCluster to.
            * type: any

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setSerCluster():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setSerCluster(value)
        """

        self.serCluster = value
        return None

    def setTopAutoCorr(self, value):
        """
            Documentation

            setTopAutoCorr() EasySQ Class function.
            =============================
            ----

            * Used to set the topAutoCorr class variable.
            * Used in the calculation of Moran's I Score.
            * calcMoransIScore()

            ----

            Parameters:
            ----------
            * value: The value you wish to set self.topAutoCorr to.
            * type: any

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setTopAutoCorr():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setTopAutoCorr(value)
        """

        self.topAutoCorr = value
        return None

    def setBotAutoCorr(self, value):
        """
            Documentation

            setBotAutoCorr() EasySQ Class function.
            =============================
            ----

            * Used to set the botAutoCorr class variable.
            * Used in the calculation of Moran's I Score.
            * calcMoransIScore()

            ----

            Parameters:
            ----------
            * value: The value you wish to set self.botAutoCorr to.
            * type: any

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling setBotAutoCorr():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * esqAn.setBotAutoCorr(value)
        """

        self.botAutoCorr = value
        return None

    # endregion

    ####################################################################################################################
    ####################################################################################################################
    # region getter functions
    """
        Documentation

        EasySQ Analysis Class Getter Functions
        =============================

        These are used mostly internally to get class variables.
    """

    def getID(self):
        """
            Documentation

            getID() EasySQ Class function.
            =============================
            ----

            * Used to get the ID class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.ID
            * type: default is string, but can be user defined.

            Examples:
            --------
            1. Calling getID():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * ID = esqAn.getID()
        """

        return self.ID

    def getDataPath(self):
        """
            Documentation

            getDataPath() EasySQ Class function.
            =============================
            ----

            * Used to get the data_path class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.data_path
            * type: default is string, but can be user defined.

            Examples:
            --------
            1. Calling getDataPath():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * data_path = esqAn.getDataPath()
        """

        return self.data_path

    def getAdata(self):
        """
            Documentation

            getAdata() EasySQ Class function.
            =============================
            ----

            * Used to get the adata class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.adata
            * type: default is AnnData object, but can be user defined. (Don't recommend.)

            Examples:
            --------
            1. Calling getAdata():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * adata = esqAn.getAdata()
        """

        return self.adata

    def getMetaGene(self):
        """
            Documentation

            getMetaGene() EasySQ Class function.
            =============================
            ----

            * Used to get the metaGene class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.metaGene
            * type: Automatically set, can be user defined.

            Examples:
            --------
            1. Calling getMetaGene():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * metaGene = esqAn.getMetaGene()
        """

        return self.metaGene

    def getCommonMarkerGenes(self):
        """
            Documentation

            getCommonMarkerGenes() EasySQ Class function.
            =============================
            ----

            * Used to get the commonMarkerGenes class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.commonMarkerGenes
            * type: Automatically set, can be user defined.

            Examples:
            --------
            1. Calling getCommonMarkerGenes():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * commonMarkerGenes = esqAn.getCommonMarkerGenes()
        """

        return self.commonMarkerGenes

    def getDfRefPanel(self):
        """
            Documentation

            getDfRefPanel() EasySQ Class function.
            =============================
            ----

            * Used to get the dfRefPanel class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.dfRefPanel
            * type: Automatically set, can be user defined.

            Examples:
            --------
            1. Calling getDfRefPanel():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * dfRefPanel = esqAn.getDfRefPanel()
        """

        return self.dfRefPanel

    def getSigLeiden(self):
        """
            Documentation

            getSigLeiden() EasySQ Class function.
            =============================
            ----

            * Used to get the sigLeiden class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.sigLeiden
            * type: Automatically set, can be user defined.

            Examples:
            --------
            1. Calling getSigLeiden():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * sigLeiden = esqAn.getSigLeiden()
        """

        return self.sigLeiden

    def getMetaLeiden(self):
        """
            Documentation

            getMetaLeiden() EasySQ Class function.
            =============================
            ----

            * Used to get the metaLeiden class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.metaLeiden
            * type: Automatically set, can be user defined.

            Examples:
            --------
            1. Calling getMetaLeiden():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * metaLeiden = esqAn.getMetaLeiden()
        """

        return self.metaLeiden

    def getLeidenClusters(self):
        """
            Documentation

            getLeidenClusters() EasySQ Class function.
            =============================
            ----

            * Used to get the leidenClusters class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.leidenClusters
            * type: Automatically set, can be user defined.

            Examples:
            --------
            1. Calling getLeidenClusters():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * leidenClusters = esqAn.getLeidenClusters()
        """

        return self.leidenClusters

    def getDfCentral(self):
        """
            Documentation

            getDfCentral() EasySQ Class function.
            =============================
            ----

            * Used to get the dfCentral class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.dfCentral
            * type: Automatically set, can be user defined.

            Examples:
            --------
            1. Calling getDfCentral():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * dfCentral = esqAn.getDfCentral()
        """

        return self.dfCentral

    def getSerCloseness(self):
        """
            Documentation

            getSerCloseness() EasySQ Class function.
            =============================
            ----

            * Used to get the serCloseness class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.serCloseness
            * type: Automatically set, can be user defined.

            Examples:
            --------
            1. Calling getSerCloseness():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * serCloseness = esqAn.getSerCloseness()
        """

        return self.serCloseness

    def getSerDegree(self):
        """
            Documentation

            getSerDegree() EasySQ Class function.
            =============================
            ----

            * Used to get the serDegree class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.serDegree
            * type: Automatically set, can be user defined.

            Examples:
            --------
            1. Calling getSerDegree():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * serDegree = esqAn.getSerDegree()
        """

        return self.serDegree

    def getSerCluster(self):
        """
            Documentation

            getSerCluster() EasySQ Class function.
            =============================
            ----

            * Used to get the serCluster class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.serCluster
            * type: Automatically set, can be user defined.

            Examples:
            --------
            1. Calling getSerCluster():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * serCluster = esqAn.getSerCluster()
        """

        return self.serCluster

    def getTopAutoCorr(self):
        """
            Documentation

            getTopAutoCorr() EasySQ Class function.
            =============================
            ----

            * Used to get the topAutoCorr class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.topAutoCorr
            * type: Automatically set, can be user defined.

            Examples:
            --------
            1. Calling getTopAutoCorr():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * topAutoCorr = esqAn.getTopAutoCorr()
        """

        return self.topAutoCorr

    def getBotAutoCorr(self):
        """
            Documentation

            getBotAutoCorr() EasySQ Class function.
            =============================
            ----

            * Used to get the botAutoCorr class variable.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * self.botAutoCorr
            * type: Automatically set, can be user defined.

            Examples:
            --------
            1. Calling getBotAutoCorr():
                * Using the EasySQ Analysis class object you created, called "esqAn":
                * botAutoCorr = esqAn.getBotAutoCorr()
        """

        return self.botAutoCorr

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

            * qcMetrics() class function for EasySQ class.
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
                * perUn = esqAn.qcMetrics()
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

    def gr_spatialAutocorr(self, adata=None, mode="moran", nPerms=None, nJobs=None):
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

    def calcMoransIScore(self, numView=12):
        """
            Documentation

            calcMoransIScore() class function docs
            =============================
            ----

            * calcMoransIScore() class function for EasySQ class.
            * Calculates the top and bottom autocorrelation using gr_spatialAutocorr()

            ----

            Parameters:
            ----------
            * numView=12: Number of top and bottom autocorr we calculate for

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling calcMoransIScore():
                * Setup an EasySQ Analysis class object called esqAn, then call with:
                * esqAn.calcMoransIScore(numView=12) or
                * esqAn.calcMoransIScore(numView=6) or
                * esqAn.calcMoransIScore(numView=n)
            ..
            2. Displaying the graphs:
                * After calling esqAn.calcMoransIScore() you can call esqAn.plotMoransIScore()
                * Calling esqAn.showPlots() will display these graphs.
        """

        self.gr_spatialAutocorr(mode='moran')
        num_view = numView
        self.setTopAutoCorr(
            self.getAdata().uns["moranI"]["I"].sort_values(ascending=False).head(num_view).index.tolist())

        self.setBotAutoCorr(
            self.getAdata().uns["moranI"]["I"].sort_values(ascending=True).head(num_view).index.tolist())

        return None

    def plotMoransIScore(self, size=None, figsize=(3, 3)):
        """
            Documentation

            plotMoransIScore() class function docs
            =============================
            ----

            * plotMoransIScore() class function for EasySQ class.
            * Plots the graphs calculated by calcMoransIScore() class function

            ----

            Parameters:
            ----------
            * figsize=(3, 3): the size of the window
            * size=None: Size of the scatter point/shape.

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling plotMoransIScore():
                * Setup an EasySQ Analysis class object called esqAn.
                * Call calcMoransIScore to set up the autocorrelation values.
                * Then, call with:
                * esqAn.plotMoransIScore(figsize=(3, 3)) or
                * esqAn.plotMoransIScore(figsize=(5, 5)) or
                * esqAn.plotMoransIScore(figsize=(n, n))
            ..
            2. Displaying the graphs:
                * After calling esqAn.calcMoransIScore() you can call esqAn.plotMoransIScore()
                * Calling esqAn.showPlots() will display these graphs.
        """

        # print("Genes with high autocorrelation")
        self.spatialScatter(graphs=self.getTopAutoCorr(), cmap="Reds", size=size, figsize=figsize)

        # print("Genes with low autocorrelation")
        self.spatialScatter(graphs=self.getBotAutoCorr(), cmap="Reds", size=size, figsize=figsize)

        return None
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
