# Sam Crouse
# scrouse2@uwyo.edu

# Class which implements the tools built in EasySQTools.py
# Allows creation of a squidpy analysis object with built in tools

# python imports
from copy import deepcopy
import pandas as pd
from scipy.cluster import hierarchy as sch
import anndata as ad

# my imports
import EasySQTools as tools


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
        * The demos can be found at (https://github.com/SamuelCrouse/EasySQ)
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
    def __init__(self, data_path=""):
        self.ID = tools.createCode(9)

        # passed value initialization
        # this lets us create an empty Analysis object by passing ""
        if data_path != "":
            self.data_path = data_path
            try:
                self.adata = self.readVizgen(self.data_path)

            except FileNotFoundError:
                self.adata = ad.read_h5ad(self.data_path)

        else:
            self.data_path = None
            self.adata = None

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
    def readVizgen(self, data_path, transform_path=None, transform_file='micron_to_mosaic_pixel_transform.csv'):
        """
            Documentation

            readVizgen() class function docs
            =============================
            ----

            * readVizgen() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools readVizgen() function.
            * Read in and create adata based on metadata files.

            ----

            Parameters:
            ----------
            * data_path="path/to/my/data/": The path to the directory containing your metadata.
            * transform_path="path/to/my/dir/": Path to the root directory which must hold an 'images' directory. This dir should contain a transform_file called, by default 'micron_to_mosaic_pixel_transform.csv.' The path format looks like this ...users/pythonproject/ and your project structure will be images/micron_to_mosaic_pixel_transform.csv. You can mostly ignore these two parameters.
            * transform_file='micron_to_mosaic_pixel_transform.csv': Can be ignored unless you want to provide your own transform file.
            * This must be the name of the transform file stored in the 'images' directory.

            Returns:
            -------
            * AnnData object that contains your meta data.

            Examples:
            --------
            1. Calling readVizgen()
                * newAdata = esqAn.readVizgen()
                I don't recommend doing this, because it is bad design to use a current class object that manages its
                own data to create a new, unassociated dataset; especially when there are methods in place to do so.
                If you need to work with a new set of data, or even if you plan to use them in tandem, then you should
                create a new 'esq' object and pass it the directory containing your data.
            ..
            2. When is readVizgen() called? During the import statement automatically.:
                * path = os.getcwd() + '\\tutorial_data_1\\'
                * esqAn = esq.Analysis(data_path=path)
                * esqAn.print()
        """

        return tools.readVizgen(data_path=data_path, transform_path=transform_path, transform_file=transform_file)

    def availableGraphs(self, log=True):
        """
            Documentation

            availableGraphs() class function docs
            =============================
            ----

            * availableGraphs() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools availableGraphs() function.
            * Log the graphs available to plot, then return them.

            ----

            Parameters:
            ----------
            * log=True: Bool determining whether to print the available graphs to the console or not.

            Returns:
            -------
            * A list of strings which are the available graphs.

            Examples:
            --------
            1. Calling availableGraphs()
                * esqAn.availableGraphs()
                * Prints the graphs available for the stored adata.
        """

        return tools.availableGraphs(self.getAdata(), log=log)

    # plot a spatial scatter. use colorInit to generate a custom color palette by default
    def spatialScatter(self, graphs, adata=None, show=False, libraryID=None, wspace=None, size=None,
                       shape=None, groups=None, cmap=None, figsize=None):
        """
            Documentation

            spatialScatter() class function docs
            =============================
            ----

            * spatialScatter() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools spatialScatter() function.
            * Produces a spatial scatter plot given a graph type and adata.

            ----

            Parameters:
            ----------
            * graphs=["leiden", "Slc17a7"]: A list of strings containing the graphs to plot.
            * adata=AnnData Object or None: An AnnData Object that contains the data you wish to plot. If None is chosen, will use the currently tracked adata. Should only supply this if you are using a subsample.
            * show=False: Whether to show the graphs now, or wait for later.
            * libraryID="id1": A string containing the libraryID. Ignore this unless you need to set one.
            * wspace=0.4: Dictates the width of space between the panels.
            * size=4: Size of the scatter point/shape.
            * shape='circle': Shape of the scatter point. 'circle', 'square', 'hex'.
            * groups=["0", "1", "3"]: Select which values to plot.
            * cmap="Reds": Colormap for continuous annotations.
            * figsize=(12, 12): Size of the figure in inches.

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling spatialScatter()
                * esqAn.spatialScatter(graphs=["leiden"])
                * Calls spatial scatter on an existing analysis class object, and runs on the stored adata.
        """

        if adata is None:
            return tools.spatialScatter(adata=self.getAdata(), graphs=graphs, show=show, libraryID=libraryID,
                                        wspace=wspace, size=size, shape=shape, groups=groups, cmap=cmap,
                                        figsize=figsize)

        else:
            return tools.spatialScatter(adata=adata, graphs=graphs, show=show, libraryID=libraryID, wspace=wspace,
                                        size=size, shape=shape, groups=groups, cmap=cmap, figsize=figsize)

    # runs the adata setup and basic analysis tools for the user
    def adataSetup(self):
        """
            Documentation

            adataSetup() class function docs
            =============================
            ----

            * adataSetup() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools adataSetup() function.
            * Runs some basic analysis functions for a generic analysis.
            * These include qcMetrics, leiden, among others.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling adataSetup()
                * esqAn.adataSetup()
                * Runs the setup on the adata stored within esqAn
        """

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
            * Acts as a wrapper for EasySQTools qcMetrics() function.
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
        """
            Documentation

            layers() class function docs
            =============================
            ----

            * layers() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools layers() function.
            * Returns the layer named "counts" and sets it to the value of the adata variables.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling layers()
                * esqAn.layers()
                * Runs layers() on the stored adata.
        """

        return tools.layers(self.getAdata())

    def highlyVariableGenes(self, nTopGenes=4000):
        """
            Documentation

            highlyVariableGenes() class function docs
            =============================
            ----

            * highlyVariableGenes() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools highlyVariableGenes() function.
            * Returns the Scanpy sc.pp.highly_variable_genes() function value.

            ----

            Parameters:
            ----------
            * nTopGenes=4000: The number of highly variable genes to keep

            Returns:
            -------
            * Returns the value of sc.pp.highly_variable_genes(). Updates .var or returns depending on "inplace"

            Examples:
            --------
            1. Calling highlyVariableGenes()
                * esqAn.highlyVariableGenes(nTopGenes=3000)
                * Runs on the stored adata.
        """

        return tools.highlyVariableGenes(self.getAdata(), nTopGenes)

    def normalizeTotal(self, inplace=True):
        """
            Documentation

            normalizeTotal() class function docs
            =============================
            ----

            * normalizeTotal() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools normalizeTotal() function.
            * Normalizes the counts per cells.
            * Returns the Scanpy sc.pp.normalize_total() function value.

            ----

            Parameters:
            ----------
            * inPlace=True: Whether to perform in place or return the data calculated.

            Returns:
            -------
            * Returns the value of sc.pp.normalize_total(). Updates or returns depending on "inplace"

            Examples:
            --------
            1. Calling normalizeTotal()
                * esqAn.normalizeTotal()
                * Runs on the stored adata.
        """

        return tools.normalizeTotal(self.getAdata(), inplace=inplace)

    def log1p(self):
        """
            Documentation

            log1p() class function docs
            =============================
            ----

            * log1p() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools log1p() function.
            * Logarithmize the data matrix.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * Returns the value of sc.pp.log1p().

            Examples:
            --------
            1. Calling log1p()
                * esqAn.log1p()
                * Runs on the stored adata.
        """

        return tools.log1p(self.getAdata())

    def pp_pca(self):
        """
            Documentation

            pp_pca() class function docs
            =============================
            ----

            * pp_pca() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools pp_pca() function.
            * Runs principle component analysis on the stored adata.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * Returns the value of sc.pp.pca().

            Examples:
            --------
            1. Calling pp_pca()
                * esqAn.pp_pca()
                * Runs on the stored adata.
        """

        return tools.pp_pca(self.getAdata())

    def tl_pca(self, svdSolver="arpack"):
        """
            Documentation

            tl_pca() class function docs
            =============================
            ----

            * tl_pca() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools tl_pca() function.
            * Runs principle component analysis on the stored adata.
            * The difference between tl_pca and pl_pca is the function that it wraps.
            * Some demos use pl and tl in different scenarios.

            ----

            Parameters:
            ----------
            * svdSolver="arpack": The SVD solver to use.

            Returns:
            -------
            * Returns the value of sc.tl.pca().

            Examples:
            --------
            1. Calling tl_pca()
                * esqAn.tl_pca()
                * Runs on the stored adata.
        """

        return tools.tl_pca(self.getAdata(), svdSolver=svdSolver)

    def pp_neighbors(self, nNeighbors=15, nPcs=None):
        """
            Documentation

            pp_neighbors() class function docs
            =============================
            ----

            * pp_neighbors() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools pp_neighbors() function.
            * Computes a neighborhood graph of observations.

            ----

            Parameters:
            ----------
            * nNeighbors=10: The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation.
            * nPcs=20: Use this many PCs. If n_pcs==0 use .X if use_rep is None.

            Returns:
            -------
            * Returns the value of sc.pp.neighbors().

            Examples:
            --------
            1. Calling pp_neighbors()
                * esqAn.pp_neighbors()
                * Runs on the stored adata.
        """

        return tools.pp_neighbors(self.getAdata(), nNeighbors=nNeighbors, nPcs=nPcs)

    def tl_umap(self):
        """
            Documentation

            tl_umap() class function docs
            =============================
            ----

            * tl_umap() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools tl_umap() function.
            * Computes a neighborhood graph of observations.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * Returns the value of sc.tl.umap().

            Examples:
            --------
            1. Calling tl_umap()
                * esqAn.tl_umap()
                * Runs on the stored adata.
        """

        return tools.tl_umap(self.getAdata())

    def leiden(self, resolution=1.0):
        """
            Documentation

            leiden() class function docs
            =============================
            ----

            * leiden() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools leiden() function.
            * Clusters cells into cell groups using the Leiden algorithm.

            ----

            Parameters:
            ----------
            * resolution=1.0: A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters.

            Returns:
            -------
            * Returns the value of sc.tl.leiden().

            Examples:
            --------
            1. Calling leiden()
                * esqAn.leiden()
                * Runs on the stored adata.
        """

        return tools.leiden(self.getAdata(), resolution=resolution)

    def filterCells(self, minCounts=50):
        """
            Documentation

            filterCells() class function docs
            =============================
            ----

            * filterCells() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools filterCells() function.
            * Filter cell outliers based on counts and numbers of genes expressed.

            ----

            Parameters:
            ----------
            * minCounts=50: The minimum number of genes expressed for a cell to pass filtering.

            Returns:
            -------
            * Returns the value of sc.pp.filter_cells().

            Examples:
            --------
            1. Calling filterCells()
                * esqAn.filterCells()
                * Runs on the stored adata.
        """

        return tools.filterCells(self.getAdata(), minCounts=minCounts)

    def filterGenes(self, minCells=10):
        """
            Documentation

            filterGenes() class function docs
            =============================
            ----

            * filterGenes() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools filterGenes() function.
            * Filter genes based on the number of cells or counts.

            ----

            Parameters:
            ----------
            * minCells=10: The minimum number of cells expressed required for a gene to pass filtering.

            Returns:
            -------
            * Returns the value of sc.pp.filter_genes().

            Examples:
            --------
            1. Calling filterGenes()
                * esqAn.filterGenes()
                * Runs on the stored adata.
        """

        return tools.filterGenes(self.getAdata(), minCells=minCells)

    # can be used to scale gene expression. IE Clip values that exceed 10 ('max value') standard deviations
    def scale(self, maxValue=10):
        """
            Documentation

            scale() class function docs
            =============================
            ----

            * scale() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools scale() function.
            * Scale data to unit variance and zero mean.
            * Can be used to scale gene expression. IE Clip values that exceed 10 ('max value') standard deviations.

            ----

            Parameters:
            ----------
            * maxValue=10: Clip (truncate) to this value after scaling. If None, do not clip.

            Returns:
            -------
            * Returns the value of sc.pp.scale().

            Examples:
            --------
            1. Calling scale()
                * esqAn.scale()
                * Runs on the stored adata.
        """

        return tools.scale(self.getAdata(), maxValue=maxValue)

    def plotTranscripts(self, show=False, figsize=(15, 4)):
        """
            Documentation

            plotTranscripts() class function docs
            =============================
            ----

            * plotTranscripts() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools plotTranscripts() function.
            * Plots calculated transcripts.

            ----

            Parameters:
            ----------
            * show=False: Whether to display graphs on function call, or later by calling esqAn.showPlots()
            * figsize=(20, 6): The dimensions of the plot (x, y).

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling plotTranscripts()
                * esqAn.plotTranscripts()
                * Runs on the stored adata.
        """

        return tools.plotTranscripts(self.getAdata(), show=show, figsize=figsize)

    # uses the sc.pl.umap function plot a umap. Use colorInit to generate a custom color palette by default
    def pl_umap(self, graphs=["leiden"], show=False, size=None, wspace=0.4):
        """
            Documentation

            pl_umap() class function docs
            =============================
            ----

            * pl_umap() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools pl_umap() function.
            * Similar to spatialScatter(). Takes graphs and then plots umap based on those graphs.

            ----

            Parameters:
            ----------
            * graphs=["leiden", "Slc17a7"]: A list of strings containing the graphs to plot.
            * show=False: Whether to show the graphs now, or wait for an esqAn.showPlots() call.
            * size=4: Size of the scatter point/shape.
            * wspace=0.4: Dictates the width of space between the panels.

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling pl_umap()
                * esqAn.pl_umap()
                * Runs on the stored adata.
        """

        return tools.pl_umap(self.getAdata(), graphs=graphs, show=show, size=size, wspace=wspace)

    def gr_spatialNeighbors(self, adata=None, coordType="generic", spatialKey="spatial", delaunay=False):
        """
            Documentation

            gr_spatialNeighbors() class function docs
            =============================
            ----

            * gr_spatialNeighbors() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools gr_spatialNeighbors() function.
            * Creates a graph using spatial coordinates.

            ----

            Parameters:
            ----------
            * adata=None: Run on a supplied adata or on the self adata. Usually used to run on a subsample. Recommend leaving this at None unless you need to run on a subsample.
            * coordType="generic": Type of coordinate system. 'grid' 'generic'
            * spatialKey="spatial": Key in anndta.AnnData.obsm where spatial coordinates are stored. You probably set this earlier.
            * delaunay=False: Whether to compute the graph from Delaunay triangulation. Only used when coord_type='generic'.

            Returns:
            -------
            * The value of sq.gr.spatial_neighbors().

            Examples:
            --------
            1. Calling gr_spatialNeighbors()
                * esqAn.gr_spatialNeighbors()
                * Runs on the stored adata.
            ..
            2. Calling gr_spatialNeighbors() with a subsample.
                * esqAn.gr_spatialNeighbors(adata=adataSubsample)
                * Runs on the provided subsample.
        """

        if adata is None:
            return tools.gr_spatialNeighbors(adata=self.getAdata(), coordType=coordType, spatialKey=spatialKey,
                                             delaunay=delaunay)

        else:
            return tools.gr_spatialNeighbors(adata=adata, coordType=coordType, spatialKey=spatialKey, delaunay=delaunay)

    def gr_nhoodEnrichment(self, clusterKey="leiden"):
        """
            Documentation

            gr_nhoodEnrichment() class function docs
            =============================
            ----

            * gr_nhoodEnrichment() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools gr_nhoodEnrichment() function.
            * Compute neighborhood enrichment by permutation test.

            ----

            Parameters:
            ----------
            * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.

            Returns:
            -------
            * The value of sq.gr.nhood_enrichment().

            Examples:
            --------
            1. Calling gr_nhoodEnrichment()
                * esqAn.gr_nhoodEnrichment()
                * Runs on the stored adata using default "leiden".
        """

        return tools.gr_nhoodEnrichment(adata=self.getAdata(), clusterKey=clusterKey)

    def pl_nhoodEnrichment(self, show=False, clusterKey="leiden", method="average", cmap="inferno", vmin=None,
                           vmax=None, figsize=(5, 5)):
        """
            Documentation

            pl_nhoodEnrichment() class function docs
            =============================
            ----

            * pl_nhoodEnrichment() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools pl_nhoodEnrichment() function.
            * Plot the neighborhood enrichment.

            ----

            Parameters:
            ----------X
            * show=False: Show now or later with esqAn.showPlots().
            * clusterKey="leiden": Key in anndata.AnnData.obs where graph data is stored. Computed from gr_nhoodEnrichment.
            * method="average": The linkage method to be used for dendrogram/clustering.
            * cmap="inferno": Continuous colormap to use.
            * vmin=-50: Lower bound for y axis on plot.
            * vmax=100: Upper bound for y axis on plot.
            * figsize=(5, 5): The dimensions of the plot (x, y).

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling pl_nhoodEnrichment()
                * esqAn.pl_nhoodEnrichment()
                * Runs on the stored adata.
            ..
            2. Calling pl_nhoodEnrichment() with show=True
                * esqAn.pl_nhoodEnrichment(show=True)
                * Runs on the stored adata and displays the plots immediately.
        """

        tools.pl_nhoodEnrichment(adata=self.getAdata(), show=show, clusterKey=clusterKey, method=method,
                                 cmap=cmap, vmin=vmin, vmax=vmax, figsize=figsize)

    def plotNHoodEnrichmentClusters(self, nClusters=[4], method="average", size=10, figsize=(5, 5)):
        """
            Documentation

            plotNHoodEnrichmentClusters() class function docs
            =============================
            ----

            * plotNHoodEnrichmentClusters() class function for EasySQ class.
            * Plot the neighborhood enrichment clusters. These were probably assigned earlier using these:
            * esqAn.assignReferenceCells()
            * esqAn.calculateLeidenExpressionSignatures()
            * esqAn.assignCellTypesOnExpression()

            ----

            Parameters:
            ----------
            * nClusters=[4]: Number of plots to display.
            * method="average": Method for calculating the distance between the newly formed clusters. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
            * size=10: Marker size of the plot.
            * figsize=(20, 10): Size of the plot.

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling plotNHoodEnrichmentClusters()
                * Assign the reference cells and run additional calculations.
                * esqAn.assignReferenceCells()
                * esqAn.calculateLeidenExpressionSignatures()
                * esqAn.assignCellTypesOnExpression()
                * esqAn.plotNHoodEnrichmentClusters()
                * Runs on the stored adata.
        """

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

            self.spatialScatter(graphs=["Cluster"], groups=inst_clusters, size=size, figsize=figsize)

        return None

    def gr_centralityScores(self, cluster_key="leiden"):
        """
            Documentation

            gr_centralityScores() class function docs
            =============================
            ----

            * gr_centralityScores() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools gr_centralityScores() function.
            * Compute centrality scores per cluster or cell type.

            ----

            Parameters:
            ----------
            * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.

            Returns:
            -------
            * The value of sq.gr.centrality_scores().

            Examples:
            --------
            1. Calling gr_centralityScores()
                * esqAn.gr_centralityScores()
                * Runs on the stored adata using default "leiden".
        """

        return tools.gr_centralityScores(adata=self.getAdata(), cluster_key=cluster_key)

    # graph the centrality scores
    def pl_centralityScores(self, cluster_key="leiden", figsize=None):
        """
            Documentation

            pl_centralityScores() class function docs
            =============================
            ----

            * pl_centralityScores() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools pl_centralityScores() function.
            * Plot the centrality scores that were calculated with gr_centralityScores()

            ----

            Parameters:
            ----------
            * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.
            * figsize=(10, 10): Size of the plot

            Returns:
            -------
            * The value of sq.pl_centralityScores().

            Examples:
            --------
            1. Calling pl_centralityScores()
                * esqAn.pl_centralityScores()
                * Runs on the stored adata using default "leiden".
                * Plots can be shown with esqAn.showPlots()
        """

        return tools.pl_centralityScores(adata=self.getAdata(), cluster_key=cluster_key, figsize=figsize)

    def adataSubsample(self, fraction=0.5, copy=True):
        """
            Documentation

            adataSubsample() class function docs
            =============================
            ----

            * adataSubsample() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools adataSubsample() function.
            * Subsample to a fraction of the number of observations.

            ----

            Parameters:
            ----------
            * fraction=0.5: Fraction to subsample to.
            * copy=True: Determines whether a copy is returned.

            Returns:
            -------
            * The value of sc.pp.subsample().
            * This will be the subsample that was computed.

            Examples:
            --------
            1. Calling adataSubsample()
                * subsample = esqAn.adataSubsample()
                * Runs on the stored adata.
        """

        return tools.adataSubsample(adata=self.getAdata(), fraction=fraction, copy=copy)

    def gr_co_occurrence(self, adata=None, cluster_key="leiden"):
        """
            Documentation

            gr_co_occurrence() class function docs
            =============================
            ----

            * gr_co_occurrence() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools gr_co_occurrence() function.
            * Calculates the co-occurrence probabilities for clusters.

            ----

            Parameters:
            ----------
            * adata=adataSubsample: A subsample of adata we can run on. If running on whole set, this should be None.
            * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.

            Returns:
            -------
            * The value of sq.gr.co_occurrence().

            Examples:
            --------
            1. Calling gr_co_occurrence()
                * esqAn.gr_co_occurrence()
                * Runs on the stored adata.
            ..
            2. Calling gr_co_occurrence()
                * subsample = esqAn.adataSubsample()
                * esqAn.gr_co_occurrence(adata=subsample)
                * Runs on the provided subsample.
        """

        if adata is None:
            return tools.gr_co_occurrence(adata=self.getAdata(), cluster_key=cluster_key)

        else:
            return tools.gr_co_occurrence(adata=adata, cluster_key=cluster_key)

    def pl_co_occurrence(self, adata=None, cluster_key="leiden", clusters="12", figsize=None):
        """
            Documentation

            pl_co_occurrence() class function docs
            =============================
            ----

            * pl_co_occurrence() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools pl_co_occurrence() function.
            * Plot the co-occurrence probability ratio for each cluster.

            ----

            Parameters:
            ----------
            * adata=adataSubsample: A subsample of adata we can run on. If running on whole set, this should be None.
            * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.
            * clusters="12": Cluster instances for which to plot conditional probability.
            * figsize=(10,10): Size of the plot.

            Returns:
            -------
            * The value of sq.pl.co_occurrence().

            Examples:
            --------
            1. Calling pl_co_occurrence()
                * esqAn.pl_co_occurrence()
                * Runs on the stored adata.
            ..
            2. Calling pl_co_occurrence()
                * subsample = esqAn.adataSubsample()
                * esqAn.pl_co_occurrence(adata=subsample)
                * Runs on the provided subsample.
        """

        if adata is None:
            return tools.pl_co_occurrence(adata=self.getAdata(), cluster_key=cluster_key, clusters=clusters,
                                          figsize=figsize)

        else:
            return tools.pl_co_occurrence(adata=adata, cluster_key=cluster_key, clusters=clusters, figsize=figsize)

    def gr_ripley(self, cluster_key="leiden", mode="L"):
        """
            Documentation

            gr_ripley() class function docs
            =============================
            ----

            * gr_ripley() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools gr_ripley() function.
            * Calculate various Ripley's statistics for point processes.

            ----

            Parameters:
            ----------
            * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.
            * mode="L": Which Ripley's statistic to compute. "F", "G", "L"

            Returns:
            -------
            * The value of sq.gr.ripley().

            Examples:
            --------
            1. Calling gr_ripley()
                * esqAn.gr_ripley()
                * Runs on the stored adata.
        """

        return tools.gr_ripley(adata=self.getAdata(), cluster_key=cluster_key, mode=mode)

    def pl_ripley(self, cluster_key="leiden", mode="L"):
        """
            Documentation

            pl_ripley() class function docs
            =============================
            ----

            * pl_ripley() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools pl_ripley() function.
            * Plot Ripley's statistics.

            ----

            Parameters:
            ----------
            * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.
            * mode="L": Which Ripley's statistic to compute. "F", "G", "L"

            Returns:
            -------
            * The value of sq.pl.ripley().

            Examples:
            --------
            1. Calling pl_ripley()
                * esqAn.pl_ripley()
                * Runs on the stored adata.
                * show the plots with esqAn.showPlots()
        """

        return tools.pl_ripley(adata=self.getAdata(), cluster_key=cluster_key, mode=mode)

    def gr_spatialAutocorr(self, adata=None, mode="moran", nPerms=None, nJobs=None):
        """
            Documentation

            gr_spatialAutocorr() class function docs
            =============================
            ----

            * gr_spatialAutocorr() class function for EasySQ class.
            * Acts as a wrapper for EasySQTools gr_spatialAutocorr() function.
            * Calculate Global Autocorrelation Statistic (Moran’s I or Geary’s C).

            ----

            Parameters:
            ----------
            * adata=subsample: Provided adata to run on. If running on self, adata=None.
            * mode="geary": Mode of score calculation. "moran" "geary"
            * nPerms=4: Number of permutations for the permutation test. If None, only p-values under normality assumption are computed.
            * nJobs=2: Number of parallel jobs.

            Returns:
            -------
            * The value of sq.gr.spatial_autocorr().

            Examples:
            --------
            1. Calling gr_spatialAutocorr()
                * esqAn.gr_spatialAutocorr()
                * Runs on the stored adata.
        """

        if adata is None:
            tools.gr_spatialAutocorr(adata=self.getAdata(), mode=mode, nPerms=nPerms, nJobs=nJobs)

        else:
            tools.gr_spatialAutocorr(adata=adata, mode=mode, nPerms=nPerms, nJobs=nJobs)

    def assignReferenceCells(self, gene_panel_url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41421-021-00266-1/MediaObjects/41421_2021_266_MOESM1_ESM.xlsx"):
        """
            Documentation

            assignReferenceCells() class function docs
            =============================
            ----

            * assignReferenceCells() class function for EasySQ class.
            * Assign marker gene metadata using reference dataset.
            * print/log The markers after computation.

            ----

            Parameters:
            ----------
            * gene_panel_url="https://content.com/reference_dataset.xlsx": Url of the panel we will use to reference cells.

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling assignReferenceCells()
                * esqAn.assignReferenceCells()
                * Runs on the stored adata.
        """

        # Assign marker gene metadata using reference dataset.

        # get the reference panel
        gene_panel = gene_panel_url
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

        return None

    def calculateLeidenExpressionSignatures(self):
        """
            Documentation

            calculateLeidenExpressionSignatures() class function docs
            =============================
            ----

            * calculateLeidenExpressionSignatures() class function for EasySQ class.
            * Calculates the leiden expression signatures.
            * Uses leiden, sig leiden, and meta leiden to calculate expression signatures.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling calculateLeidenExpressionSignatures()
                * esqAn.calculateLeidenExpressionSignatures()
                * Runs on the stored adata.
        """

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

        return None

    def assignCellTypesOnExpression(self):
        """
            Documentation

            assignCellTypesOnExpression() class function docs
            =============================
            ----

            * assignCellTypesOnExpression() class function for EasySQ class.
            * Assigns cell types based on expression.
            * This uses the reference panel set in assignReferenceCells()

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling assignCellTypesOnExpression()
                * esqAn.assignCellTypesOnExpression()
                * Runs on the stored adata.
        """

        meta_gene = pd.DataFrame(index=self.getSigLeiden().index.tolist())
        meta_gene["info"] = pd.Series("", index=meta_gene.index.tolist())
        meta_gene["Markers"] = pd.Series("N.A.", index=self.getSigLeiden().index.tolist())
        meta_gene.loc[self.getCommonMarkerGenes(), "Markers"] = self.getDfRefPanel().loc[
            self.getCommonMarkerGenes(), "Function"
        ]

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

        return None

    def calcSerCloseness(self):
        """
            Documentation

            calcSerCloseness() class function docs
            =============================
            ----

            * calcSerCloseness() class function for EasySQ class.
            * Calculates ser closeness value and sets the class variable.
            * closeness centrality - measure of how close the group is to other nodes.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling calcSerCloseness()
                * esqAn.calcSerCloseness()
                * Runs on the stored adata.
        """

        self.setSerCloseness(self.getDfCentral()["closeness_centrality"].sort_values(ascending=False))
        return None

    def calcSerDegree(self):
        """
            Documentation

            calcSerDegree() class function docs
            =============================
            ----

            * calcSerDegree() class function for EasySQ class.
            * Calculates ser degree value and sets the class variable.
            * The degree centrality for a node v is the fraction of nodes it is connected to.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling calcSerDegree()
                * esqAn.calcSerDegree()
                * Runs on the stored adata.
        """

        self.setSerDegree(self.getDfCentral()["degree_centrality"].sort_values(ascending=False))
        return None

    def calcSerCluster(self):
        """
            Documentation

            calcSerCluster() class function docs
            =============================
            ----

            * calcSerCluster() class function for EasySQ class.
            * Calculates ser cluster value and sets the class variable.
            * clustering coefficient - measure of the degree to which nodes cluster together.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling calcSerCluster()
                * esqAn.calcSerCluster()
                * Runs on the stored adata.
        """

        self.setSerCluster(self.getDfCentral()["average_clustering"].sort_values(ascending=False))
        return None

    def highClosenessScore(self):
        """
            Documentation

            highClosenessScore() class function docs
            =============================
            ----

            * highClosenessScore() class function for EasySQ class.
            * Plots the high end of the ser closeness score.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling highClosenessScore()
                * esqAn.highClosenessScore()
                * Runs on the stored adata.
                * Show the plot with esqAn.showPlots().
        """

        inst_clusters = self.getSerCloseness().index.tolist()[:5]
        self.spatialScatter(groups=inst_clusters, graphs="Cluster")
        return None

    def lowClosenessScore(self):
        """
            Documentation

            lowClosenessScore() class function docs
            =============================
            ----

            * lowClosenessScore() class function for EasySQ class.
            * Plots the low end of the ser closeness score.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling lowClosenessScore()
                * esqAn.lowClosenessScore()
                * Runs on the stored adata.
                * Show the plot with esqAn.showPlots().
        """

        inst_clusters = self.getSerCloseness().index.tolist()[-5:]
        self.spatialScatter(graphs="Cluster", groups=inst_clusters)
        return None

    def highDegreeCentrality(self):
        """
            Documentation

            highDegreeCentrality() class function docs
            =============================
            ----

            * highDegreeCentrality() class function for EasySQ class.
            * Plots the high end of the ser degree score.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling highDegreeCentrality()
                * esqAn.highDegreeCentrality()
                * Runs on the stored adata.
                * Show the plot with esqAn.showPlots().
        """

        inst_clusters = self.getSerDegree().index.tolist()[:5]
        self.spatialScatter(graphs="Cluster", groups=inst_clusters)
        return None

    def lowDegreeCentrality(self):
        """
            Documentation

            lowDegreeCentrality() class function docs
            =============================
            ----

            * lowDegreeCentrality() class function for EasySQ class.
            * Plots the low end of the ser degree score.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling lowDegreeCentrality()
                * esqAn.lowDegreeCentrality()
                * Runs on the stored adata.
                * Show the plot with esqAn.showPlots().
        """

        inst_clusters = self.getSerDegree().index.tolist()[-5:]
        self.spatialScatter(graphs="Cluster", groups=inst_clusters)
        return None

    def highClusteringCoefficient(self):
        """
            Documentation

            highClusteringCoefficient() class function docs
            =============================
            ----

            * highClusteringCoefficient() class function for EasySQ class.
            * Plots the high end of the ser cluster score.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling highClusteringCoefficient()
                * esqAn.highClusteringCoefficient()
                * Runs on the stored adata.
                * Show the plot with esqAn.showPlots().
        """

        inst_clusters = self.getSerCluster().index.tolist()[:5]
        self.spatialScatter(graphs="Cluster", groups=inst_clusters)
        return None

    def lowClusteringCoefficient(self):
        """
            Documentation

            lowClusteringCoefficient() class function docs
            =============================
            ----

            * lowClusteringCoefficient() class function for EasySQ class.
            * Plots the low end of the ser cluster score.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling lowClusteringCoefficient()
                * esqAn.lowClusteringCoefficient()
                * Runs on the stored adata.
                * Show the plot with esqAn.showPlots().
        """

        inst_clusters = self.getSerCluster().index.tolist()[-5:]
        self.spatialScatter(graphs="Cluster", groups=inst_clusters)
        return None

    def calcMoransIScore(self, numView=12):
        """
            Documentation

            calcMoransIScore() class function docs
            =============================
            ----

            * calcMoransIScore() class function for EasySQ class.
            * Calculates the top and bottom autocorrelation using gr_spatialAutocorr()
            * Sets the "topAutoCorr" and "botAutoCorr" class variables based on numView.
            * These sort the adata "moranI" values in descending and ascending order respectively.

            ----

            Parameters:
            ----------
            * numView=12: Number of top and bottom autocorr we calculate for.

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
        """
            Documentation

            setLeidenColors() class function docs
            =============================
            ----

            * setLeidenColors() class function for EasySQ class.
            * Sets the leiden colors used during spatial scatters and umaps. These are based on a file containing a list of hex values.

            ----

            Parameters:
            ----------
            * color_file="leiden_generated_random_3.txt": The color_file contained within the local "colors" directory.

            Returns:
            -------
            * leidenColors: A list of colors that leiden colors was set to.

            Examples:
            --------
            1. Calling setLeidenColors()
                * color_file = "leiden_generated_random_3.txt"
                * esqAn.setLeidenColors(color_file=color_file)
                * Runs on the stored adata.
            ..
            2. Calling setLeidenColors()
                * color_file = "colors/leiden_generated_random_3.txt"
                * esqAn.setLeidenColors(color_file=color_file)
                * Does the same thing as example 1.
        """

        leidenColors = tools.getColors(color_file=color_file)
        tools.setLeidenColors(adata=self.getAdata(), colors=leidenColors)
        return leidenColors

    def showPlots(self):
        """
            Documentation

            showPlots() class function docs
            =============================
            ----

            * showPlots() class function for EasySQ class.
            * Used to show the plots calculated by other methods / functions.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling showPlots()
                * esqAn.showPlots()
        """

        tools.showPlots()
        return None

    ####################################################################################################################
    ####################################################################################################################
    # region print functions
    def print(self):
        """
            Documentation

            print() class function docs
            =============================
            ----

            * print() class function for EasySQ class.
            * Prints some basic information about the AnnData object.

            ----

            Parameters:
            ----------
            * None

            Returns:
            -------
            * None

            Examples:
            --------
            1. Calling print()
                * esqAn.print()
        """

        print("\nOBJECT INFO")
        print("object: {}".format(self))
        print("object ID: {}".format(self.getID()))

        print("\nDATA")
        print(" data path: {}".format(self.getDataPath()))
        print(" adata:\n {}".format(self.getAdata()))

        print("")
        return None
    # endregion


if __name__ == "__main__":
    pass

    """
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
    # """
