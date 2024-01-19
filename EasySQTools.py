# Sam Crouse
# scrouse2@uwyo.edu
import copy
# Tools which make it easier to work with SquidPy

# python imports
import logging
import random
import string
import os
import matplotlib.pyplot as plt
import warnings
import pandas as pd
from copy import deepcopy
import seaborn as sns

logging.captureWarnings(True)
import squidpy as sq
import scanpy as sc

logging.captureWarnings(False)

"""
    Documentation

    EasySQTools() Library Function Docs
    =============================
    ----

    * Holds many of the non-class functions used in EasySQ.
    * Most of the class functions act as wrappers or abstraction layers for these functions.
    * If you don't wish (for some reason) to use EasySQ Analysis class functions you can use these directly.
    * Many of these functions act as a squidpy abstraction layer. However, many of these functions provide additional functionality. Such as the ability to produce custom color palettes and check for errors.

    ----

    Examples:
    --------
    1. Importing EasySQTools.py:
        * import EasySQTools as tools
        * tools.createCode(6)
        * tools.spatialScatter(adata=adata, graphs=["leiden"])
"""


def createCode(codeLength):
    """
        Documentation

        createCode() function docs
        =============================
        ----

        * createCode() tools function
        * Produce a random nth digit code. Can be used as an extra identifier to track classes, objects, or sms messages, etc.

        ----

        Parameters:
        ----------
        * codeLength=6: The number of characters the code will have.

        Returns:
        -------
        * String: "ab5lS3L": The code that was created.

        Examples:
        --------
        1. Calling createCode():
            * code = tools.createCode(6)
            * print(code) => "ab5lS3L"
    """

    code = ""

    i = 0
    while i < codeLength:
        numberOrLetter = random.randint(0, 1)
        if numberOrLetter == 0:  # add number
            code += str(random.randint(0, 9))

        else:  # add letter
            string.ascii_letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
            code += random.choice(string.ascii_letters)

        i += 1

    return code


def readVizgen(data_path, transform_path=None, transform_file='micron_to_mosaic_pixel_transform.csv'):
    """
        Documentation

        readVizgen() function docs
        =============================
        ----

        * readVizgen() lower level tools function for EasySQ
        * Read in and create adata based on metadata files.
        This looks through your data directory to find files that contain 'cell_metadata' and 'cell_by_gene' in their
        names. If they do, then it decides that they are your metadata files. If this is not the case,
        you can import using squidpy and set the class object using esqAn.setAdata(). Or, you can rename your files by
        adding 'cell_metadata' at the end of your cell metadata and 'cell_by_gene' at the end of your cell by gene file.

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
            * adata = tools.readVizgen(data_path="path/to/my/data/directory/")
    """

    meta_data_path = 'cell_metadata'
    cell_by_gene_path = 'cell_by_gene'

    meta_data_file = searchFiles(data_path, meta_data_path, loose=True)
    cell_by_gene_file = searchFiles(data_path, cell_by_gene_path, loose=True)

    foundMetaData = False
    foundCellByGene = False
    if meta_data_file is not None:
        foundMetaData = True

    if cell_by_gene_file is not None:
        foundCellByGene = True

    if meta_data_file is None:
        execString = "no cell_metadata file found"
        raise FileNotFoundError(execString)

    if cell_by_gene_file is None:
        execString = "no cell_by_gene file found"
        raise FileNotFoundError(execString)

    # print(meta_data_file)
    # print(cell_by_gene_file)

    if foundMetaData and foundCellByGene:
        if transform_path is None:
            adata = sq.read.vizgen(
                path=data_path,
                counts_file=cell_by_gene_file,
                meta_file=meta_data_file
            )

            return adata

        else:
            adata = sq.read.vizgen(
                path=transform_path,
                counts_file=cell_by_gene_file,
                meta_file=meta_data_file,
                transformation_file=transform_file
            )

            return adata


def adataSetup(adata):
    """
        Documentation

        adataSetup() tools function docs
        =============================
        ----

        * adataSetup() lower level tools function for EasySQ
        * Run basic analysis on the passed adata object.
        This includes printing the available graphs before the calculation and then running, qcMetrics(), layers(),
        highlyVariableGenes(), normalizeTotal(), log1p(), pp_pca(), pp_neighbors(), tl_umap(), leiden(), and then
        prints the availableGraphs() after finishing. Runs all with default parameters on the provided adata.

        ----

        Parameters:
        ----------
        * adata: AnnData object: The AnnData object you wish to do the setup analysis on.

        Returns:
        -------
        * None

        Examples:
        --------
        1. Calling adataSetup()
            * tools.setupAdata(adata=adata)
    """

    print("start setup calculations")
    print("BEFORE:")
    availableGraphs(adata)

    print("\ncalculating QC metrics...")
    perUnassigned = qcMetrics(adata)
    print("percent unassigned {}%".format(perUnassigned))

    print("\nrunning a bunch of calculations...")
    print("layers...")
    layers(adata)
    print("\nhighly variable genes...")
    highlyVariableGenes(adata)
    print("\nnormalize total...")
    normalizeTotal(adata)
    print("\nlog1p...")
    log1p(adata)
    print("\nprinciple component analysis...")
    pp_pca(adata)
    print("\nneighbors...")
    pp_neighbors(adata)
    print("\ncalculate UMAP...")
    tl_umap(adata)
    print("\nleiden...")
    leiden(adata)

    print("AFTER:")
    availableGraphs(adata)
    print("finished setup calculations\n")

    return None


# region adataSetup functions
# calculates quality control metrics and returns the percent of genes unassigned
def qcMetrics(adata, percentTop=(50, 100)):
    """
        Documentation

        qcMetrics() tools function docs
        =============================
        ----

        * Tools function for EasySQ class function qcMetrics.
        * Acts as a wrapper for the sc.pp.calculate_qc_metrics.
        * Calculates a number of qc metrics for an AnnData object.

        ----

        Parameters:
        ----------
        * adata: The AnnData object to calculate on.
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
            * "import EasySQTools as tools"
            * tools.qcMetrics(adata)
            * **If possible call this from its wrapper using EasySQ class "Analysis."**
    """

    try:
        makeVarNamesUnique(adata=adata)

        sc.pp.calculate_qc_metrics(adata, percent_top=percentTop, inplace=True, qc_vars=["mt"])
        perUnassigned = adata.obsm["blank_genes"].to_numpy().sum() / adata.var["total_counts"].sum() * 100

        return perUnassigned

    except Exception as e:  # yes, I know it's bad to do this.
        print("unknown error, running without qc_vars")
        print("ERROR: {}".format(e))

    sc.pp.calculate_qc_metrics(adata, percent_top=percentTop, inplace=True)
    perUnassigned = adata.obsm["blank_genes"].to_numpy().sum() / adata.var["total_counts"].sum() * 100

    return perUnassigned


def makeVarNamesUnique(adata):
    """
        Documentation

        makeVarNamesUnique() tools function docs
        =============================
        ----

        * Lower level tools function for EasySQ
        * Makes variable names unique.

        ----

        Parameters:
        ----------
        * adata: The AnnData object whose variable names are to be made unique.

        Returns:
        -------
        * None

        Examples:
        --------
        1. Calling makeVarNamesUnique(adata):
            * tools.makeVarNamesUnique(adata)
    """

    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("mt-")

    return None


def layers(adata):
    """
        Documentation

        layers() tools function docs
        =============================
        ----

        * Tools function for EasySQ class function layers().
        * Returns the layer named "counts" and sets it to the value of the adata variables.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run layers() on.

        Returns:
        -------
        * None

        Examples:
        --------
        1. Calling layers()
            * tools.layers(adata)
            * Runs layers() on the provided adata.
    """

    adata.layers["counts"] = adata.X.copy()
    return None


def highlyVariableGenes(adata, nTopGenes=4000):
    """
        Documentation

        highlyVariableGenes() tools function docs
        =============================
        ----

        * highlyVariableGenes() function for EasySQ class function highlyVariableGenes().
        * Acts as a wrapper for the Scanpy sc.pp.highly_variable_genes().
        * Returns the Scanpy sc.pp.highly_variable_genes() function value.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * nTopGenes=4000: The number of highly variable genes to keep

        Returns:
        -------
        * Returns the value of sc.pp.highly_variable_genes(). Updates .var or returns depending on "inplace"

        Examples:
        --------
        1. Calling highlyVariableGenes()
            * tools.highlyVariableGenes(adata, nTopGenes=3000)
            * Runs on the provided adata.
    """

    return sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=nTopGenes)


def normalizeTotal(adata, inplace=True):
    """
        Documentation

        normalizeTotal() tools function docs
        =============================
        ----

        * normalizeTotal() function for EasySQ class function normalizeTotal().
        * Acts as a wrapper for the Scanpy sc.pp.normalize_total() function.
        * Normalizes the counts per cells.
        * Returns the Scanpy sc.pp.normalize_total() function value.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * inPlace=True: Whether to perform in place or return the data calculated.

        Returns:
        -------
        * Returns the value of sc.pp.normalize_total(). Updates or returns depending on "inplace."

        Examples:
        --------
        1. Calling normalizeTotal()
            * tools.normalizeTotal(adata)
            * Runs on the provided adata.
    """

    return sc.pp.normalize_total(adata, inplace=inplace)


def log1p(adata):
    """
        Documentation

        log1p() tools function docs.
        =============================
        ----

        * log1p() function for EasySQ class function log1p().
        * Acts as a wrapper for the Scanpy sc.pp.log1p() function.
        * Logarithmize the data matrix.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.

        Returns:
        -------
        * Returns the value of sc.pp.log1p().

        Examples:
        --------
        1. Calling log1p()
            * tools.log1p(adata)
            * Runs on the stored adata.
    """

    return sc.pp.log1p(adata)


def pp_pca(adata):
    """
        Documentation

        pp_pca() tools function docs
        =============================
        ----

        * pp_pca() function for EasySQ class function pp_pca().
        * Acts as a wrapper for the Scanpy sc.pp.pca() function.
        * Runs principle component analysis on the provided adata.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.

        Returns:
        -------
        * Returns the value of sc.pp.pca().

        Examples:
        --------
        1. Calling pp_pca()
            * tools.pp_pca(adata)
            * Runs on the provided adata.
    """

    return sc.pp.pca(adata)


def tl_pca(adata, svdSolver="arpack"):
    """
        Documentation

        tl_pca() tools function docs
        =============================
        ----

        * tl_pca() function for EasySQ class function tl_pca().
        * Acts as a wrapper for the Scanpy sc.tl.pca() function.
        * Runs principle component analysis on the stored adata.
        * The difference between tl_pca and pl_pca is the function that it wraps.
        * Some demos use pl and tl in different scenarios.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * svdSolver="arpack": The SVD solver to use.

        Returns:
        -------
        * Returns the value of sc.tl.pca().

        Examples:
        --------
        1. Calling tl_pca()
            * tools.tl_pca(adata)
            * Runs on the provided adata.
    """

    return sc.tl.pca(adata, svd_solver=svdSolver)


def pp_neighbors(adata, nNeighbors=15, nPcs=None):
    """
        Documentation

        pp_neighbors() tools function docs
        =============================
        ----

        * pp_neighbors() function for EasySQ class function pp_neighbors().
        * Acts as a wrapper for the Scanpy sc.pp.neighbors() function.
        * Computes a neighborhood graph of observations.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * nNeighbors=10: The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation.
        * nPcs=20: Use this many PCs. If n_pcs==0 use .X if use_rep is None.

        Returns:
        -------
        * Returns the value of sc.pp.neighbors().

        Examples:
        --------
        1. Calling pp_neighbors()
            * tools.pp_neighbors(adata)
            * Runs on the provided adata.
    """

    return sc.pp.neighbors(adata, n_neighbors=nNeighbors, n_pcs=nPcs)


# UMAP calculations functions
def tl_umap(adata):
    """
        Documentation

        tl_umap() tools function docs
        =============================
        ----

        * tl_umap() function for EasySQ class function tl_umap().
        * Acts as a wrapper for the Scanpy sc.tl.umap() function.
        * Computes a neighborhood graph of observations.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.

        Returns:
        -------
        * Returns the value of sc.tl.umap().

        Examples:
        --------
        1. Calling tl_umap()
            * tools.tl_umap(adata)
            * Runs on the provided adata.
    """

    return sc.tl.umap(adata)


# UMAP plot function
def pl_umap(adata, graphs=["leiden"], show=False, size=None, wspace=0.4):
    """
        Documentation

        pl_umap() tools function docs
        =============================
        ----

        * pl_umap() function for EasySQ class function pl_umap().
        * Acts as a wrapper for the Scanpy sc.pl.umap() function.
        * Similar to spatialScatter(). Takes graphs and then plots umap based on those graphs.
        * Will autogenerate a color palette if the current one is too small. Automatically reduces size of palette if too large.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
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
            * tools.pl_umap(adata, graphs=["leiden", "Slc17a7"])
            * Runs on the provided adata.
            * Plots can be shown with tools.showPlots()
    """

    tryCount = 0  # tracks the number of retries so we can limit the number of retries
    totalTries = 10  # default of 2 allowed tries
    loopForColors = True
    while loopForColors:
        try:
            # check for retries
            if tryCount >= totalTries:
                break

            tryCount += 1

            sc.pl.umap(adata, color=graphs, size=size, wspace=wspace, show=False)

            containsLeiden = False  # see if the graphs contain leiden clusters
            # check if leiden is being graphed
            for graph in graphs:
                if graph.lower() == "leiden":  # check if leiden is in graphs to plot
                    containsLeiden = True

                    # If the colors have been set, then break and return
                    if getLeidenColors(adata) is not None:
                        loopForColors = False

                    # otherwise, close the grey graph and do a color init.
                    else:
                        print("leiden colors are None, closing graph and running an init.")

                        # clear the empty graph
                        figToClose = plt.get_fignums()[-1]
                        plt.close(figToClose)

                        leidenColorInit(adata=adata)

                    break  # break after leiden has been found

            if not containsLeiden:
                loopForColors = False

        # this prevents the error where we expected 32 items and got 33 or similar palette error.
        except ValueError as e:
            print(e)
            # default leiden color initialization
            for graph in graphs:
                if graph.lower() == "leiden":  # check if leiden is in graphs to plot
                    print("VALUE ERROR: Palette size mismatch, correcting now")

                    # clear the empty graph
                    figToClose = plt.get_fignums()[-1]
                    plt.close(figToClose)

                    leidenColorInit(adata=adata)

                    break

            # check for retries
            if tryCount >= totalTries:
                break

            tryCount += 1

        # this usually occurs when leiden graph is trying to be plotted, but leiden() hasn't been run yet.
        except KeyError as e:
            if str(e).find("Could not find key leiden in .var_names or .obs.columns.") != -1:
                raise KeyError("leiden not found in adata. Are you trying to plot without running leiden() first?")

            elif str(e).find("Could not find 'umap' or 'X_umap' in .obsm") != -1:
                raise KeyError("Could not find 'umap' in .obsm, have you run tl_umap() first?")

            else:
                raise KeyError(e)

    if show:
        plt.show()

    return None


# calculate spatial neighbors data
def gr_spatialNeighbors(adata, coordType="generic", spatialKey="spatial", delaunay=False):
    """
        Documentation

        gr_spatialNeighbors() tools function docs
        =============================
        ----

        * gr_spatialNeighbors() function for EasySQ class function gr_spatialNeighbors().
        * Acts as a wrapper for the Squidpy sq.gr.spatial_neighbors() function.
        * Creates a graph using spatial coordinates.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * coordType="generic": Type of coordinate system. 'grid' 'generic'
        * spatialKey="spatial": Key in anndta.AnnData.obsm where spatial coordinates are stored. You probably set this earlier.
        * delaunay=False: Whether to compute the graph from Delaunay triangulation. Only used when coord_type='generic'.

        Returns:
        -------
        * The value of sq.gr.spatial_neighbors().

        Examples:
        --------
        1. Calling gr_spatialNeighbors()
            * tools.gr_spatialNeighbors(adata=adata)
            * Runs on the provided adata.
    """

    return sq.gr.spatial_neighbors(adata=adata, coord_type=coordType, spatial_key=spatialKey, delaunay=delaunay)


# calculate nhoodEnrichment
def gr_nhoodEnrichment(adata, clusterKey="leiden"):
    """
        Documentation

        gr_nhoodEnrichment() tools function docs
        =============================
        ----

        * gr_nhoodEnrichment() function for EasySQ class function gr_nhoodEnrichment().
        * Acts as a wrapper for the Squidpy sq.gr.nhood_enrichment() function.
        * Compute neighborhood enrichment by permutation test.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.

        Returns:
        -------
        * The value of sq.gr.nhood_enrichment().

        Examples:
        --------
        1. Calling gr_nhoodEnrichment()
            * tools.gr_nhoodEnrichment(adata)
            * Runs on the provided adata using default "leiden".
    """

    return sq.gr.nhood_enrichment(adata=adata, cluster_key=clusterKey)


# plot nhoodEnrichment data: No return value
def pl_nhoodEnrichment(adata, show=False, clusterKey="leiden", method="average", cmap="inferno", vmin=-50, vmax=100,
                       figsize=(5, 5)):
    """
        Documentation

        pl_nhoodEnrichment() tools function docs
        =============================
        ----

        * pl_nhoodEnrichment() function for EasySQ class function pl_nhoodEnrichment().
        * Acts as a wrapper for the Squidpy sq.pl.nhood_enrichment() function.
        * Plot the neighborhood enrichment.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
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
            * tools.pl_nhoodEnrichment(adata)
            * Runs on the provided adata.
        ..
        2. Calling pl_nhoodEnrichment() with show=True
            * tools.pl_nhoodEnrichment(adata, show=True)
            * Runs on the provided adata and displays the plots immediately.
    """

    sq.pl.nhood_enrichment(adata=adata, cluster_key=clusterKey, method=method, cmap=cmap, vmin=vmin,
                           vmax=vmax, figsize=figsize)

    if show:
        plt.show()

    return None


# compute the centrality scores
def gr_centralityScores(adata, cluster_key="leiden"):
    """
        Documentation

        gr_centralityScores() tools function docs
        =============================
        ----

        * gr_centralityScores() function for EasySQ class function gr_centralityScores().
        * Acts as a wrapper for the Squidpy sq.gr.centrality_scores() function.
        * Compute centrality scores per cluster or cell type.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.

        Returns:
        -------
        * The value of sq.gr.centrality_scores().

        Examples:
        --------
        1. Calling gr_centralityScores()
            * tools.gr_centralityScores(adata)
            * Runs on the provided adata using default "leiden".
    """

    return sq.gr.centrality_scores(adata=adata, cluster_key=cluster_key)


# graph the centrality scores
def pl_centralityScores(adata, cluster_key="leiden", figsize=None):
    """
        Documentation

        pl_centralityScores() tools function docs
        =============================
        ----

        * pl_centralityScores() function for EasySQ class function pl_centralityScores().
        * Acts as a wrapper for the Squidpy sq.pl.centrality_scores() function.
        * Plot the centrality scores that were calculated with gr_centralityScores()

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.
        * figsize=(10, 10): Size of the plot

        Returns:
        -------
        * The value of sq.pl_centralityScores().

        Examples:
        --------
        1. Calling pl_centralityScores()
            * tools.pl_centralityScores(adata)
            * Runs on the provided adata using default "leiden".
            * Plots can be shown with esqAn.showPlots()
    """

    return sq.pl.centrality_scores(adata=adata, cluster_key=cluster_key, figsize=figsize)


# Subsample to a fraction of the number of observations.
def adataSubsample(adata, fraction=0.5, copy=True):
    return sc.pp.subsample(adata, fraction=fraction, copy=copy)


# compute the co-occurrence probability
def gr_co_occurrence(adata, cluster_key="leiden"):
    """
        Documentation

        gr_co_occurrence() tools function docs
        =============================
        ----

        * gr_co_occurrence() function for EasySQ class function gr_co_occurrence().
        * Acts as a wrapper for the Squidpy sq.gr.co_occurrence() function.
        * Calculates the co-occurrence probabilities for clusters.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.

        Returns:
        -------
        * The value of sq.gr.co_occurrence().

        Examples:
        --------
        1. Calling gr_co_occurrence()
            * tools.gr_co_occurrence(adata)
            * Runs on the provided adata.
    """

    return sq.gr.co_occurrence(adata=adata, cluster_key=cluster_key)


# graph the co-occurrence probability
def pl_co_occurrence(adata, cluster_key="leiden", clusters="12", figsize=None):
    """
        Documentation

        pl_co_occurrence() tools function docs
        =============================
        ----

        * pl_co_occurrence() function for EasySQ class function pl_co_occurrence().
        * Acts as a wrapper for the Squidpy sq.pl.co_occurrence() function.
        * Plot the co-occurrence probability ratio for each cluster.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.
        * clusters="12": Cluster instances for which to plot conditional probability.
        * figsize=(10,10): Size of the plot.

        Returns:
        -------
        * The value of sq.pl.co_occurrence().

        Examples:
        --------
        1. Calling pl_co_occurrence()
            * tools.pl_co_occurrence(adata)
            * Runs on the provided adata.
    """

    return sq.pl.co_occurrence(adata=adata, cluster_key=cluster_key, clusters=clusters, figsize=figsize)


# compute Ripley's statistics
def gr_ripley(adata, cluster_key="leiden", mode="L"):
    """
        Documentation

        gr_ripley() tools function docs
        =============================
        ----

        * gr_ripley() function for EasySQ class function gr_ripley().
        * Acts as a wrapper for the Squidpy sq.gr.ripley() function.
        * Calculate various Ripley's statistics for point processes.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.
        * mode="L": Which Ripley's statistic to compute. "F", "G", "L"

        Returns:
        -------
        * The value of sq.gr.ripley().

        Examples:
        --------
        1. Calling gr_ripley()
            * tools.gr_ripley(adata)
            * Runs on the provided adata.
    """

    return sq.gr.ripley(adata, cluster_key=cluster_key, mode=mode)


# graph Ripley's statistics
def pl_ripley(adata, cluster_key="leiden", mode="L"):
    """
        Documentation

        pl_ripley() tools function docs
        =============================
        ----

        * pl_ripley() function for EasySQ class function pl_ripley().
        * Acts as a wrapper for the Squidpy sq.pl.ripley() function.
        * Plot Ripley's statistics.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * clusterKey="leiden": Key in anndata.AnnData.obs where clustering is stored. This almost always runs on leiden.
        * mode="L": Which Ripley's statistic to compute. "F", "G", "L"

        Returns:
        -------
        * The value of sq.pl.ripley().

        Examples:
        --------
        1. Calling pl_ripley()
            * tools.pl_ripley(adata)
            * Runs on the provided adata.
            * show the plots with esqAn.showPlots()
    """

    return sq.pl.ripley(adata, cluster_key=cluster_key, mode=mode)


def gr_spatialAutocorr(adata, mode="moran", nPerms=None, nJobs=None):
    """
        Documentation

        gr_spatialAutocorr() tools function docs
        =============================
        ----

        * gr_spatialAutocorr() function for EasySQ class function pl_ripley().
        * Acts as a wrapper for the Squidpy sq.gr.spatial_autocorr() function.
        * Calculate Global Autocorrelation Statistic (Moran’s I or Geary’s C).

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * mode="geary": Mode of score calculation. "moran" "geary"
        * nPerms=4: Number of permutations for the permutation test. If None, only p-values under normality assumption are computed.
        * nJobs=2: Number of parallel jobs.

        Returns:
        -------
        * The value of sq.gr.spatial_autocorr().

        Examples:
        --------
        1. Calling gr_spatialAutocorr()
            * tools.gr_spatialAutocorr(adata)
            * Runs on the provided adata.
    """

    return sq.gr.spatial_autocorr(adata, mode=mode, n_perms=nPerms, n_jobs=nJobs)


def leiden(adata, resolution=1.0, ignore=False):
    """
        Documentation

        leiden() tools function docs
        =============================
        ----

        * leiden() function for EasySQ class function leiden().
        * Acts as a wrapper for the Scanpy sc.tl.leiden() function.
        * Clusters cells into cell groups using the Leiden algorithm.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * resolution=1.0: A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters.
        * ignore=False: Parameter for determining error throwing if neighbors hasn't been run. Will run neighbors if neighbors hasn't been run, if this still produces an error it will throw it.

        Returns:
        -------
        * Returns the value of sc.tl.leiden().

        Examples:
        --------
        1. Calling leiden()
            * tools.leiden(adata)
            * Runs on the provided adata.
    """

    try:
        return sc.tl.leiden(adata, resolution=resolution)

    except KeyError as e:
        e = str(e).strip().replace("'", '')
        catchStr = 'No "neighbors" in .uns'

        if e == catchStr and not ignore:
            warnNeighborsStr = "No neighbors found, running 'neighbors()' default now! " \
                               "Run neighbors() with your args first if this is not what you want!"
            warnings.warn("Warning................................\n{}".format(warnNeighborsStr))

            pp_neighbors(adata)
            return leiden(adata, ignore=True)  # call leiden again with ignore, so if we get this error we will raise

        else:
            raise


# can be used to filter cells with low expression
def filterCells(adata, minCounts=50):
    """
        Documentation

        filterCells() tools function docs
        =============================
        ----

        * filterCells() function for EasySQ class function filterCells().
        * Acts as a wrapper for the Scanpy sc.pp.filter_cells() function.
        * Filter cell outliers based on counts and numbers of genes expressed.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * minCounts=50: The minimum number of genes expressed for a cell to pass filtering.

        Returns:
        -------
        * Returns the value of sc.pp.filter_cells().

        Examples:
        --------
        1. Calling filterCells()
            * tools.filterCells(adata)
            * Runs on the provided adata.
    """

    return sc.pp.filter_cells(adata, min_counts=minCounts)


# can be used to filter genes that are expressed in too few cells
def filterGenes(adata, minCells=10):
    """
        Documentation

        filterGenes() tools function docs
        =============================
        ----

        * filterGenes() function for EasySQ class function filterGenes().
        * Acts as a wrapper for the Scanpy sc.pp.filter_genes() function.
        * Filter genes based on number of cells or counts.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * minCells=10: The minimum number of cells expressed required for a gene to pass filtering.

        Returns:
        -------
        * Returns the value of sc.pp.filter_genes().

        Examples:
        --------
        1. Calling filterGenes()
            * tools.filterGenes(adata)
            * Runs on the provided adata.
    """

    return sc.pp.filter_genes(adata, min_cells=minCells)


def scale(adata, maxValue=10):
    """
        Documentation

        scale() tools function docs
        =============================
        ----

        * scale() function for EasySQ class function scale().
        * Acts as a wrapper for the Scanpy sc.pp.scale() function.
        * Scale data to unit variance and zero mean.
        * Can be used to scale gene expression. IE Clip values that exceed 10 ('max value') standard deviations.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * maxValue=10: Clip (truncate) to this value after scaling. If None, do not clip.

        Returns:
        -------
        * Returns the value of sc.pp.scale().

        Examples:
        --------
        1. Calling scale()
            * tools.scale(adata)
            * Runs on the provided adata.
    """

    return sc.pp.scale(adata, max_value=maxValue)


# plots transcript data
# requires running of qc metrics
def plotTranscripts(adata, show=False, figsize=(15, 4)):
    """
        Documentation

        plotTranscripts() tools function docs
        =============================
        ----

        * plotTranscripts() function for EasySQ class function plotTranscripts().
        * Plots calculated transcripts based on the provided AnnData object.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * show=False: Whether to display graphs on function call, or later by calling esqAn.showPlots()
        * figsize=(20, 6): The dimensions of the plot (x, y).

        Returns:
        -------
        * None

        Examples:
        --------
        1. Calling plotTranscripts()
            * tools.plotTranscripts(adata)
            * Runs on the provided adata.
    """

    fig, axs = plt.subplots(1, 4, figsize=figsize)

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

    if show:
        plt.show()

    return None


# endregion


# easier to use version of spatial_scatter than sq default
# takes in adata and a list of colors / filters IE ["leiden", "n_counts"]
#
# must setup adata by running adataSetup before this will work
def spatialScatter(adata, graphs, show=False, libraryID=None, wspace=None, size=None, shape=None,
                   groups=None, cmap=None, figsize=None):
    """
        Documentation

        spatialScatter() tools function docs
        =============================
        ----

        * spatialScatter() function for EasySQ class function spatialScatter().
        * Acts as a wrapper for the Squidpy sq.pl.spatial_scatter() function.
        * Produces a spatial scatter plot given a graph type and adata.
        * Contains a bit more functionality that the Squidpy spatial_scatter()
        This function deals with color palettes automatically for the user. This allows the user to plot things with any
        resolution. If the palette contains to few colors, it will generate a new palette automatically.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
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
            * tools.spatialScatter(adata, graphs=["leiden"])
            * Plots a spatial scatter of the type given by graphs based on the provided adata.
    """

    tryCount = 0  # tracks the number of retries so we can limit the number of retries
    totalTries = 10  # default of 2 allowed tries
    loopForColors = True
    while loopForColors:
        try:
            # check for retries
            if tryCount >= totalTries:
                break

            tryCount += 1

            sq.pl.spatial_scatter(
                adata,
                shape=shape,
                color=graphs,
                library_id=libraryID,
                wspace=wspace,
                size=size,
                groups=groups,
                cmap=cmap,
                figsize=figsize,
            )

            containsLeiden = False  # see if the graphs contain leiden clusters
            # check if leiden is being graphed
            for graph in graphs:
                if graph.lower() == "leiden":  # check if leiden is in graphs to plot
                    containsLeiden = True

                    # If the colors have been set, then break and return
                    if getLeidenColors(adata) is not None:
                        loopForColors = False

                    # otherwise, close the grey graph and do a color init.
                    else:
                        print("leiden colors are None, closing graph and running an init.")

                        # clear the empty graph
                        figToClose = plt.get_fignums()[-1]
                        plt.close(figToClose)

                        leidenColorInit(adata=adata)

                    break  # break after leiden has been found

            if not containsLeiden:
                loopForColors = False

        # this prevents the error where we expected 32 items and got 33 or similar palette error.
        except ValueError as e:
            print(e)
            # default leiden color initialization
            for graph in graphs:
                if graph.lower() == "leiden":  # check if leiden is in graphs to plot
                    print("VALUE ERROR: Palette size mismatch, correcting now")

                    # clear the empty graph
                    figToClose = plt.get_fignums()[-1]
                    plt.close(figToClose)

                    leidenColorInit(adata=adata)

                    break

            # check for retries
            if tryCount >= totalTries:
                break

            tryCount += 1

        # this usually occurs when leiden graph is trying to be plotted, but leiden() hasn't been run yet.
        except KeyError as e:
            if str(e).find("Could not find key leiden in .var_names or .obs.columns.") != -1:
                raise KeyError("leiden not found in adata. Are you trying to plot without running leiden() first?")

            else:
                raise KeyError(e)

    if show:
        plt.show()

    return None


# default leiden color initialization
def leidenColorInit(adata):
    """
        Documentation

        leidenColorInit() tools function docs
        =============================
        ----

        * Initializes the color palette for the adata.
        This function gets the current leiden colors, and if they are none it sets them to the default. Which is found
        in the colors directory under "leiden_generated_random_3.txt". If that palette is too small, it will generate a
        new palette.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.

        Returns:
        -------
        * list: Strings: A list of strings that the color palette was set to.

        Examples:
        --------
        1. Calling leidenColorInit()
            * tools.leidenColorInit(adata)
            * Runs on the provided adata.
        ..
        2. When is leidenColorInit() called?
            * Called in tools.spatialScatter() to set color palettes.
            * Called in tools.pl_umap() to set color palettes.
    """
    print("running leiden color init")

    # check if leiden colors have already been set
    currentLeidenColors = getLeidenColors(adata=adata)
    colors = currentLeidenColors

    colorLength = len(adata.obs["leiden"].value_counts())

    if colors is None:
        # get the default color set
        leidenColors = getColors('leiden_generated_random_5.txt')

    else:
        leidenColors = colors

    if len(leidenColors) < colorLength:
        print("Too few colors in palette! Generating and saving a longer one now!")
        leidenColors = createPalette(colorLength, save=True, log=True)

    if len(leidenColors) > colorLength:
        print("Too many colors in palette! Automatically reducing size now!")

    while len(leidenColors) > colorLength:
        leidenColors.pop(-1)

    return setLeidenColors(adata=adata, colors=leidenColors)


# sets the leiden colors based on given colors. If no colors, we generate them ourselves.
def setLeidenColors(adata, colors=None, length=200):
    """
        Documentation

        setLeidenColors() tools function docs
        =============================
        ----

        * Sets the color palette based on a list of string hex colors.
        * If no colors are provided, it generates a new palette of default length 200 colors.
        * It automatically saves the newly created palette in the colors directory.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata to run on.
        * colors=["#546326", "#33B1A8", "#F91316"]: A list of string hex colors defining the color palette.
        * length=200: The length of the newly generated color palette if no colors are provided.

        Returns:
        -------
        * list: Strings: A list of string hex colors that the color palette was set to.

        Examples:
        --------
        1. Calling setLeidenColors()
            * tools.setLeidenColors(adata)
            * Runs on the provided adata.
    """
    print("setting leiden colors")

    if colors is None:
        print("No colors given, generating and saving now!")
        leidenColors = createPalette(length=length, save=True, log=True)

        adata.uns['leiden_colors'] = leidenColors
        return leidenColors

    else:
        adata.uns['leiden_colors'] = colors
        return colors


def getLeidenColors(adata):
    """
        Documentation

        getLeidenColors() tools function docs
        =============================
        ----

        * Gets the currently set leiden color palette.
        * If the default palette is too small, it returns None.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata whose color palette we will get.

        Returns:
        -------
        * list: Strings: A list of string hex colors that defines the current leiden color palette.
        * None: If the default color palette is all grey. (Makes it useless)

        Examples:
        --------
        1. Calling getLeidenColors()
            * tools.getLeidenColors(adata)
            * Runs on the provided adata.
    """

    try:
        # check for default leiden_colors of '#808080' for each item
        for color in adata.uns['leiden_colors']:
            if str(color) != "#808080":
                # print(adata.uns['leiden_colors'])
                return adata.uns['leiden_colors']

        # print(None)
        return None

    except KeyError:
        # print(None)
        return None


def showPlots():
    """
        Documentation

        showPlots() tools function docs
        =============================
        ----

        * showPlots() function for EasySQ class function showPlots().
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
            * tools.showPlots()
    """

    plt.show()

    return None


def getColors(color_file):
    """
        Documentation

        getColors() tools function docs
        =============================
        ----

        * Search for a 'colors' directory and see if it contains the "color_file."
        * Defaults to the current working dir.
        * Defaults to searching in the "colors" directory which should be present in your file structure.

        ----

        Parameters:
        ----------
        * color_file="leiden_generated_random_3.txt": The file containing your palette. Should look like:
        #1DA309
        #AB5443
        #441AA5
        #EF6CE8
        #29E43A
        #36A271
        #8B025A
        #40E677
        #6DF7B9
        #3B06DA
        #50E5DA
        #FAA0E5
        #C34B73
        #C5FBEB
        #634BC8
        #1DBE46

        Returns:
        -------
        * list: Strings: A list of string hex colors that defines the current leiden color palette.

        Examples:
        --------
        1. Calling getColors()
            * tools.getColors("leiden_generated_random_3.txt")
    """

    # search for the given file
    if color_file.find("colors/") != -1:
        color_file = color_file.replace("colors/", '')

    color_path = str(os.getcwd()).replace('\\', '/') + '/'
    color_path = searchFiles(color_path, 'colors')

    # read the colors from the file and return them
    current_color_path = color_path + '/' + color_file
    with open(current_color_path, "r") as color_file:
        colors = color_file.read().split('\n')

    return colors


def createPalette(length, save=False, log=True):
    """
        Documentation

        createPalette() tools function docs
        =============================
        ----

        * Creates a random color palette of given length.
        * Optionally saves this palette to the "colors" directory.
        * Optionally logs the produced colors to the console.

        ----

        Parameters:
        ----------
        * length=200: The size of the color palette to produce.
        * save=False: Whether to save the produced palette to the "colors" directory or not.
        * log=True: Whether to log the returned colors or not.

        Returns:
        -------
        * list: Strings: A list of string hex colors that defines the current leiden color palette.

        Examples:
        --------
        1. Calling createPalette()
            * tools.createPalette(200)
            * Creates a random color palette containing 200 colors. Doesn't save it by default.
    """

    if type(length) != int:
        raise TypeError("length must be of type int")

    colors = []
    i = 0
    while i < length:
        color = random.randrange(0, 2 ** 24)
        hexColor = hex(color)

        if len(str(hexColor)) != 8:
            continue

        stdColor = "#" + str(hexColor[2:]).upper()
        colors.append(stdColor)

        i += 1

    # write to colors file
    if save:
        colors_dir = searchFiles(os.getcwd() + '\\', 'colors')

        if colors_dir is not None:
            # get number for save file
            color_files = os.listdir(colors_dir)
            # print(color_files)

            # search through the files so that we can give the generated colors the proper number
            file_num = 1
            save_file_template = 'leiden_generated_random_num.txt'
            for i in range(len(color_files)):
                # check through these files to see if any are auto generated, if so get the biggest file num
                if color_files[i].find(save_file_template.replace("num.txt", '')) != -1:
                    # get the file number 'num' from the currently generated files
                    current_file_num = color_files[i].replace("leiden_generated_random_", '')
                    current_file_num = int(current_file_num.replace(".txt", ''))

                    # prevents file duplication and means we create a file_num bigger than existing one
                    if current_file_num >= file_num:
                        file_num = current_file_num
                        file_num += 1

            save_file_name = save_file_template.replace("num", str(file_num))

            # save the file
            save_file = colors_dir + '\\' + save_file_name
            if log: print("saving to: {}".format(save_file))

            file_content = ""
            for color in colors:
                file_content += str(color) + '\n'

            file_content = file_content[:-1]  # remove the last '\n' from the file_content
            # print(file_content)

            with open(save_file, 'x') as f:
                f.write(file_content)

        else:
            raise FileNotFoundError('"colors" directory not found in current working directory. Can not save!')

    return colors


# searches a provided directory and subdirectories for a file or dir
# loose=bool lets us search for a name part instead of a full file name
def searchFiles(data_path, fileToFind, loose=False):
    """
        Documentation

        searchFiles() tools function docs
        =============================
        ----

        * Searches a provided directory and subdirectories for a file or directory.
        * Returns the path to the file or directory.

        ----

        Parameters:
        ----------
        * data_path="path/to/my/directoryToSearch/
        * fileToFind="theFileToFind.txt" or fileToFind="theDirToFind": The file or directory to return the path to.
        * loose=False: Should file we are looking for be the full name (false), or just part of the name of the file (true)?

        Returns:
        -------
        * String=".../users/myPythonProject/myDir/fileWeWantedToFind/: The path to the file or directory we searched for.

        Examples:
        --------
        1. Calling searchFiles()
            * path = os.getcwd().replace('\\', '/') + '/tutorial_data_1/'
            * cell_by_gene_path_test = 'cell_by_gene'
            * result = searchFiles(path, fileToFind=cell_by_gene_path_test, loose=True)
            * Returns the path to the file which has 'cell_by_gene' in its name within the cwd/tutorial_data_1/ directory.
    """

    data_path = copy.deepcopy(data_path)
    fileToFind = copy.deepcopy(fileToFind)
    dirList = os.listdir(data_path)

    dirs = []
    for file in dirList:
        file = str(file)

        if loose:
            if file.find(str(fileToFind)) != -1:
                return data_path + file

        else:
            if file == str(fileToFind):
                return data_path + file

        if file.find('.') == -1:  # this is a directory
            dirs.append(file)

    # print(dirs)

    for dir in dirs:
        # print(dir)
        new_path = data_path + dir + '/'

        returnVal = searchFiles(new_path, fileToFind, loose=loose)

        if returnVal is not None:
            return returnVal


# these are values that can be passed into spatialScatter() colors for graphing
def availableGraphs(adata, log=True):
    """
        Documentation

        availableGraphs() tools function docs
        =============================
        ----

        * availableGraphs() function for EasySQ class function availableGraphs().
        * Optionally log the graphs available to plot, then return them.

        ----

        Parameters:
        ----------
        * adata=AnnData object: The anndata object whose available graphs we would like to get.
        * log=True: Bool determining whether to print the available graphs to the console or not.

        Returns:
        -------
        * A list of strings which are the available graphs.

        Examples:
        --------
        1. Calling availableGraphs()
            * tools.availableGraphs(adata)
            * Logs the graphs available for the provided adata.
    """

    var_list = list(adata.var_names)
    obs_columns_list = list(adata.obs.columns)

    if log:
        print("AVAILABLE GRAPHS:\n Pass any of these to spatialScatter() to plot!\n Genes: {}\n Analyses: {}"
              "".format(var_list, obs_columns_list))

    graphs = var_list + obs_columns_list
    return graphs


# tests for EasySQ.py
if __name__ == "__main__":
    pass

    """
    # note: createCode(codeLength) tests
    testCode = createCode(8)
    print(testCode)
    # """

    """
    # note: readVizgen(data_path) tests
    path = 'F:/sunlabmerfishdata/QSFL01222023/'
    transformation_path = 'C:/Users/scrouse2/Crouse_Work_Files/PyCharm_Projects/EasySQ/'
    testAdata = readVizgen(path, transform_path=transformation_path)
    print(testAdata)

    # note: spatialScatter() tests and supporting funcs like getColors() and adataSetup()
    print("spatialScatter() and analysis tests")
    print(availableGraphs(testAdata))  # graphs available before setup
    spatialScatter(testAdata, ['volume', 'Th'])  # plot some of them
    adataSetup(testAdata)  # setup adata: produces additional graphs
    spatialScatter(testAdata, ['leiden'])  # run some spatial scatter tests
    # """

    """
    # note: searchFiles(data_path, fileToFind) tests
    path = 'F:/sunlabmerfishdata/QSFL01222023/'

    # shouldn't find this file
    print("\n\nun-found test")
    result = searchFiles(path, 'should not find this file.unfound')
    print("result: {}\n".format(result))

    # should find this test file
    print("found test")
    result = searchFiles(path, 'roi_metadata0_cat1.bin')
    print("result: {}\n\n".format(result))

    # loose testing
    path = os.getcwd().replace('\\', '/') + '/tutorial_data_1/'
    print(path)

    meta_data_path_test = 'cell_metadata'
    cell_by_gene_path_test = 'cell_by_gene'

    # should find this test file
    print("found test")
    result = searchFiles(path, fileToFind=meta_data_path_test, loose=True)
    print("result: {}\n\n".format(result))

    print("found test")
    result = searchFiles(path, fileToFind=cell_by_gene_path_test, loose=True)
    print("result: {}\n\n".format(result))
    # """

    """
    # note: color palette creation testing
    print("create palette")
    palette = createPalette(200, save=True)
    print(palette)

    # note: get colors testing
    # note: I like random_1, random_2, random_3
    print("get colors")
    testColors = getColors("leiden_generated_random_1.txt")
    print(testColors)

    # note: get and set LeidenColors testing
    path = 'F:/sunlabmerfishdata/QSFL01222023/'
    transformation_path = 'C:/Users/scrouse2/Crouse_Work_Files/PyCharm_Projects/EasySQ/'
    testAdata = readVizgen(path, transform_path=transformation_path)
    print(testAdata)

    setLeidenColors(adata=testAdata, colors=testColors)
    # setLeidenColors(adata=testAdata)
    print(getLeidenColors(adata=testAdata))

    # """
