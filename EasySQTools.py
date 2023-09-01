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


# produce a random nth digit code
# can be used as an extra identifier to track classes, objects, or sms messages
def createCode(codeLength):
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


# takes in a data_path: file with vizgen data, metadata, cellbygene data, etc.
# returns an anndata object based on this data
# optionally takes in a transform_path and file which are set to the sq.read.vizgen('transformation_file=') parameter.
#  transform_path is the path to the root directory which must hold an 'images' file which should contain a
#  transform_file called by default 'micron_to_mosaic_pixel_transform.csv'
#  path looks like this ...users/pythonproject/ images /micron_to_mosaic_pixel_transform.csv
#                         transform_path>leads to images>.csv
# 'images' should not be listed in any of the passed arguments
def readVizgen(data_path, transform_path=None, transform_file='micron_to_mosaic_pixel_transform.csv'):
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
        execString = "no meta_file found"
        raise FileNotFoundError(execString)

    if cell_by_gene_file is None:
        execString = "no cell_by_gene_file found"
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


# runs some basic scanpy functions to setup adata for graphing
def adataSetup(adata):
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

        qcMetrics() ".tools" function docs
        =============================
        ----

        * Tools function for EasySQ class function qcMetrics.
        * Acts as a wrapper for sc.pp.calculate_qc_metrics.
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

    except:  # yes, I know it's bad to do this.
        print("unknown error, running without qc_vars")

    sc.pp.calculate_qc_metrics(adata, percent_top=percentTop, inplace=True)
    perUnassigned = adata.obsm["blank_genes"].to_numpy().sum() / adata.var["total_counts"].sum() * 100

    return perUnassigned


# todo figure out a better method for this
def makeVarNamesUnique(adata):
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("mt-")


def layers(adata):
    adata.layers["counts"] = adata.X.copy()
    return


def highlyVariableGenes(adata, nTopGenes=4000):
    return sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=nTopGenes)


def normalizeTotal(adata, inplace=True):
    return sc.pp.normalize_total(adata, inplace=inplace)


def log1p(adata):
    return sc.pp.log1p(adata)


def pp_pca(adata):
    return sc.pp.pca(adata)


def tl_pca(adata, svdSolver="arpack"):
    return sc.tl.pca(adata, svd_solver=svdSolver)


def pp_neighbors(adata):
    return sc.pp.neighbors(adata)


# UMAP calculations functions
def tl_umap(adata):
    return sc.tl.umap(adata)


# UMAP plot function
def pl_umap(adata, graphs=["leiden"], show=False, size=None, wspace=0.4):
    loopForColors = True
    while loopForColors:
        # try catch for leiden existence
        try:
            sc.pl.umap(adata, color=graphs, size=size, wspace=wspace, show=False)

        except KeyError as e:
            e = str(e).strip().replace("'", '')
            catchStr = 'Could not find key leiden in .var_names or .obs.columns.'

            if e == catchStr:
                raise KeyError("Could not find key leiden in .var_names or .obs.columns! Please run leiden() first!")

        # set the leiden colors and regenerate if no colors have been set
        if getLeidenColors(adata) is not None:
            loopForColors = False

            if show:
                plt.show()

        else:
            # clear the none colored graph
            plt.close()

            # default leiden color initialization
            for graph in graphs:
                if graph.lower() == "leiden":  # check if leiden is in colors
                    leidenColorInit(adata=adata)

    return


def clustering(adata):
    pass


# calculate spatial neighbors data
def gr_spatialNeighbors(adata, coordType="generic", spatialKey="spatial", delaunay=False):
    return sq.gr.spatial_neighbors(adata=adata, coord_type=coordType, spatial_key=spatialKey, delaunay=delaunay)


# calculate nhoodEnrichment
def gr_nhoodEnrichment(adata, clusterKey="leiden"):
    return sq.gr.nhood_enrichment(adata=adata, cluster_key=clusterKey)


# plot nhoodEnrichment data: No return value
def pl_nhoodEnrichment(adata, show=False, clusterKey="leiden", method="average", cmap="inferno", vmin=-50, vmax=100,
                       figsize=(5, 5)):
    sq.pl.nhood_enrichment(adata=adata, cluster_key=clusterKey, method=method, cmap=cmap, vmin=vmin,
                           vmax=vmax, figsize=figsize)

    if show:
        plt.show()


# compute the centrality scores
def gr_centrality_scores(adata, cluster_key="leiden"):
    return sq.gr.centrality_scores(adata=adata, cluster_key=cluster_key)


# graph the centrality scores
def pl_centrality_scores(adata, cluster_key="leiden", figsize=None):
    return sq.pl.centrality_scores(adata=adata, cluster_key=cluster_key, figsize=figsize)


# compute the co-occurrence probability
def gr_co_occurrence(adata, cluster_key="leiden"):
    return sq.gr.co_occurrence(adata=adata, cluster_key=cluster_key)


# graph the co-occurrence probability
def pl_co_occurrence(adata, cluster_key="leiden", clusters="12", figsize=None):
    return sq.pl.co_occurrence(adata=adata, cluster_key=cluster_key, clusters=clusters, figsize=figsize)


# compute Ripley's statistics
def gr_ripley(adata, cluster_key="leiden", mode="L"):
    return sq.gr.ripley(adata, cluster_key=cluster_key, mode=mode)


# graph Ripley's statistics
def pl_ripley(adata, cluster_key="leiden", mode="L"):
    return sq.pl.ripley(adata, cluster_key=cluster_key, mode=mode)


def gr_spatialAutocorr(adata, mode="moran", nPerms=100, nJobs=1):
    return sq.gr.spatial_autocorr(adata, mode=mode, n_perms=nPerms, n_jobs=nJobs)


def leiden(adata, resolution=1.0, ignore=False):
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
    return sc.pp.filter_cells(adata, min_counts=minCounts)


# can be used to filter genes that are expressed in too few cells
def filterGenes(adata, minCells=10):
    return sc.pp.filter_genes(adata, min_cells=minCells)


# can be used to scale gene expression. IE Clip values that exceed 10 ('max value') standard deviations
def scale(adata, maxValue=10):
    return sc.pp.scale(adata, max_value=maxValue)


# plots transcript data
# requires running of qc metrics
def plotTranscripts(adata, show=False, figsize=(15, 4)):
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


# endregion


# easier to use version of spatial_scatter than sq default
# takes in adata and a list of colors / filters IE ["leiden", "n_counts"]
#
# must setup adata by running adataSetup before this will work
def spatialScatter(adata, graphs, show=False, libraryID=None, wspace=0.4, size=None, shape=None,
                   groups=None, cmap=None, figsize=None):
    loopForColors = True
    while loopForColors:
        try:
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

            # set the leiden colors and regenerate if no colors have been set
            if getLeidenColors(adata) is not None:
                loopForColors = False

                if show:
                    plt.show()

            else:
                # clear the none colored graph
                plt.close()

                # default leiden color initialization
                for graph in graphs:
                    if graph.lower() == "leiden":  # check if leiden is in colors
                        leidenColorInit(adata=adata)

        except ValueError:
            # default leiden color initialization
            for graph in graphs:
                if graph.lower() == "leiden":  # check if leiden is in colors
                    leidenColorInit(adata=adata)

    return


# default leiden color initialization
def leidenColorInit(adata):
    # check if leiden colors have already been set
    currentLeidenColors = getLeidenColors(adata=adata)
    colors = currentLeidenColors

    colorLength = len(adata.obs["leiden"].value_counts())

    if colors is None:
        # get the default color set
        leidenColors = getColors('leiden_generated_random_3.txt')

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
    if colors is None:
        print("No colors given, generating and saving now!")
        leidenColors = createPalette(length=length, save=True, log=True)

        adata.uns['leiden_colors'] = leidenColors
        return leidenColors

    else:
        adata.uns['leiden_colors'] = colors
        return colors


def getLeidenColors(adata):
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
    plt.show()


# search for a 'colors' directory which contains the following file names
# defaults to the current working dir.
# can pass in a custom file
# defaults to searching in the 'colors' directory
def getColors(color_file):
    # search for the given file
    color_path = str(os.getcwd()).replace('\\', '/') + '/'
    color_path = searchFiles(color_path, 'colors')

    # read the colors from the file and return them
    current_color_path = color_path + '/' + color_file
    with open(current_color_path, "r") as color_file:
        colors = color_file.read().split('\n')

    return colors


# creates a random color palette of given length
# saves this color palette in 'colors' directory
def createPalette(length, save=False, log=True):
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
    path = os.getcwd().replace('\\', '/') + '/tutorial_data/'
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

    # todo implement clustering
