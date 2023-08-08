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
    meta_data_path = 'cell_metadata.csv'
    cell_by_gene_path = 'cell_by_gene.csv'

    meta_data_file = searchFiles(data_path, meta_data_path)
    cell_by_gene_file = searchFiles(data_path, cell_by_gene_path)

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
    printAvailableGraphs(adata)

    print("\ncalculating QC metrics...")
    perUnassigned = qc_metrics(adata)
    print("percent unassigned {}%".format(perUnassigned))

    print("\nrunning a bunch of data calculations...")
    print("layers...")
    adata.layers["counts"] = adata.X.copy()
    print("\nhighly variable genes...")
    highlyVariableGenes(adata)
    print("\nnormalize total...")
    normalizeTotal(adata)
    print("\nlog1p...")
    log1p(adata)
    print("\nprinciple component analysis...")
    pca(adata)
    print("\nneighbors...")
    neighbors(adata)
    print("\ncalculate UMAP...")
    calcUMAP(adata)
    print("\nleiden...")
    leiden(adata)

    print("AFTER:")
    printAvailableGraphs(adata)
    print("finished setup calculations\n")


# region adataSetup functions
def qc_metrics(adata, percentTop=(50, 100)):
    try:
        adata.var_names_make_unique()
        adata.var["mt"] = adata.var_names.str.startswith("mt-")

        sc.pp.calculate_qc_metrics(adata, percent_top=percentTop, inplace=True, qc_vars=["mt"])
        perUnassigned = adata.obsm["blank_genes"].to_numpy().sum() / adata.var["total_counts"].sum() * 100

        return perUnassigned

    except:  # yes, I know it's bad to do this.
        print("unknown error, running without qc_vars")

    sc.pp.calculate_qc_metrics(adata, percent_top=percentTop, inplace=True)
    perUnassigned = adata.obsm["blank_genes"].to_numpy().sum() / adata.var["total_counts"].sum() * 100

    return perUnassigned


def highlyVariableGenes(adata, nTopGenes=4000):
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=nTopGenes)


def normalizeTotal(adata):
    sc.pp.normalize_total(adata)


def log1p(adata):
    sc.pp.log1p(adata)


def pca(adata):
    sc.pp.pca(adata)


def neighbors(adata):
    sc.pp.neighbors(adata)


def calcUMAP(adata):
    sc.tl.umap(adata)


def leiden(adata, ignore=False):
    try:
        sc.tl.leiden(adata)

    except KeyError as e:
        e = str(e).strip().replace("'", '')
        catchStr = 'No "neighbors" in .uns'

        if e == catchStr and not ignore:
            warnNeighborsStr = "No neighbors found, running 'neighbors()' default now! " \
                               "Run neighbors() with your args first if this is not what you want!"
            warnings.warn("Warning................................\n{}".format(warnNeighborsStr))

            neighbors(adata)
            leiden(adata, ignore=True)  # call leiden again with ignore, so if we get this error we will raise
            return

        else:
            raise

# endregion


# easier to use version of spatial_scatter than sq default
# takes in adata and a list of colors / filters IE ["leiden", "n_counts"]
#
# must setup adata by running adataSetup before this will work
def spatialScatter(adata, graphs, show=True, colors=None):
    # default leiden color initialization
    for color in graphs:
        if color.lower() == "leiden":  # check if leiden is in colors
            try:  # see if the leiden_colors are already set
                defaultColors = adata.uns['leiden_colors']  # does nothing, as this is an error test
                break

            except KeyError:
                pass

            # if not, set them
            colorLength = len(testAdata.obs["leiden"].value_counts())

            if colors is None:
                leidenColors = getColors()[0]

            else:
                leidenColors = colors

            while len(leidenColors) > colorLength:
                leidenColors.pop(-1)

            adata.uns['leiden_colors'] = leidenColors

    sq.pl.spatial_scatter(
        adata,
        shape=None,
        color=graphs,
    )

    if show:
        plt.show()


# search for a 'colors' directory which contains the following file names
# defaults to the current working dir.
# can pass in a custom file
def getColors(color_file=None):
    # note: import color files that I have created
    color_path = str(os.getcwd()).replace('\\', '/') + '/'

    color_path = searchFiles(color_path, 'colors')

    if color_file is not None:
        color_path = color_file

    color_file_1 = 'leiden_color_set_1_gradient.csv'
    color_file_2 = 'leiden_color_set_1_random.csv'
    color_file_3 = 'leiden_color_set_2_random.csv'
    color_file_4 = 'leiden_color_set_3_random.csv'

    color_files = [color_file_1, color_file_2, color_file_3, color_file_4]
    colorData = []

    for file in color_files:
        current_color_path = color_path + '/' + file
        with open(current_color_path, "r") as color_file:
            colors = color_file.read().split('\n')

        colorData.append(colors)

    return colorData


# searches a provided directory and subdirectories for a file or dir
def searchFiles(data_path, fileToFind):
    data_path = copy.deepcopy(data_path)
    fileToFind = copy.deepcopy(fileToFind)
    dirList = os.listdir(data_path)

    dirs = []
    for file in dirList:
        file = str(file)
        # print(file)

        if file == str(fileToFind):
            return data_path + file

        elif file.find('.') == -1:  # this is a directory
            dirs.append(file)

    # print(dirs)

    for dir in dirs:
        # print(dir)
        new_path = data_path + dir + '/'

        returnVal = searchFiles(new_path, fileToFind)

        if returnVal is not None:
            return returnVal


# these are values that can be passed into spatialScatter() colors for graphing
def printAvailableGraphs(adata):
    var_list = list(adata.var_names)
    obs_columns_list = list(adata.obs.columns)

    print("AVAILABLE GRAPHS:\n {}\n {}".format(var_list, obs_columns_list))


# tests for EasySQ.py
if __name__ == "__main__":
    """
    # note: createCode(codeLength) tests
    testCode = createCode(8)
    print(testCode)
    # """

    # """
    # note: readVizgen(data_path) tests
    path = 'F:/sunlabmerfishdata/QSFL01222023/'
    transformation_path = 'C:/Users/scrouse2/Crouse_Work_Files/PyCharm_Projects/EasySQ/'
    testAdata = readVizgen(path, transform_path=transformation_path)
    print(testAdata)

    # note: spatialScatter() tests and supporting funcs like getColors() and adataSetup()
    print("spatialScatter() and analysis tests")
    printAvailableGraphs(testAdata)  # graphs available before setup
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
    # """

    # todo implement clustering
