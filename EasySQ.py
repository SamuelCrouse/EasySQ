# Sam Crouse
# scrouse2@uwyo.edu

# Class which implements the tools built in EasySQTools.py
# Allows creation of a squidpy analysis object with built in tools

# python imports


# my imports
import EasySQTools as tools


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

    def spatialScatter(self, graphs, show=True, colors=None):
        return tools.spatialScatter(adata=self.getAdata(), graphs=graphs, show=show, colors=colors)

    # endregion

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
    esqAnalysis.spatialScatter('Th')

