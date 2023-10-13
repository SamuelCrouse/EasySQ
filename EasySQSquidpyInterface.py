import scanpy as sc
import squidpy as sq
import EasySQ as esq

# other needed libraries
import time
import os


if __name__ == "__main__":
    t0 = time.time()  # time the execution of the program

    path = os.getcwd() + '/tutorial_data_1/'
    esqAn = esq.Analysis(data_path=path)
    esqAn.print()

    esqAn.qcMetrics()

    esqAn.availableGraphs()
    esqAn.spatialScatter(graphs=["Oxgr1"])
    # sq.pl.spatial_scatter(esqAn.getAdata(), color=["Oxgr1"])
