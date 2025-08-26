import matplotlib.pyplot as plt

from pywellgeo.well_tree.well_tree_tno import *

import os

def plot_wells(wells, filename_noext="well"):
    """
     plot the well trees in 3D and 2D projections
    :param wells: list of WellTreeTNO objects or a single WellTreeTNO object
    :param filename_noext: name of the file to save the plot
    """

    if (filename_noext!= None):
        directory = os.path.dirname(filename_noext)
        if not os.path.exists(directory):
            os.makedirs(directory)
        tofile = filename_noext + "_welltree_3D.png"
    else:
        tofile = None

    for i, well in enumerate(wells):
        if not isinstance(well, WellTreeTNO):
            raise TypeError(f"Expected WellTreeTNO object, got {type(well)} at index {i}.")
        if (i == 0):
            fig, ax = well.plotTree3D(doplot=(i == len(wells) - 1), tofile= tofile)
            ax.set_aspect('equal')
        else:
            well.plotTree3D(fig=fig, ax=ax, doplot=(i == len(wells)-1), tofile=tofile)
            ax.set_aspect('equal')
    plt.close()

    axes = [0, 1, 2]
    axes_s = ["xy", "xz", "yz"]
    for axi in axes:
        for i, well in enumerate(wells):
            if (i == 0):
                fig, ax = well.plotTree(axis=axi)
                ax.set_aspect('equal')
            else:
                well.plotTree(fig=fig, ax=ax, axis=axi)
                ax.set_aspect('equal')
        if (filename_noext!=None):
            plt.savefig(filename_noext + "_welltree_%s.png" % axes_s[axi])
        else:
            plt.show()
        plt.close()
