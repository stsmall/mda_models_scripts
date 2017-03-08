import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import colors
import numpy as np


def plot_allele_frequency(dfworm, outfile="./temp/testfig.png"):
    """Plots site frequency spectrum and joint SFS
    Returns
    -------
    SFS : fig
        plot of SFS with each a different color on same plot
    jSFS : fig
        joint SFS between 2 villages
    """
    try:
        fig = sns.distplot(df_allele).get_figure()
    except np.linalg.LinAlgError:
        from IPython import embed
        embed()
    fig.savefig(outfile)

def sel_allele_trace(dfworm, outfile="./temp/testfig.png"):
    """Plots trace of selected alleles
    """
    #dfworm.sel['1S'] > 0
    #dfworm.sel['2S'] > 0
    #dfworm.sel['1F'] > 0
    #np.where(dfworm.sel['2F'] > 0)
    #fig.savefig(outfile)

def plot_coordinates_host(dfHost, height=7, width=7):
    """ Plots host coordinates
    """
    cmap = colors.ListedColormap(['k','g','y','r','b'])
    bounds=[0,1,2,3,4]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    fig, ax = plt.subplots()
    ax.scatter([i[0] for i in dfHost.coordinates],
            [i[1] for i in dfHost.coordinates],
            c = dfHost.village, cmap=cmap)
    ax.set_ylabel('y position')
    ax.set_xlabel('x position')
    fig.savefig('test.png')
    plt.tight_layout()
    return(fig)
