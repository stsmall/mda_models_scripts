import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def plot_allele_frequency(df_allele, outfile="./temp/testfig.png"):
    """Plots allele frequency
    """

    fig = sns.distplot(df_allele).get_figure()
    fig.savefig(outfile)


def plot_coordinates_host(dfHost):
    """ Plots host coordinates
    """
    pass
