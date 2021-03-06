"""
Functions for plotting locations of given clusters from sequencing data CPseq files.
"""

import matplotlib.pyplot as plt

def plot_species_per_tile(file_name, contains=None, color = 'b'):
    """
    Args:
        file_name - str, name of 1 filtered tiles CPseq file
    """
    x_s = []
    y_s = []

    with open(file_name, 'r') as fil:

        for line in fil:
            linesplit = line.split("\t")
            if (contains in linesplit[1]) or (contains in linesplit[2]) or (contains in linesplit[4]):
                locations = linesplit[0].split(":")
                x = int(locations[5])
                y = int(locations[6])*-1
                x_s.append(x)
                y_s.append(y)
    
    plt.scatter(x_s, y_s, s=1, c = color)


def plot_fiducial_in_tile(CPseq_name, fig_name):
    """
    Args:
        CPseq_name - str, name of 1 filtered tiles CPseq file
        fig_name - str, output
    """
    plt.figure(figsize=(11,11))

    # plot one or multiple species here
    plot_species_per_tile(CPseq_name, 'FID', 'k')

    plt.axis('off')
    plt.savefig(fig_name, dpi=300, bbox_inches='tight')


def plot_fiducial_all_tiles(CPseq_names, fig_names):
    """
    Args:
        CPseq_names - List[str]
        fig_names - List[str]
    """
    for seqfile,figfile in zip(CPseq_names, fig_names):
        plot_fiducial_in_tile(seqfile, figfile)


print('\nPlotting fiducials from input:\n', snakemake.input)
plot_fiducial_all_tiles(snakemake.input, snakemake.output)
