"""
Create boxplots using matplotlib
"""

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def boxplot(boxplot_name, title, x_axis, y_axis, lists, x_ticks):

    """Create a boxplot from a list of data and save the graph to an
    output file. To make a plot with multiple boxes, pass a list of lists.

    Parameters
    ___________

    boxplot_name:
        The file name the graph will be saved under. Must use a supported file
        extension, such as .png. Input as a string.

    title: Title of boxplot, displayed at top of plot. A string.

    x_axis: label for x axis. A string.

    y_axis: label for y axis. A string.

    lists: list of ints or floats
        Data to be graphed. May pass a list of lists to plot multiple samples

    x_ticks: The x-axis labels for each individual box. Pass a list of strings,
            corresponding to the index of sample in lists

    Returns
    _________

    Box plot of the input lists saved as out_file_name in the current directory

    """


    out_file = boxplot_name
    figure_title = title
    figure_xaxis = x_axis
    figure_yaxis = y_axis
    data_to_plot = lists
    sample_labels = x_ticks

    width = 10
    height = 3

    fig = plt.figure(figsize=(width, height), dpi=300)
    ax = fig.add_subplot(1, 1, 1)

    ax.boxplot(data_to_plot)
    ax.set_title(figure_title)
    # set custom labels with the names of each list
    ax.set_xticklabels(sample_labels)
    plt.xticks(rotation=90)
    ax.set_xlabel(figure_xaxis)
    ax.set_ylabel(figure_yaxis)

    plt.savefig(out_file, bbox_inches='tight')
