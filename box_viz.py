import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt




def boxplot(file, groups, group_col_name, gene_name, boxplot_name):


    data = []
    for line in file:
        data = data.append(line)

    num_ticks = [i for i in range(1, len(groups)+2)]

    fig = plt.figure(figsize=(10, 3), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    ax.boxplot(data)
    plt.xticks(num_ticks, groups, rotation='vertical')
    plt.title(gene_name)
    plt.xlabel(group_col_name)
    plt.ylabel('gene read #')
    plt.savefig(boxplot_name, bbox_inches='tight')
