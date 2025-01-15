import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from matplotlib.patches import ConnectionPatch

def Plot_graph(g, file, name, centro):  # this function plot the graph
    vertices = g.vertices
    fig = plt.figure(figsize=(20, 25))
    plt.title(name)
    # fig.suptitle(args.name, fontsize=48)
    rows = 6
    column = 4
    grid = plt.GridSpec(rows, column, wspace=.25, hspace=.25)
    axes = []
    for i in range(24):
        exec(f"plt.subplot(grid{[i]})")
        axes.append(plt.gca())
        if centro != None:
            plt.axvline(min(centro['chr' + str(i + 1)]), color='green')
            plt.axvline(max(centro['chr' + str(i + 1)]), color='green')
        max_m = 0
        for v in vertices:
            if int(v.chromosome) == i + 1:
                plt.plot(int(v.pos), int(v.cn), marker="o", markersize=2, color="green", alpha=0.7)
                if v.cn > max_m:
                    max_m = v.cn
        for e in g.edges:
            # print(e[0])
            node1 = g.return_node(e[0])
            node2 = g.return_node(e[1])
            if int(node1.chromosome) == i + 1 and int(node2.chromosome) == i + 1 and e[3] == 'S':
                plt.plot([node1.pos, node2.pos], [node1.cn, node2.cn], markersize=0.5, color="red", alpha=0.3)
            elif int(node1.chromosome) == i + 1 and int(node2.chromosome) == i + 1 and e[3] == 'R':
                plt.plot([node1.pos, node2.pos], [node1.cn, node2.cn], markersize=0.5, color="blue", alpha=0.3)
            elif int(node1.chromosome) == i + 1 and int(node2.chromosome) == i + 1 and e[
                3] == 'SV' and node1.pos != node2.pos:
                x = [node1.pos, (node1.pos + node2.pos) / 2, node2.pos]
                y = [node1.cn, max(node1.cn, node2.cn) + 1, node2.cn]
                if node2.pos < node1.pos:
                    x = [node2.pos, (node1.pos + node2.pos) / 2, node1.pos]
                    y = [node2.cn, max(node1.cn, node2.cn) + 1, node1.cn]
                # print(x)
                x2 = np.linspace(x[0], x[-1], 100)
                y2 = interpolate.pchip_interpolate(x, y, x2)
                plt.plot(x2, y2, markersize=0.5, color="black", alpha=0.3)
        plt.title('chr' + str(i + 1))
        plt.ylim([-0.5, max_m + 3])
    for e in g.edges:
        node1 = g.return_node(e[0])
        node2 = g.return_node(e[1])
        if int(node1.chromosome) != int(node2.chromosome) and e[3] == 'SV':
            xy1 = (node1.pos, node1.cn)
            xy2 = (node2.pos, node2.cn)
            con1 = ConnectionPatch(xyA=xy2, xyB=xy1, coordsA="data", coordsB="data",
                                   axesA=axes[int(node2.chromosome) - 1], axesB=axes[int(node1.chromosome) - 1],
                                   color="black", alpha=0.6)
            axes[int(node2.chromosome) - 1].add_artist(con1)
    i = -0.5
    j = -0.5
    prev = 0
    for e in g.edges:
        node1 = g.return_node(e[0])
        node2 = g.return_node(e[1])
        if int(node1.chromosome) != int(node2.chromosome) and e[3] == 'SV':
            if int(node1.chromosome) != prev:
                prev = int(node1.chromosome)
                i = 0
            else:
                i = i + 0.5
            exec(f"plt.subplot(grid{[int(node1.chromosome) - 1]})")
            xy1 = (node1.pos, node1.cn)
            xy2 = (node2.pos, node2.cn)
            if node1.type == 'H' and node2.type == 'H':
                plt.annotate('(_,+)', xy=(node1.pos, node1.cn + 0.5 + i))
            elif node1.type == 'H' and node2.type == 'T':
                plt.annotate('(_,_)', xy=(node1.pos, node1.cn + 0.5 + i))
            elif node1.type == 'T' and node2.type == 'T':
                plt.annotate('(+,_)', xy=(node1.pos, node1.cn + 0.5 + i))
            elif node1.type == 'T' and node2.type == 'H':
                plt.annotate('(+,+)', xy=(node1.pos, node1.cn + 0.5 + i))
            exec(f"plt.subplot(grid{[int(node2.chromosome) - 1]})")
            j += 0.5
            if node1.type == 'H' and node2.type == 'H':
                plt.annotate('(_,+)', xy=(node2.pos, node2.cn + 0.5 + j))
            elif node1.type == 'H' and node2.type == 'T':
                plt.annotate('(_,_)', xy=(node2.pos, node2.cn + 0.5 + j))
            elif node1.type == 'T' and node2.type == 'T':
                plt.annotate('(+,_)', xy=(node2.pos, node2.cn + 0.5 + j))
            elif node1.type == 'T' and node2.type == 'H':
                plt.annotate('(+,+)', xy=(node2.pos, node2.cn + 0.5 + j))
    plt.savefig(file, dpi=200)