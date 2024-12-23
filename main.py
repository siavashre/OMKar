from collections import defaultdict
import os
import sys
import shutil
from parsers import *
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import pulp as p
import argparse
import itertools
from copy import deepcopy
import math
from matplotlib import rcParams
from utill import *
from bionano_metadata import *
import pandas as pd
import csv
from KarReporter.KarUtils.read_OMKar_output import *
from pathlib import Path
rcParams['pdf.fonttype'] = 42


class Vertices:  # Class of vertices in out graph
    def __init__(self):
        self.type = ''  # Type of vertices it can be Head or Tail
        self.id = 0  # Id of vertices
        self.chromosome = 0  # which chromosome this vertices belong chrX = 23 and chrY = 24
        self.pos = 0  # in terms of bp the position of this vertices
        self.cn = 0  # integer CN of segment this vertices represent
        self.edges = []  # its like adjacenty matrix it contains all adjacent vertices

    def append_edges(self, i):  # add new vertices to adjancty list
        self.edges.append(i)
        self.edges = list(set(self.edges))

    def remove_edges(self, i):  # remove vertices from adjacenty list
        if i in self.edges:
            self.edges.remove(i)
        else:
            print('Not here baba')

    def print_v(self):  # print vertices info
        print(self.id, self.chromosome, self.pos, self.cn, self.type, self.edges)


class Graph:  # Class of graph
    def __init__(self):
        self.vertices = []  # list of Object of vertices
        self.edges = []  # list of edges. each edges is tuple of (u,v,M,type) u and v are vertices connected to each other M represent edge multiplicities and type can be "S", "SV", "R", "D"
        self.terminal_vertices_ids = []  # subset of vertices that begin/end a chromosome

    def print_node(self):  # print all nodes
        for v in self.vertices:
            print(v.id, v.chromosome, v.pos, v.cn, v.type, v.edges)

    def output_node(self, out_file):
        out_str = ""
        for v in self.vertices:
            out_str += "{}\t{}\t{}\t{}\t{}\t{}\n".format(v.id, v.chromosome, v.pos, v.cn, v.type, v.edges)
        with open(out_file, 'w') as fp_write:
            fp_write.write(out_str)

    def output_edges(self, out_file):
        out_str = ""
        for e in self.edges:
            out_str += str(e) + '\n'
        with open(out_file, 'w') as fp_write:
            fp_write.write(out_str)

    def return_node(self, id):  # return node by id
        for v in self.vertices:
            if id == v.id:
                return v
        return None

    def return_edges(self, a, b):  # return all type of edge between two nodes. it return a list of edges or None
        ans = []
        for e in self.edges:
            if (a == e[0] and b == e[1]) or (a == e[1] and b == e[0]):
                ans.append(e)
        if len(ans) > 0:
            return ans
        return None

    def add_dummy_edges(self, u, v, cnn):  # add dummy edges with type "D" between two nodes if already edge exist between them increase the CN
        e = self.return_edges(u, v)
        if max(u, v) % 2 == 1 and abs(u - v) == 1:
            if len(e) == 1:
                # if e[0][2] == 2:
                #      e = e[0]
                #      cn = e[2] - 1
                #      self.edges.remove(e)
                #      self.edges.append((e[0], e[1], cn, e[3]))
                # else:
                self.edges.append((u, v, cnn, 'D'))
                self.return_node(u).append_edges(v)
                self.return_node(v).append_edges(u)
        elif e == None:
            self.edges.append((u, v, cnn, 'D'))
            self.return_node(u).append_edges(v)
            self.return_node(v).append_edges(u)
        else:
            if len(e) == 1:
                e = e[0]
                cn = e[2] + cnn
                self.edges.remove(e)
                self.edges.append((e[0], e[1], cn, e[3]))
            else:
                e = e[0]
                cn = e[2] + cnn
                self.edges.remove(e)
                self.edges.append((e[0], e[1], cn, e[3]))

    def update_edges(self, a, b, count, e_type):  # update an edge between node a nad b with type = e_type to the new number
        e_list = self.return_edges(a, b)
        for e in e_list:
            if e[3] == e_type:
                self.edges.remove(e)
                if count > 0:
                    self.edges.append((e[0], e[1], count, e[3]))

    def label_terminal_vertices(self):
        ## mark the first and the last vertices on this graph, for each chromosome
        chr_min_max = {}
        for v in self.vertices:
            if v.chromosome not in chr_min_max:
                chr_min_max[v.chromosome] = {'min': v.pos, 'max': v.pos}
            else:
                if v.pos < chr_min_max[v.chromosome]['min']:
                    chr_min_max[v.chromosome]['min'] = v.pos
                if v.pos > chr_min_max[v.chromosome]['max']:
                    chr_min_max[v.chromosome]['max'] = v.pos
        for v in self.vertices:
            if v.pos == chr_min_max[v.chromosome]['min'] or v.pos == chr_min_max[v.chromosome]['max']:
                self.terminal_vertices_ids.append(v.id)

    # def previous_vertex_of_same_chrom(self, current_v):
    #     current_v_idx = self.vertices.index(current_v)
    #     if current_v_idx == 0:
    #         return None
    #     previous_v = self.vertices[current_v_idx - 1]
    #     if previous_v.chromosome == v.chromosome:
    #         return previous_v
    #     else:
    #         return None
    #
    # def next_vertex_of_same_chrom(self, current_v):
    #     current_v_idx = self.vertices.index(current_v)
    #     if current_v_idx == len(self.vertices) - 1:
    #         return None
    #     next_v = self.vertices[current_v_idx + 1]
    #     if previous_v.chromosome == v.chromosome:
    #         return previous_v
    #     else:
    #         return None


def find_bp_in_segment(chromosome, point,
                       segments):  # this function get chromosome and position as input(from sv call) and add this point to segment bp list. It helps us to then spliting segments to new one.
    for i in segments:
        if str(i.chromosome) == str(chromosome):
            if i.start <= point <= i.end:
                i.bp.append(point)


def merge_list(l):  # this function merge all breakpoints in a segment within a window of 50Kbp
    l = list(set(l))
    l = sorted(l)
    if len(l) == 2:
        return l
    WINDOW = 50000
    ans = []
    prev = None
    group = []
    for item in l[1:-1]:
        if prev is None or abs(item - prev) <= WINDOW:
            group.append(item)
        else:
            ans.append(np.mean(group))
            group = [item]
        prev = item
    if len(group) > 0:
        ans.append(np.mean(group))
    return [l[0]] + ans + [l[-1]]


def rev_dir(a):
    if a == 'H':
        return 'T'
    return 'H'





def find_start_end(prev_point, start, label_list):
    ans = []
    for i in label_list:
        if prev_point <= i < start:
            ans.append(i)
    if len(ans) < 2:
        return 0, 0
    else:
        return min(min(ans)+1, prev_point+1), start-1


def next_prev_label(chromosome, pos):  # not used can be deleted
    for k in rcop.keys():
        label_list = list(rcop[k].keys())
        if k == chromosome:
            for index, i in enumerate(label_list):
                if i == pos:
                    return label_list[index - 1], label_list[(index + 1) % (len(label_list))]


def detect_overlap_map(chromosome, pos, xmap):  # this function detect if there is a map(alignment that overlap this region)
    # maybe it is good idea instead of 1bp assume a window here.
    min_x = float(pos) - 20000
    max_x = float(pos) + 20000
    min_seen = 1000000000
    max_seen = 0
    for xmapid in xmap.keys():
        x = xmap[xmapid]
        if int(chromosome) == int(x['RefContigID']):
            if min_seen > x['RefStartPos']:
                min_seen = x['RefStartPos']
            if max_seen < x['RefEndPos']:
                max_seen = x['RefEndPos']
            if x['RefStartPos'] <= min_x <= x['RefEndPos'] and x['RefStartPos'] <= max_x <= x['RefEndPos']:
                return True
    if min_seen <= float(pos):
        return True
    if max_seen >= float(pos):
        return True
    return False

def Plot_graph(g, file, name):  # this function plot the graph
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



def dfs(i, temp, g, visited):  # run dfs alg on graph g on vertices i. visited is a list of seen vertices before visiting this node.
    visited.append(i)
    temp.append(i)
    for v in g.return_node(i).edges:
        if v not in visited:
            temp = dfs(v, temp, g, visited)
    return temp

def check_contains_masked_region(l_m, start, end):
    contains = False
    border_start = False
    border_end = False
    for i in l_m:
        if i[0] - 0 <= start < i[1] + 0:
            border_start = True
            contains = True
        if i[0] - 0 <= end <= i[1] + 0:
            border_end = True
            contains = True
        if start <= i[0] and end >= i[1]:
            contains = True
    return contains, border_start, border_end

def find_connected_components(g):  # find connected components in a graph g
    cc = []
    visited = []
    for i in range(len(g.vertices)):
        if i not in visited:
            temp = []
            a = dfs(i, temp, g, visited)
            if len(a) > 1:
                cc.append(a)
    return cc




def calculate_seg_length(e, g):  # calculate segment length of segment edge e
    l = abs(g.return_node(e[0]).pos - g.return_node(e[1]).pos)
    return l, 1 + math.ceil(l / 5000000)


def estimating_edge_multiplicities_in_CC(component, g, xmap):
    Lp_prob = p.LpProblem('Problem', p.LpMinimize)  # this line initiate the ILP problem
    objective = 0  # variable for some of our objective
    sv_sum = 0  # sum all variables for SV edges
    cn_tune = 0  # sum all variables for changing CN (tuning CN)
    odd_vertices_pen = 0
    component_edges = return_all_edges_in_cc(component, g)
    print('Siavash', component_edges)
    for i in range(len(component_edges)):
        e = component_edges[i]
        if e[3] != 'S':  # in it is not Segment edge
            component_edges[i] = [e, p.LpVariable('X' + str(i), lowBound=0, cat=p.LpInteger)]  # create an ILP variable Xi for >= 0
            print('X' + str(i), e)
            if e[3] == 'SV':
                sign_func = p.LpVariable('T' + str(i), cat=p.LpBinary)
                Lp_prob += component_edges[i][1] <= 10 * sign_func
                Lp_prob += component_edges[i][1] >= sign_func
                if g.return_node(e[0]).chromosome != g.return_node(e[1]).chromosome:  # I want to give more weight to the translocation SVs
                    Lp_prob += component_edges[i][1] >= 0  # sum in the LP_probe
                    # sv_sum += 8 * component_edges[i][1] #updating sv_sum with weight of 8 of this variable
                    sv_sum += 8 * sign_func  # updating sv_sum with weight of 8 of this variable
                else:
                    Lp_prob += component_edges[i][1] >= 0
                    # sv_sum += component_edges[i][1]
                    sv_sum += sign_func
            elif e[
                3] == 'R':  # if it is reference edge we want to have another constraint if a map overlap it it should be traverse at least one so updating the constraints
                node = g.return_node(e[0])
                if detect_overlap_map(node.chromosome, node.pos, xmap):
                    Lp_prob += component_edges[i][1] >= 1
        else:  # if it is Segment edge
            component_edges[i] = [e, p.LpVariable('Y' + str(i), cat=p.LpInteger),
                                  p.LpVariable('Z' + str(i), cat=p.LpInteger)]  # variable Y_i is for tuning CN of each segment if we need to change the CN. Z = absolute value of Y
            print('Y' + str(i), e)
            # For having abselute value in ILP for we need to define variable Z as well.
            Lp_prob += component_edges[i][2] >= -component_edges[i][1]  # this is used for defining abselute value
            Lp_prob += component_edges[i][2] >= component_edges[i][1]
            # the following if is a condition for setting the limit on how much the CN can be changed for segment with length less than 1Mbp it is one
            if calculate_seg_length(e, g)[0] <= 1000000:
                Lp_prob += component_edges[i][2] <= 1
            else:  # for rest it is 25% of their CN
                Lp_prob += component_edges[i][2] <= math.ceil(e[2] / 4)
    for v in component:
        v_edges = []  # all edges conntigin to node v
        for e in component_edges:
            if e[0][0] == v or e[0][1] == v:
                v_edges.append(e)
        if len(v_edges) > 1:
            cond = 0  # this is a CN balance condition
            for e in v_edges:
                if e[0][3] == 'SV' and e[0][0] == e[0][1]:
                    cond = cond + 2 * e[1]  # it there is a SV edge from node v to v
                    objective = objective - 2 * e[1]
                elif e[0][3] != 'S':
                    cond = cond + e[1]  # sum all non Segment edge
                    objective = objective - e[1]
            for e in v_edges:
                if e[0][3] == 'S':
                    cond = cond + e[1]  # summing tuning variable
                    Lp_prob += cond <= e[0][2]  # this is important line wich we apply copy number balance condition
                    if calculate_seg_length(e[0], g)[0] > 50000:  # for segment less than 50 Kbp no penalty applied for tuning segment CN
                        if e[0][2] != 2:
                            cn_tune += e[2] * calculate_seg_length(e[0], g)[1] * 3  # more coeeficient for something has called already
                        else:
                            cn_tune += e[2] * calculate_seg_length(e[0], g)[1]
                    aa = p.LpVariable('A' + str(v), lowBound=0, cat=p.LpInteger)
                    print('A' + str(v), aa)
                    bb = p.LpVariable('B' + str(v), lowBound=0, cat=p.LpBinary) #This variable shows if sum dgree this vertex is odd or no
                    print('B' + str(v), bb)
                    Lp_prob += cond + e[0][2] - 2 * aa == bb
                    objective = objective + e[0][2] - e[1]  # updating the Objective Function
                    odd_vertices_pen = odd_vertices_pen + bb
                    # objective = objective + bb
        else:
            objective = objective + v_edges[0][0][2] - v_edges[0][1]
            aa = p.LpVariable('A' + str(v), lowBound=0, cat=p.LpInteger)
            print('A' + str(v), aa)
            bb = p.LpVariable('B' + str(v), lowBound=0, cat=p.LpBinary) #This variable shows if sum dgree this vertex is odd or no
            print('B' + str(v), bb)
            Lp_prob += v_edges[0][0][2] - v_edges[0][1] - 2 * aa == bb
            # objective = objective +  bb
            odd_vertices_pen = odd_vertices_pen + bb
            if calculate_seg_length(v_edges[0][0], g)[0] > 50000:  # for segment less than 50 Kbp no penalty applied for tuning segment CN
                if v_edges[0][2] != 2:
                    cn_tune += v_edges[0][2] * calculate_seg_length(v_edges[0][0], g)[1] * 3
                else:
                    cn_tune += v_edges[0][2] * calculate_seg_length(v_edges[0][0], g)[1]
    # Just for debug
    print('obj', objective)
    print('Sv_sum', sv_sum)
    print('CN_tune', cn_tune)
    # objective = 10 * objective  - 9 * sv_sum + 15 * cn_tune
    objective = 10 * objective - 9 * sv_sum + 10 * cn_tune
    objective =  objective -  sv_sum +  cn_tune + odd_vertices_pen
    print('obj', objective)
    Lp_prob += objective
    print(Lp_prob)
    status = Lp_prob.solve()
    print(p.LpStatus[status])  # The solution status
    print(Lp_prob.variables())
    for i in Lp_prob.variables():
        # Set edges multiplicity based on variiable values
        if i.name != '__dummy' and i.name.startswith('X'):
            index = int(i.name[1:])
            component_edges[index][0] = list(component_edges[index][0])
            component_edges[index][0][2] = int(i.varValue)
            component_edges[index] = tuple(component_edges[index][0])
            print(component_edges[index][0], component_edges[index][1], i.varValue, component_edges[index][3])
            g.update_edges(component_edges[index][0], component_edges[index][1], int(i.varValue), component_edges[index][3])
        if i.name != '__dummy' and i.name.startswith('Y'):
            index = int(i.name[1:])
            component_edges[index][0] = list(component_edges[index][0])
            valve = int(component_edges[index][0][2]) - int(i.varValue)
            component_edges[index][0][2] = int(component_edges[index][0][2]) - int(i.varValue)
            component_edges[index] = tuple(component_edges[index][0])
            g.update_edges(component_edges[index][0], component_edges[index][1], valve, component_edges[index][3])
    for v in Lp_prob.variables():
        print(v.name, "=", v.varValue)
    for i in range(len(component_edges)):
        if component_edges[i][1] == None:
            component_edges[i] = component_edges[i][0]
    print(p.value(Lp_prob.objective), 'Objective Answer')
    ans = []
    for i in component_edges:
        if i[2] > 0:
            ans.append(i)
        elif i[2] == 0:  # if estimate of an edge is 0, update the graph  if two different type of edge connect same node it need to be check
            # CHECK CHECK CHECK
            if g.return_edges(i[0], i[1]) == None:
                g.return_node(i[0]).remove_edges(i[1])
                if i[0] != i[1]:
                    g.return_node(i[1]).remove_edges(i[0])
    return ans


def remove_edge(g2, e):  # this function remove edge e from graph g2 and if edge e is not in g print is not in graph
    g = deepcopy(g2)
    component_edges = g.edges
    for i in range(len(component_edges)):
        if e == component_edges[i]:
            if component_edges[i][2] == 1:
                component_edges = component_edges[:i] + component_edges[i + 1:]
                if len(g.return_edges(e[0], e[1])) == 1:  # if only one type of edge connect them to each other.
                    g.return_node(e[0]).remove_edges(e[1])
                    if e[0] != e[1]:
                        g.return_node(e[1]).remove_edges(e[0])
            else:
                component_edges[i] = list(component_edges[i])
                component_edges[i][2] -= 1
                component_edges[i] = tuple(component_edges[i])
            g.edges = component_edges
            return g
    print('Edge is not in graph', e)


def dfs_count(g, v, visited):  # count number of available visired edge from v in graph g
    # for more info look at here https://www.geeksforgeeks.org/fleurys-algorithm-for-printing-eulerian-path/
    count = 1
    visited.append(v)
    for i in g.return_node(v).edges:
        if i not in visited:
            count += dfs_count(g, i, visited)
    return count


def isValidNextEdge(g, u, v, e):  # check this edge is bridge or not. Can we remove the edge?
    if len(g.return_node(u).edges) == 1 and g.return_node(u).edges[0] == v:  # if there is no other way just traverse it
        return True
    else:
        # check if it is bridge or not
        visited = []
        count1 = dfs_count(g, u, visited)
        g2 = remove_edge(g, e)
        visited = []
        count2 = dfs_count(g2, u, visited)
        return False if count1 > count2 else True


def check_traverse_segment(prev, u, next):  # this function check that after a any SV or R edge we traverse Segment edge #this one prevent none meaningfull path
    if prev // 2 == u // 2:
        return True
    elif u // 2 == next // 2:
        return True
    return False

def check_traverse_two_consecutive_segment(prev, u , next,g, chrom):
    if prev == next and g.return_edges(u, next) != None:
        if len( g.return_edges(u, next)) ==1 :
            if g.return_edges(u, next)[0][3]=='S':
                if g.return_node(max(u , next) + 1) != None  and g.return_node(max(min(u , next) - 1, 0)) != None:
                    return False
    return  True


def printEulerUtil(g, u, prev, chrom):  # find Eulerian path or circuts in graph g starting with node u and previouse seen node is prev.
    valid = []
    next_find = False
    next_node = -1
    back_to_chrom = False
    back_to_chrom_node = -1
    for v in g.return_node(u).edges:
        e_list = g.return_edges(u, v)  # all edges between u and v
        for e in e_list:
            if isValidNextEdge(g, u, v, e) and check_traverse_segment(prev, u, v):# and check_traverse_two_consecutive_segment(prev, u, v, g, chrom) :  # check traversing edge between u and v is valid and meaningfull
                valid.append((v, e))
                if abs(v - u) == 1 and v != prev:  # We want to force if next node is availanle traverse it at first between all options.
                    next_find = True
                    next_node = v
                elif g.return_node(u) != g.return_node(v) and g.return_node(v) == chrom:
                    back_to_chrom = True
                    back_to_chrom_node = v
    valid = sorted(valid, key=lambda tup: (tup[0], tup[1][3]))  # this is all valid options for traverse
    # print('haa', u, prev, valid,next_find, next_node)
    if len(valid) == 1:  # if length is equal to one Just traverse it.
        g = remove_edge(g, valid[0][1])
        return [u] + printEulerUtil(g, valid[0][0], u, chrom)
    elif len(valid) > 1:
        # if back_to_chrom:
        #     for i in valid:
        #         if i[0] == back_to_chrom_node:
        #             g = remove_edge(g, i[1])
        #             return [u] + printEulerUtil(g, i[0], u,chrom)
        if next_find:  # if next is available just traverse it
            for i in valid:
                if i[0] == next_node:
                    g = remove_edge(g, i[1])
                    return [u] + printEulerUtil(g, i[0], u, chrom)
        find = False
        for i in valid:
            if i[0] != prev and i[0] != u:  # first if next node is not previouse and current one jump to it
                g = remove_edge(g, i[1])
                find = True
                # print('reza')
                return [u] + printEulerUtil(g, i[0], u, chrom)
        if not find:
            for i in valid:
                if i[0] != prev:
                    g = remove_edge(g, i[1])
                    # print('mamad')
                    return [u] + printEulerUtil(g, i[0], u, chrom)
            g = remove_edge(g, valid[0][1])
            return [u] + printEulerUtil(g, valid[0][0], u, chrom)
    return [u]





def detect_segment_odd_degree(component, component_edges):  # detect vertices with odd degree
    ans = []
    d = {}
    for c in component:
        d[c] = 0
    for e in component_edges:
        d[e[0]] += e[2]
        d[e[1]] += e[2]
    for i in d.keys():
        if d[i] % 2 != 0:
            ans.append(i)
    return ans

def detect_residue_dgree(component, component_edges, terminal_v_ids):
    ans = []
    d = {}
    for c in component:
        d[c] = 0
    for e in component_edges:
        if e[3]=='S':
            d[e[0]] += e[2]
            d[e[1]] += e[2]
        else:
            d[e[0]] -= e[2]
            d[e[1]] -= e[2]
    for i in d.keys():
        if d[i] > 0:
            if i not in terminal_v_ids:
                ans.append((i, d[i]))
    return ans
def scoring_paths(path_list, segment_vertices, g, centro):
    best_score = 99999
    best_path = ''
    for p in path_list:
        ans = []
        temp = [p[0]]
        score = 0
        for i in range(1, len(p) - 1):
            if p[i] in segment_vertices and p[i - 1] == p[i + 1]:
                temp.append(p[i])
                score = score + 10 * abs(check_non_centromeric_path(temp, g, centro) - 1)
                ans.append(temp)
                temp = [p[i]]
            else:
                temp.append(p[i])
        temp.append(p[-1])
        score = score + 10 * abs(check_non_centromeric_path(temp, g, centro) - 1)
        ans.append(temp)
        print('paths_score', score, ans)
        if score < best_score:
            best_score = score
            best_path = p
    return best_path


def printEulerTour(component, component_edges, g, output_file):  # Find Eulerian path/circuts in connected components in graph g
    def vertices_distance(v1, v2):
        ## assumes checkes of intra-chr is done prior
        return abs(g.vertices[v1].pos - g.vertices[v2].pos)

    def group_residual_verticeds(rv):
        ## debugging use
        ## group RVs by Chr and HT-type, hopefully within each Chr, #H == #T
        tally = {}
        for v in rv:
            v_idx = v[0]
            v_count = v[1]
            entry = (g.vertices[v_idx].chromosome, g.vertices[v_idx].type)
            if entry in tally:
                tally[entry] += v_count
            else:
                tally[entry] = v_count
        ## check if HT in each Chr pairs up
        paired_up = True
        for entry, count in tally.items():
            other = 'T' if entry[1] == 'H' else 'H'
            if (entry[0], other) not in tally or tally[((entry[0], other))] != count:
                paired_up = False
                break
        return paired_up, tally

    def segment_edge_is_amplified(v1, v2):
        for e in component_edges:
            if (e[0] == v1 and e[1] == v2) or (e[1] == v1 and e[0] == v2):
                if e[3] == 'S':
                    if e[2] > 2:
                        return True
        return False

    def self_edge_multiplicity(v1):
        for e in component_edges:
            if e[0] == v1 and e[1] == v1:
                if e[3] == 'SV':
                    if e[2] > 0:
                        return e[2]

    g2 = Graph()  # create new graph g2
    g2.edges = component_edges
    for v in component:
        g2.vertices.append(g.return_node(v))
    odd_vertices = []
    print('Component', component)
    print('Component edges', component_edges)

    segment_vertices = detect_segment_vertices(component, component_edges)

    ## ZJ: heuristically connecting residual vertices with closest node (must be intra-chr, must be H-T or T-H, as they represent loss or focal amp.)
    g2.label_terminal_vertices()
    residue_vertices = detect_residue_dgree(component, component_edges, g2.terminal_vertices_ids)
    print('RESIDUE vertices', residue_vertices)
    rv_paired, rv_tally = group_residual_verticeds(residue_vertices)
    if not rv_paired:
        print(f'RV not paired: {rv_tally}')
    else:
        print(f'RV paired: {rv_tally}')

    ## resolve potential self-edge first
    ## only add self-edge if: CN-amp is connecting to another self-edge
    residue_vertices = detect_residue_dgree(component, component_edges, g2.terminal_vertices_ids)
    print('RESIDUE vertices (before self-edge): ', residue_vertices)
    # H node means the other self-edge must be after, T node means its before
    for i1, i1_residual in residue_vertices:
        if i1_residual <= 1:
            continue
        rv1 = g.return_node(i1)
        other_foldback_found = 0
        if rv1.type == 'H':
            # extend while section has CN gain
            i_i, i_j = i1, i1+1
            while True:
                v_i, v_j = g.return_node(i_i), g.return_node(i_j)
                if v_i is None or v_j is None:
                    break
                if not segment_edge_is_amplified(i_i, i_j):
                    break
                foundback_multiplicity = self_edge_multiplicity(i_j)
                if foundback_multiplicity:
                    other_foldback_found = foundback_multiplicity
                    break
                i_i, i_j = i_i+2, i_j+2
        else:  # 'T'
            i_i, i_j = i1-1, i1
            while True:
                v_i, v_j = g.return_node(i_i), g.return_node(i_j)
                if v_i is None or v_j is None:
                    break
                if not segment_edge_is_amplified(i_i, i_j):
                    break
                foundback_multiplicity = self_edge_multiplicity(i_i)
                if foundback_multiplicity:
                    other_foldback_found = foundback_multiplicity
                    break
                i_i, i_j = i_i-2, i_j-2
        if other_foldback_found:
            g2.add_dummy_edges(i1, i1, other_foldback_found)

    residue_vertices = detect_residue_dgree(component, component_edges, g2.terminal_vertices_ids)
    print('RESIDUE vertices (after self-edge): ', residue_vertices)
    if len(residue_vertices) > 0:
        # compute pair wise distances for all pairs
        rv_distances = []
        for idx1, (i1, _) in enumerate(residue_vertices):
            rv1 = g.return_node(i1)
            for (i2, _) in residue_vertices[idx1+1:]:
                rv2 = g.return_node(i2)
                if rv1.chromosome == rv2.chromosome:
                    if (rv1.type, rv2.type) == ('H', 'T') or (rv1.type, rv2.type) == ('T', 'H'):
                        rv_distances.append((abs(rv1.pos - rv2.pos), i1, i2))
        # heuristic pairing
        rv_distances.sort()
        remaining_residuals = {i[0]: i[1] for i in residue_vertices}
        for d in rv_distances:
            edge_available_residuals = min(remaining_residuals[d[1]], remaining_residuals[d[2]])
            if edge_available_residuals:
                g2.add_dummy_edges(d[1], d[2], edge_available_residuals)
                remaining_residuals[d[1]] -= edge_available_residuals
                remaining_residuals[d[2]] -= edge_available_residuals

    odd_vertices = detect_segment_odd_degree(component, component_edges)
    print('ODD vertices', odd_vertices)
    # for i in range(0,len(segment_vertices), 2):
    #     g2.add_dummy_edges(segment_vertices[i], segment_vertices[i+1])
    # print(g2.edges)
    print('TOUR')
    print(segment_vertices)
    if len(odd_vertices) == 0:  # Eulerian circuites exist
        # a = printEulerUtil(g2, segment_vertices[0], -1)
        a = []
        for i in segment_vertices:
            a.append(printEulerUtil(g2, i, -1, g.return_node(i).chromosome))
        a = scoring_paths(a, segment_vertices, g, centro)
    elif len(odd_vertices) == 2:  # Eulerian path exists
        if odd_vertices[0] in segment_vertices:  # it is better to start path finding from telomere regions
            a = printEulerUtil(g2, odd_vertices[0], -1, g.return_node(odd_vertices[0]).chromosome)
        elif odd_vertices[1] in segment_vertices:
            a = printEulerUtil(g2, odd_vertices[1], -1, g.return_node(odd_vertices[0]).chromosome)
        else:  # If not exist add a dummy edges to make it Eulerian and then find from segment vertices
            g2.add_dummy_edges(odd_vertices[0], odd_vertices[1], 1)
            # a = printEulerUtil(g2, segment_vertices[0], -1)
            a = []
            for i in segment_vertices:
                a.append(printEulerUtil(g2, i, -1, g.return_node(i).chromosome))
            a = scoring_paths(a, segment_vertices, g, centro)
    else:  # if more than two vertices with odd degree. Connect them to each other to make the graph Eulerian
        if len(set(odd_vertices).intersection(
                set(segment_vertices))) == 0:  # no telomere nodes with odd degree connect all of them to gether like previouse setp
            for i in range(0, len(odd_vertices), 2):
                g2.add_dummy_edges(odd_vertices[i], odd_vertices[i + 1], 1)
            a = []
            for i in segment_vertices:
                a.append(printEulerUtil(g2, i, -1, g.return_node(i).chromosome))
            a = scoring_paths(a, segment_vertices, g, centro)
            # a = printEulerUtil(g2, segment_vertices[-1], -1)#Siavash in chaneed shode 
        else:
            count = 0
            save_index = 0
            for i in range(0, len(odd_vertices), 2):  # expect the segment vertices connect all to gether and then find Eulerian path from it
                if (odd_vertices[i] in segment_vertices or odd_vertices[i + 1] in segment_vertices) and count == 0:
                    count += 1
                    save_index = i
                else:
                    g2.add_dummy_edges(odd_vertices[i], odd_vertices[i + 1] , 1)
            if odd_vertices[save_index] in segment_vertices:
                a = printEulerUtil(g2, odd_vertices[save_index], -1, g.return_node(odd_vertices[save_index]).chromosome)
            else:
                a = printEulerUtil(g2, odd_vertices[save_index + 1], -1, g.return_node(odd_vertices[save_index]).chromosome)

    print('all edges with dummy: ', g2.edges)
    g2.print_node()
    print('Answer', a)

    with open(output_file, 'w') as fp_write:
        for edge in g2.edges:
            fp_write.write(str(edge) + '\n')

    return a, g2.edges


def detect_del_dup_cn(chromosome, start, end, segments):  # this function detect that for a deletion or duplication, do we have CNV call as well or no
    # chromosome, start, end position of deletion call are input
    for i, s in enumerate(segments):  # search in segments
        if int(s.chromosome) == int(chromosome):
            overlap = max(0, min(end, s.end) - max(start, s.start))
            # if these two windows are overlapping each other 0.9 and 1.2 of theirlength and the segment have CN different from one of it adjacent segments it return true.
            if overlap > 0.9 * (s.end - s.start):
                if (end - start) < 1.2 * (s.end - s.start):
                    if (s.int_cn != segments[(i + 1) % len(segments)].int_cn and s.chromosome == segments[(i + 1) % len(segments)].chromosome) or (
                            s.int_cn != segments[i - 1].int_cn and s.chromosome == segments[i - 1].chromosome):
                        # print('Kir', chromosome, start, end, overlap)
                        # print('Kir', start , end, s.start,s.end)
                        return True, s.start, s.end
            if overlap > 0.9 * (end - start):
                if (s.end - s.start) < 1.2 * (end - start):
                    if (s.int_cn != segments[(i + 1) % len(segments)].int_cn and s.chromosome == segments[(i + 1) % len(segments)].chromosome) or (
                            s.int_cn != segments[i - 1].int_cn and s.chromosome == segments[i - 1].chromosome):
                        return True, s.start, s.end
    return False, None, None


def detect_duplicatioon_inversion_cn(sv, xmap, segments):  # same as above for duplication calls.
    window_lim = 300000
    node_dir1, node_dir2 = detect_sv_directions(sv, xmap)
    if node_dir1 == 'T' and node_dir2 == 'T':  # right fold back
        for i, s in enumerate(segments):
            if int(s.chromosome) == int(sv.ref_c_id1):  # only compared with the next contigs
                if abs(s.end - sv.ref_end) < window_lim:
                    if (s.chromosome == segments[(i + 1) % len(segments)].chromosome):
                        if (s.int_cn != segments[(i + 1) % len(segments)].int_cn ):
                            return True, s.end
                    else:
                        if (s.int_cn != segments[(i - 1) % len(segments)].int_cn ) and (s.chromosome == segments[(i - 1) % len(segments)].chromosome):
                            return True, s.end
    elif node_dir1 == 'H' and node_dir2 == 'H':  # left foldback
        for i, s in enumerate(segments):  # compared with prev contigs
            if int(s.chromosome) == int(sv.ref_c_id1):
                if abs(s.start - sv.ref_start) < window_lim:
                    if s.chromosome == segments[i - 1].chromosome:
                        if(s.int_cn != segments[i - 1].int_cn ):
                            return True, s.start
                    else:
                        if (s.int_cn != segments[(i + 1) % len(segments)].int_cn ) and (s.chromosome == segments[(i + 1) % len(segments)].chromosome):
                            return True, s.start
    return False, None


def detect_receprical_translocation(sv, xmap, smap):  # sometimes one of these reciprocal translocation have low confidence but this function we retrive it
    window_lim = 200000
    for i in smap:
        if i.ref_c_id1 == sv.ref_c_id1 and i.ref_c_id2 == sv.ref_c_id2 and sv.smap_id != i.smap_id:
            if abs(i.ref_start - sv.ref_start) < window_lim and abs(i.ref_end - sv.ref_end) < window_lim:
                i_dir1, i_dir2 = detect_sv_directions(i, xmap)
                sv_dir1, sv_dir2 = detect_sv_directions(sv, xmap)
                if i_dir1 != sv_dir1 and i_dir2 != sv_dir2:
                    return True, i
    return False, None


def close_gaps_between_segments(segments):
    # Sort segments by chromosome and start position to ensure they are consecutive
    segments.sort(key=lambda seg: (int(seg.chromosome), seg.start))

    for i in range(len(segments) - 1):
        current_segment = segments[i]
        next_segment = segments[i + 1]

        # Check if both segments are on the same chromosome
        if current_segment.chromosome == next_segment.chromosome:
            # Calculate the gap between the end of the current segment and the start of the next segment
            gap = next_segment.start - current_segment.end

            if gap > 1:
                # Adjust the end of the current segment to close the gap
                current_segment.end = next_segment.start - 1
                # Update the bp list with the new start and end values
                current_segment.bp[-1] = current_segment.end

    return segments



def reverse_path(path):
    # Split the path into individual movements
    movements = path.split()
    # Reverse the list of movements to reverse the order
    reversed_movements = movements[::-1]
    # Initialize an empty list to store the reversed movements with reversed directions
    reversed_with_directions = []
    # Reverse the directions in each movement and add it to the reversed_with_directions list
    for movement in reversed_movements:
        direction = movement[-1]  # Get the last character (either '+' or '-')
        reversed_direction = '+' if direction == '-' else '-'  # Reverse the direction
        movement_with_direction = movement[:-1] + reversed_direction  # Create the reversed movement with direction
        reversed_with_directions.append(movement_with_direction)
    # Join the reversed movements with reversed directions to form the reversed path
    reversed_path = ' '.join(reversed_with_directions)
    return reversed_path


def find_indices_with_sum_of_2(numbers):
    n = len(numbers)
    ans = []
    for i in range(n):
        for j in range(i + 1, n):
            if numbers[i] + numbers[j] == 2:
                ans.append((i, j))
    return ans


def share_same_segments(path1, path2):
    # Split the paths into individual movements
    movements1 = path1.split()
    movements2 = path2.split()
    for i in movements1:
        for j in movements2:
            if i == j:
                return True
    return False


def convert_segment_to_path(p):
    ans = []
    for i in p:
        direction = i[-1]
        seg_number = int(i[:-1])
        numbers1, numbers2 = seg_number * 2 - 2, seg_number * 2 - 1
        if direction == '+':
            ans.append(numbers1)
            ans.append(numbers2)
        else:
            ans.append(numbers2)
            ans.append(numbers1)
    return ans


def swap_segment(p1, p2, g, centro):
    movements1 = p1.split()
    movements2 = p2.split()
    for i in range(len(movements1)):
        for j in range(len(movements2)):
            if movements1[i] == movements2[j]:
                new_p1 = movements2[:j] + movements1[i:]
                new_p2 = movements1[:i] + movements2[j:]
                if check_non_centromeric_path(convert_segment_to_path(new_p1), g, centro) == 1 and check_non_centromeric_path(convert_segment_to_path(new_p2),
                                                                                                                              g, centro) == 1:
                    return True, ' '.join(new_p1)+' ', ' '.join(new_p2)+' '
    return False, p1, p2


def fix_dicentric(paths, scores, g, centro):
    bad_path = []
    bad_scores = []
    ans_path = []
    ans_score = []
    for i in range(len(paths)):
        if scores[i] != 1:
            bad_path.append(paths[i])
            bad_scores.append(scores[i])
        else:
            ans_score.append(scores[i])
            ans_path.append(paths[i])
    list_of_indices = find_indices_with_sum_of_2(bad_scores)
    seen = []
    for pair_indices in list_of_indices:
        i, j = pair_indices[0], pair_indices[1]
        if i not in seen and j not in seen:
            condition = False
            if share_same_segments(bad_path[i], bad_path[j]):
                condition, p1, p2 = swap_segment(bad_path[i], bad_path[j], g, centro)
            elif share_same_segments(reverse_path(bad_path[i]), bad_path[j]):
                condition, p1, p2 = swap_segment(reverse_path(bad_path[i]), bad_path[j], g, centro)
            if condition:
                ans_path.append(p1)
                ans_score.append(1)
                ans_path.append(p2)
                ans_score.append(1)
                seen.append(i)
                seen.append(j)
    for i in range(len(bad_scores)):
        if i not in seen:
            ans_path.append(bad_path[i])
            ans_score.append(bad_scores[i])
    return ans_path, ans_score
def convert_string_to_path_direction(path):
    segments = path.strip().split()  # Split the string by spaces and remove any leading/trailing whitespace
    numbers = []
    directions = []
    for segment in segments:
        number = int(segment[:-1])  # Convert the number part to integer
        direction = segment[-1]  # Get the last character as the direction
        numbers.append(number)
        directions.append(direction)
    return numbers, directions
def calculate_path_length_from_centro(path, g, centro, arm):
    l = 0
    for p in path:
        id1 = int(p) *2 -2
        id2 = int(p)*2 -1
        l = l + (g.vertices[id2].pos - g.vertices[id1].pos)
        chrom = g.vertices[id1].chromosome
    if arm == 'p':
        if l / max(centro['chr'+str(chrom)]) < 0.9:
            return True
        else:
            return False
    if arm == 'q':
        end = 0
        for v in g.vertices:
            if v.chromosome == chrom and v.pos > end:
                end = v.pos
        if l / (end - max(centro['chr'+str(chrom)])) < 0.9:
            return True
        else:
            return False
    return False
def merge_path_with_0_score(ans2 , scores_path, g, centro):
    removed = set()
    for i in range(len(scores_path)):
        if scores_path[i] == 0:
            p = ans2[i]
            s1,d1 = convert_string_to_path_direction(p)
            for j in range(len(scores_path)):
                if j != i and scores_path[j] != 0:
                    s2, d2 = convert_string_to_path_direction(ans2[j])
                    print('Ger', ans2, scores_path, removed, i, s1, d1, s2 , d2,j)
                    if s1[0] == s2[-1] and d1[0]!=d2[-1] and calculate_path_length_from_centro(s1,g,centro,'p'):
                        ans2[j] = ans2[j]  + ans2[i]
                        scores_path[j]+= scores_path[i]
                        removed.add(i)
                        break
                    elif s1[-1] == s2[0] and d1[-1]!=d2[0] and calculate_path_length_from_centro(s1,g,centro,'q') :
                        ans2[j] = ans2[i] + ans2[j]
                        scores_path[j] += scores_path[i]
                        removed.add(i)
                        break

    ans2 =  [ans2[i] for i in range(len(ans2)) if i not in removed]
    scores_path = [scores_path[i] for i in range(len(scores_path)) if i not in removed]
    return  ans2 , scores_path

def convert_path_to_segment(p, component_edges,g):  # this is important function that convert Eulerion path with vertices ID to segment path.
    component = list(set(p))

    ## form dict for all cc edges to be efficient in search
    component_edge_dict = {}  # {(node1, node2, type)}: multiplicity
    for edge_itr in component_edges:
        component_edge_dict[(edge_itr[0], edge_itr[1], edge_itr[3])] = edge_itr[2]

    print(component_edge_dict)

    ## relabel the edges
    edge_labels = []

    def append_and_reduce_multiplicity(node1, node2, edge_type):
        edge_labels.append(edge_type)
        if component_edge_dict[(node1, node2, edge_type)] == 1:
            component_edge_dict.pop((node1, node2, edge_type))
        else:
            component_edge_dict[(node1, node2, edge_type)] -= 1

    # skip the first one because it is always Segment edge
    if (p[0], p[1], 'S') in component_edge_dict:
        append_and_reduce_multiplicity(p[0], p[1], 'S')
    elif (p[1], p[0], 'S') in component_edge_dict:
        append_and_reduce_multiplicity(p[1], p[0], 'S')
    else:
        raise RuntimeError('initialization error')

    for node_idx in range(1, len(p) - 1):
        current_node = p[node_idx]
        next_node = p[node_idx + 1]
        previous_edge_type = edge_labels[node_idx - 1]
        # if previous is transition (SV/REF), next must be SEG; if previous is SEG, prefer next to be transition
        if previous_edge_type in ['SV', 'R', 'D']:
            if (current_node, next_node, 'S') not in component_edge_dict and (next_node, current_node, 'S') not in component_edge_dict:
                # print(component_edge_dict)
                # print(edge_labels)
                raise RuntimeError(f'illegal follow up of SV/R, no S present: {component_edge_dict}; {edge_labels}; {p}; {component_edges}')
            else:
                if (current_node, next_node, 'S') in component_edge_dict:
                    append_and_reduce_multiplicity(current_node, next_node, 'S')
                elif (next_node, current_node, 'S') in component_edge_dict:
                    append_and_reduce_multiplicity(next_node, current_node, 'S')
        elif previous_edge_type in 'S':
            if (current_node, next_node, 'R') in component_edge_dict:
                append_and_reduce_multiplicity(current_node, next_node, 'R')
            elif (next_node, current_node, 'R') in component_edge_dict:
                append_and_reduce_multiplicity(next_node, current_node, 'R')
            elif (current_node, next_node, 'SV') in component_edge_dict:
                append_and_reduce_multiplicity(current_node, next_node, 'SV')
            elif (next_node, current_node, 'SV') in component_edge_dict:
                append_and_reduce_multiplicity(next_node, current_node, 'SV')
            elif (current_node, next_node, 'D') in component_edge_dict:
                append_and_reduce_multiplicity(current_node, next_node, 'D')
            elif (next_node, current_node, 'D') in component_edge_dict:
                append_and_reduce_multiplicity(next_node, current_node, 'D')
            elif (current_node, next_node, 'S') in component_edge_dict:
                append_and_reduce_multiplicity(current_node, next_node, 'S')
            elif (next_node, current_node, 'S') in component_edge_dict:
                append_and_reduce_multiplicity(next_node, current_node, 'S')
            else:
                raise RuntimeError('illegal follow up of S, no SV/R/D/S edge available (ie. no edge available)')

    ## split into paths: split whenever two adjacent edge types are both S
    # collect all node_idx that have both edges as S
    split_node_indices = []
    for edge_idx in range(1, len(edge_labels)):
        current_label = edge_labels[edge_idx]
        previous_label = edge_labels[edge_idx - 1]

        if current_label == 'S' and previous_label == 'S':
            split_node_indices.append(edge_idx)

    # split
    print('split node indices: ', split_node_indices)
    paths = []
    previous_node_idx_cutoff = 0
    for split_node_idx in split_node_indices:
        current_path = [p[node_idx] for node_idx in range(previous_node_idx_cutoff, split_node_idx + 1)]
        paths.append(current_path)
        previous_node_idx_cutoff = split_node_idx  # we want this node to be present in both the splitted paths, end of path1 + start of path2
    last_path = [p[node_idx] for node_idx in range(previous_node_idx_cutoff, len(p))]
    paths.append(last_path)

    print('debug paths separation: ', paths)
    # merge first and last path for cycle
    if p[0] == p[-1]:
        if edge_labels[0] != 'S' or edge_labels[-1] != 'S':
            # the start node is meant to be split if both are S
            new_path = paths[-1]
            for node_idx in range(1, len(paths[0])):
                # do not repeat the annealed node
                new_path.append(paths[0][node_idx])
            paths.pop(0)
            paths.pop(len(paths) - 1)
            paths.append(new_path)

    # convert to segment names
    scores_path = []
    ans2 = []
    for p in paths:
        scores_path.append(check_non_centromeric_path(p,g, centro))
        temp = ''
        for i in range(0, len(p) - 1, 2):
            # print('as', i, len(p), p)
            seg_number = int((max(p[i], p[i + 1]) + 1) / 2)
            direction = '+'
            if p[i] > p[i + 1]:
                direction = '-'
            temp = temp + str(seg_number) + direction + ' '
        ans2.append(temp)

    # return ans2, [-1 for _ in range(len(ans2))]
    ans2 , scores_path = fix_dicentric(ans2, scores_path, g, centro)
    ans2 , scores_path = merge_path_with_0_score(ans2 , scores_path, g, centro)
    return ans2 , scores_path
def check_exiest_call(chromosome, start, end, type, all_seg):  # if we have a call like gain or loss  but in CNV it is filtered retrive it
    for s in all_seg:
        if s.chromosome == chromosome and s.start >= start and s.end <= end and s.type[:4] == type[:4]:
            return True
    return False


def extend_segments_cn(segments,
                       all_seg):  # this function check that if between two cnv call gap is less than 400Kbp and there is a call in CNV call but it marked or filtered we assume it is true and extend the segment length
    start = True
    for i in range(0, len(segments) - 1):
        s = segments[i]
        next_seg = segments[i + 1]
        if s.chromosome == next_seg.chromosome:
            if abs(next_seg.start - s.end) < 400000:
                if check_exiest_call(s.chromosome, s.end, next_seg.start, s.type, all_seg):
                    s.end = next_seg.start - 1
                    s.bp = [s.start, s.end]
                    segments[i] = s
    return segments


def average_values_greater_than_a(input_dict, a):
    total = 0
    count = 0
    for key, value in input_dict.items():
        if key > a:
            total += value
            count += 1
    return total / count if count > 0 else None


def cn_in_mask_N_region(chromosome, start, end, cop, centro):
    if chromosome not in ['22', '21', '15', '14', '13']:  # Andy shared this chromosome with me that Bionano CNV call in P arm is not reliable
        # return 2
        return round(np.average(list(cop.values()))) #instead of assumption that it is always diploid return the average cn of that chromosome
    else:
        if chromosome == '22' and start < max(centro['chr22']) and end < max(centro['chr22']):
            return round(average_values_greater_than_a(cop, max(centro['chr22'])))
        if chromosome == '21' and start < max(centro['chr21']) and end < max(centro['chr21']):
            return round(average_values_greater_than_a(cop, max(centro['chr21'])))
        if chromosome == '15' and start < max(centro['chr15']) and end < max(centro['chr15']):
            return round(average_values_greater_than_a(cop, max(centro['chr15'])))
        if chromosome == '14' and start < max(centro['chr14']) and end < max(centro['chr14']):
            # print(chromosome,round(average_values_greater_than_a(cop, max(centro['chr14']))), 'ponzi')
            return round(average_values_greater_than_a(cop, max(centro['chr14'])))
        if chromosome == '13' and start < max(centro['chr13']) and end < max(centro['chr13']):
            return round(average_values_greater_than_a(cop, max(centro['chr13'])))
    return 2




def merge_segments_all_seg_smap(segments, all_seg, smap, centro):
    ans = []
    # limit = 50000
    limit = 200000
    for sv in smap:
        if sv.sv_type == 'deletion':  # and sv.ref_c_id1=='17':# and sv.ref_start > 20400000 and sv.ref_start < 21000000:
            a = 0
            for s in all_seg:
                if sv.ref_c_id1 == s.chromosome and s.type.startswith('loss') and s.width > limit:# and not s.type.endswith('masked'): # why do we have 200,000 here?  # should be same number az search Min Seg Length #adserqa
                    if abs(s.start - min(sv.ref_start, sv.ref_end)) < limit and abs(s.end - max(sv.ref_start, sv.ref_end)) < limit:
                        if not is_overlapping(min(centro['chr' + str(sv.ref_c_id1)]), max(centro['chr' + str(sv.ref_c_id1)]), sv.ref_start, sv.ref_end):
                            ans.append(s)
        elif sv.sv_type.startswith('dup'):
            for s in all_seg:
                if sv.ref_c_id1 == s.chromosome and s.type.startswith('gain') and s.width > limit:
                    if not sv.sv_type.endswith('inverted'):
                        if abs(sv.ref_start - s.start) < limit and abs(sv.ref_end - s.end) < limit:
                            if not is_overlapping(min(centro['chr' + str(sv.ref_c_id1)]), max(centro['chr' + str(sv.ref_c_id1)]), sv.ref_start, sv.ref_end):
                                ans.append(s)
                    elif sv.sv_type.endswith('inverted'):
                        if not is_overlapping(min(centro['chr' + str(sv.ref_c_id1)]), max(centro['chr' + str(sv.ref_c_id1)]), sv.ref_start, sv.ref_end):
                            if abs(s.start - np.mean([sv.ref_start ,sv.ref_end])) < limit:
                                ans.append(s)
                            elif abs(s.end - np.mean([sv.ref_start ,sv.ref_end])) < limit:
                                ans.append(s)
    for s in ans:
        if s not in segments:
            segments.append(s)
    return segments


def adjust_and_remove_overlapping_segments(all_segments):
    # Sort segments by chromosome and start position
    all_segments.sort(key=lambda seg: (int(seg.chromosome), seg.start))

    adjusted_segments = []
    current_chromosome = None
    prev_segment = None

    for segment in all_segments:
        # Check if we are still on the same chromosome
        if current_chromosome != segment.chromosome:
            current_chromosome = segment.chromosome
            prev_segment = None  # Reset previous segment when chromosome changes

        if prev_segment:
            # Check if the current segment is fully contained within the previous segment
            if prev_segment.start <= segment.start and prev_segment.end >= segment.end:
                # Skip adding this segment, as it's fully contained within the previous one
                continue
            elif segment.start <= prev_segment.end:
                # Adjust the start of the current segment to prevent overlap
                segment.start = prev_segment.end + 1  # Ensure no overlap by adding 1 to prev_segment.end
                segment.bp = [segment.start, segment.end]

        # Add the current segment to the adjusted list (either non-overlapping or adjusted)
        adjusted_segments.append(segment)
        # Update the previous segment reference
        prev_segment = segment

    return adjusted_segments


def fix_coordinate(segments, all_seg , smap):
    limit = 200000
    for s1 in range(len(all_seg)):
        dist = 999999999
        if all_seg[s1].type.startswith('loss'):
            for i in smap:
                if i.sv_type.startswith('dele'):
                    if i.ref_c_id1 == all_seg[s1].chromosome :
                        if abs(i.ref_start - all_seg[s1].start) < limit and abs(i.ref_end - all_seg[s1].end) < min(limit, dist):
                            dist = abs(i.ref_start - all_seg[s1].start) < limit and abs(i.ref_end - all_seg[s1].end)
                            all_seg[s1].end = i.ref_end
                            all_seg[s1].start = i.ref_start
                            all_seg[s1].bp = [i.ref_start, i.ref_end]
        elif all_seg[s1].type.startswith('gain'):
            for i in smap:
                if i.sv_type.startswith('dup'):
                    if i.ref_c_id1 == all_seg[s1].chromosome :#and str(i.ref_c_id1) == '10' :
                        if not i.sv_type.endswith('inverted'):
                            if abs(i.ref_start - all_seg[s1].start) < limit and abs(i.ref_end - all_seg[s1].end) < limit and dist > min(abs(i.ref_start - all_seg[s1].start),abs(i.ref_end - all_seg[s1].end)) :
                                dist = min(abs(i.ref_start - all_seg[s1].start),abs(i.ref_end - all_seg[s1].end))
                                all_seg[s1].end = i.ref_end
                                all_seg[s1].start = i.ref_start
                                all_seg[s1].bp = [i.ref_start, i.ref_end]
                        elif i.sv_type.endswith('inverted'): # need to check direction as well
                            if abs(all_seg[s1].start - np.mean([i.ref_start ,i.ref_end])) < limit and abs(all_seg[s1].start - np.mean([i.ref_start ,i.ref_end])) < dist:
                                dist = abs(all_seg[s1].start - np.mean([i.ref_start ,i.ref_end]))
                                all_seg[s1].start = min([i.ref_start ,i.ref_end])
                                all_seg[s1].bp = [min([i.ref_start ,i.ref_end]), all_seg[s1].bp[1]]
                                if s1 > 0 and all_seg[s1].start < all_seg[s1-1].end and all_seg[s1].chromosome == all_seg[s1-1].chromosome:
                                    if all_seg[s1-1].start < all_seg[s1].start:
                                        all_seg[s1-1].end = all_seg[s1].start -1
                                        all_seg[s1-1].bp = [all_seg[s1-1].start, all_seg[s1-1].end]
                                    else:
                                        all_seg[s1].start = all_seg[s1-1].end + 1
                                        all_seg[s1].bp = [all_seg[s1].start, all_seg[s1].bp[1]]
                            elif abs(all_seg[s1].end - np.mean([i.ref_start ,i.ref_end])) < limit and abs(all_seg[s1].end - np.mean([i.ref_start ,i.ref_end])) < dist:
                                dist = abs(all_seg[s1].end - np.mean([i.ref_start ,i.ref_end]))
                                all_seg[s1].end = max([i.ref_start, i.ref_end])
                                all_seg[s1].bp = [all_seg[s1].bp[0], max([i.ref_start, i.ref_end])]
                                if s1 + 1 < len(all_seg) and all_seg[s1].end > all_seg[s1+1].start and all_seg[s1].chromosome == all_seg[s1+1].chromosome:
                                    if all_seg[s1].end < all_seg[s1+1].end :
                                        all_seg[s1+1].start = all_seg[s1].end + 1
                                        all_seg[s1+1].bp = [all_seg[s1+1].start, all_seg[s1+1].end]
                                    else:
                                        all_seg[s1].end = all_seg[s1+1].start-1
                                        all_seg[s1].bp = [all_seg[s1].bp[0], all_seg[s1].end]
    for s1 in range(len(segments)):
        dist = 99999999
        if segments[s1].type.startswith('loss'):
            for i in smap:
                if i.sv_type.startswith('dele'):
                    if i.ref_c_id1 == segments[s1].chromosome :
                        if abs(i.ref_start - segments[s1].start) < limit and abs(i.ref_end - segments[s1].end) < min(limit, dist):
                            dist = abs(i.ref_start - segments[s1].start) < limit and abs(i.ref_end - segments[s1].end)
                            segments[s1].end = i.ref_end
                            segments[s1].start = i.ref_start
                            segments[s1].bp = [i.ref_start , i.ref_end]
        elif segments[s1].type.startswith('gain'):
            for i in smap:
                if i.sv_type.startswith('dup'):
                    if i.ref_c_id1 == segments[s1].chromosome:# and str(i.ref_c_id1) == '6' :
                        if not i.sv_type.endswith('inverted'):
                            if abs(i.ref_start - segments[s1].start) < limit and abs(i.ref_end - segments[s1].end) < limit and dist > min(abs(i.ref_start - segments[s1].start),abs(i.ref_end - segments[s1].end)):
                                segments[s1].end = i.ref_end
                                segments[s1].start = i.ref_start
                                segments[s1].bp = [i.ref_start , i.ref_end]
                        elif i.sv_type.endswith('inverted'):
                            if abs(segments[s1].start - np.mean([i.ref_start ,i.ref_end])) < limit and abs(segments[s1].start - np.mean([i.ref_start ,i.ref_end])) < dist:
                                dist = abs(segments[s1].start - np.mean([i.ref_start ,i.ref_end]))
                                segments[s1].start = min([i.ref_start ,i.ref_end])
                                segments[s1].bp = [segments[s1].start, segments[s1].end]
                                if s1 > 0 and segments[s1].start < segments[s1-1].end and segments[s1].chromosome == segments[s1-1].chromosome:
                                    if segments[s1 - 1].start < segments[s1].start:
                                        segments[s1-1].end = segments[s1].start -1
                                        segments[s1-1].bp = [segments[s1-1].start, segments[s1-1].end]
                                    else:
                                        segments[s1].start = segments[s1 - 1].end + 1
                                        segments[s1].bp = [segments[s1].start, segments[s1].bp[1]]
                            elif abs(segments[s1].end - np.mean([i.ref_start ,i.ref_end])) < limit and abs(segments[s1].end - np.mean([i.ref_start ,i.ref_end])) < dist:
                                dist = abs(segments[s1].end - np.mean([i.ref_start ,i.ref_end]))
                                segments[s1].end = max([i.ref_start, i.ref_end])
                                segments[s1].bp = [segments[s1].bp[0], max([i.ref_start, i.ref_end])]
                                if s1 + 1 < len(segments) and segments[s1].end > segments[s1+1].start and segments[s1].chromosome == segments[s1+1].chromosome:
                                    if segments[s1].end < segments[s1 + 1].end:
                                        segments[s1+1].start = segments[s1].end + 1
                                        segments[s1+1].bp = [segments[s1+1].start, segments[s1+1].end]
                                    else:
                                        segments[s1].end = segments[s1+1].start-1
                                        segments[s1].bp = [segments[s1].bp[0], segments[s1].end]

    # segments.sort(key=lambda x: (int(x.chromosome), x.start, x.end))
    # all_seg.sort(key=lambda x: (int(x.chromosome), x.start, x.end))
    # for s1 in range(1,len(all_seg)):
    #     if all_seg[s1-1].end > all_seg[s1].start:
    #         all_seg[s1].start = all_seg[s1-1].end + 1
    # for s1 in range(1,len(segments)):
    #     if segments[s1-1].end > segments[s1].start:
    #         segments[s1].start = segments[s1].end + 1
    return  adjust_and_remove_overlapping_segments(segments), adjust_and_remove_overlapping_segments(all_seg)

######################################################################################################################################
def run_omkar(cnv_path, smap_path, rcmap_path, xmap_path, centro_path, name, output_path):
    segments, all_seg = parse_cnvcall(cnv_path)
    smap = parse_smap(smap_path)
    segments, all_seg = fix_coordinate(segments, all_seg, smap)
    segments.sort(key=lambda x: (int(x.chromosome), x.start, x.end))
    all_seg.sort(key=lambda x: (int(x.chromosome), x.start, x.end))
    global centro
    centro = parse_centro(centro_path)
    segments = merge_segments_all_seg_smap(segments, all_seg, smap, centro)  # Need to debug this function
    segments.sort(key=lambda x: (int(x.chromosome), x.start))
    rcov, rcop = parse_rcmap(rcmap_path)
    chrY_cn = int(np.average(list(rcop['24'].values())) + 0.5)
    chrX_cn = round(np.average(list(rcop['23'].values())))
    print('dashag', np.average(list(rcop['23'].values())), chrX_cn, chrY_cn)
    if chrY_cn > 0:
        chrX_cn = max(1, chrX_cn)
    print('dashag', np.average(list(rcop['23'].values())), chrX_cn, chrY_cn)
    xmap = parse_xmap(xmap_path)
    repo_dir = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(f'{output_path}/{name}/', exist_ok=True)
    output = f"{output_path}/{name}/{name}_raw.txt"
    output2 = f"{output_path}/{name}/{name}_SV.txt"
    output3 = f"{output_path}/{name}/{name}_flagged.txt"
    output4 = f"{output_path}/{name}/{name}.txt"
    file = f"{output_path}/{name}/{name}.pdf"
    file2 = f"{output_path}/{name}/{name}_2.png"

    svs = []
    segments = extend_segments_cn(segments, all_seg)  # fill the gap between calls.
    segments_contains_masked_regions = []
    segments_border_masked_regions_start = []
    segments_border_masked_regions_end = []
    for k in rcop.keys():
        seg_list = []
        label_list = list(rcop[k].keys())
        for s in segments:
            if s.chromosome == k:
                if s.width > 200000:  # if call has length greater than 200Kbp assume a segment # should be same number az search Min Seg Length #adserqa
                    seg_list.append(s)
        prev_point = list(rcop[k].keys())[0]
        if len(seg_list) == 0:  # create segment for start of chromosme
            new_seg = Segments()
            new_seg.start = 0
            new_seg.end = list(rcop[k].keys())[-1]
            new_seg.width = list(rcop[k].keys())[-1]
            new_seg.chromosome = k
            new_seg.fractional_cn = cn_in_mask_N_region(k, new_seg.start, new_seg.end, rcop[k], centro)
            new_seg.int_cn = cn_in_mask_N_region(k, new_seg.start, new_seg.end, rcop[k], centro)  # assumption that sample is diploide. default CN = 2
            if int(k) == 23:
                new_seg.fractional_cn = chrX_cn
                new_seg.int_cn = chrX_cn
            if int(k) == 24:
                new_seg.fractional_cn = chrY_cn
                new_seg.int_cn = chrY_cn
            new_seg.bp = [0, list(rcop[k].keys())[-1]]
            segments.append(new_seg)
    for k in rcop.keys():
        seg_list = []
        label_list = list(rcop[k].keys())
        for s in segments:
            if s.chromosome == k:
                if s.width > 200000:  # if call has length greater than 200Kbp assume a segment
                    seg_list.append(s)
        prev_point = list(rcop[k].keys())[0]
        if len(seg_list) == 0:  # create segment for start of chromosme
            new_seg = Segments()
            new_seg.start = 0
            new_seg.end = list(rcop[k].keys())[-1]
            new_seg.width = list(rcop[k].keys())[-1]
            new_seg.chromosome = k
            new_seg.fractional_cn = cn_in_mask_N_region(k, new_seg.start, new_seg.end, rcop[k], centro)
            new_seg.int_cn = cn_in_mask_N_region(k, new_seg.start, new_seg.end, rcop[k], centro)  # assumption that sample is diploide. default CN = 2
            if int(k) == 23:
                new_seg.fractional_cn = chrX_cn
                new_seg.int_cn = chrX_cn
            if int(k) == 24:
                new_seg.fractional_cn = chrY_cn
                new_seg.int_cn = chrY_cn
            new_seg.bp = [0, list(rcop[k].keys())[-1]]
            segments.append(new_seg)
        else:
            seg_list.sort(key=lambda x: x.start)
            for s in seg_list:
                start, end = find_start_end(prev_point, s.start,
                                            label_list)  # there are labels between two segments then create segment with CN =2 between them
                if start != 0 and end != 0:  #
                    new_seg = Segments()
                    new_seg.start = start
                    new_seg.end = end
                    new_seg.width = end - start
                    new_seg.chromosome = k
                    new_seg.fractional_cn = cn_in_mask_N_region(k, new_seg.start, new_seg.end, rcop[k], centro)
                    new_seg.int_cn = cn_in_mask_N_region(k, new_seg.start, new_seg.end, rcop[k], centro)  # assumption that sample is diploide. default CN = 2
                    if int(k) == 23:
                        new_seg.fractional_cn = chrX_cn
                        new_seg.int_cn = chrX_cn
                    if int(k) == 24:
                        new_seg.fractional_cn = chrY_cn
                        new_seg.int_cn = chrY_cn
                    new_seg.bp = [start, end]
                    segments.append(new_seg)
                prev_point = s.end
            start, end = find_start_end(s.end, label_list[-1], label_list)
            if start != 0 and end != 0:
                new_seg = Segments()
                new_seg.start = start
                new_seg.end = end
                new_seg.width = end - start
                new_seg.chromosome = k
                new_seg.fractional_cn = cn_in_mask_N_region(k, new_seg.start, new_seg.end, rcop[k], centro)
                new_seg.int_cn = cn_in_mask_N_region(k, new_seg.start, new_seg.end, rcop[k], centro)
                if int(k) == 23:
                    new_seg.fractional_cn = chrX_cn
                    new_seg.int_cn = chrX_cn
                if int(k) == 24:
                    new_seg.fractional_cn = chrY_cn
                    new_seg.int_cn = chrY_cn
                new_seg.bp = [start, end]
                segments.append(new_seg)

    segments.sort(key=lambda x: (int(x.chromosome), x.start))
    segments = close_gaps_between_segments(segments)
    # for s in segments:
    #     print('asli', s.chromosome, s.start, s.end, s.int_cn, sorted(s.bp))
    with open(output2.replace('txt', 'bed'), 'w', newline='') as bed:
        bed_writer = csv.writer(bed, delimiter='\t')
        bed_writer.writerow(['chromosome', 'start', 'end', 'SV_type', 'smap_id', 'confidence'])
        with open(output2, 'w') as f:
            f.write(
                '#h\tSmapEntryID\tQryContigID\tRefcontigID1\tRefcontigID2\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tConfidence\tType\tXmapID1\tXmapID2\tLinkID\tQryStartIdx\tQryEndIdx\tRefStartIdx\tRefEndIdx\tZygosity\tGenotype\tGenotypeGroup\tRawConfidence\tRawConfidenceLeft\tRawConfidenceRight\tRawConfidenceCenter\tSVsize\tSVfreq\tOrientation\tVAF\n')
            for i in smap:
                if i.sv_type.startswith('dele') and not i.sv_type.endswith('nbase') and not i.sv_type.endswith(
                        'tiny') and i.size > 50000 and i.size < 1000000 and i.confidence > 0.8:
                    f.write(i.line)
                    bed_writer.writerow([i.ref_c_id1, i.ref_start, i.ref_end, i.sv_type, i.smap_id, i.confidence])
                elif i.sv_type == 'duplication' or i.sv_type == 'duplication_split' or i.sv_type == 'duplication_inverted':
                    f.write(i.line)
                    bed_writer.writerow([i.ref_c_id1, i.ref_start, i.ref_end, i.sv_type, i.smap_id, i.confidence])
                elif i.sv_type == 'inversion' and i.confidence >= 0.7:
                    f.write(i.line)
                    bed_writer.writerow([i.ref_c_id1, i.ref_start, i.ref_end, i.sv_type, i.smap_id, i.confidence])
        f.close()
    for i in smap:
        # translocation applied filters.
        if i.sv_type.startswith('trans') and i.confidence >= 0.05 and not i.sv_type.endswith('common') and not i.sv_type.endswith('segdupe') and (
                i.ref_c_id1 != i.ref_c_id2 or abs(i.ref_end - i.ref_start) > 300000):
            svs.append(i)  # fixed
            exist, s = detect_receprical_translocation(i, xmap, smap)
            if exist:
                svs.append(s)
                print(s.line.strip())
            print(i.line.strip())
        # indels
        elif i.sv_type.startswith('inse') or i.sv_type.startswith('delet'):
            # elif i.sv_type.startswith('delet'):
            if not i.sv_type.endswith('nbase') and not i.sv_type.endswith('tiny') and i.confidence >= 0:
                if i.sv_type.startswith('delet') and i.size > 200000:
                    print(i.line.strip())
                    if detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end, segments)[0]:
                        _, i.ref_start, i.ref_end = detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end, segments)
                        svs.append(i)
                        print(i.line.strip())
                    elif i.confidence >= 0.98 and i.size > 2000000:  # This is a limit for detection thos misslabeld translocation as deletion in Bionano Pipeline
                        #     # _, i.ref_start, i.ref_end = detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end, segments)
                        svs.append(i)
                    #     print(i.line.strip())
                    # print('Gerye', i.line.strip())
                elif i.size > 500000 and abs(i.ref_end - i.ref_start) > 500000:  # this would be for insertion length more than 500Kbp
                    svs.append(i)
                    print(i.line.strip())
        # if we have inversion SV
        elif i.sv_type == 'inversion' and i.confidence >= 0.7:  # filter low confidance
            start, end = 0, 0
            dir = ''
            dir1, dir2 = detect_sv_directions(i, xmap)
            s = find_in_smap(i.linkID, smap)  # inversion has two rwo in smap file. we find them with Link ID
            if dir1 == 'H':  # update inversion call. it is to complicated but baisically calculate inversion start and endpoint
                dir = 'left'
                start = s.ref_start
                end = i.ref_end
            else:
                dir = 'right'
                start = i.ref_start
                end = s.ref_start
            start, end = min(start, end), max(start, end)
            i.ref_start = start
            i.ref_end = end
            if abs(end - start) > 400000:  # apply filter on size of inversion
                svs.append(i)
                print(i.line.strip(), start, end, dir)
                print(s.line.strip())
        elif i.sv_type == 'inversion_paired' and i.confidence >= 0.7:  # if it is full inversion
            s = find_in_smap(i.linkID, smap)
            if abs(i.ref_start - s.ref_end) > 500000:
                i.ref_end = s.ref_end
                print(i.line.strip())
                svs.append(i)
        elif i.sv_type == 'duplication_inverted':
            if detect_duplicatioon_inversion_cn(i, xmap, segments)[0]:
                _, fold_point = detect_duplicatioon_inversion_cn(i, xmap, segments)  # because CN is changed
                i.ref_start = fold_point
                i.ref_end = fold_point
                svs.append(i)
                print(i.line.strip())
        elif i.sv_type == 'duplication' or i.sv_type == 'duplication_split':
            if detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end, segments)[0]:
                _, i.ref_start, i.ref_end = detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end, segments)
                svs.append(i)
                print(i.line.strip())
        # elif i.sv_type == 'duplication_split': #This maybe deleted, and get back to duplocation_split
        #     # if detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end)[0]:
        #     #     _, i.ref_start, i.ref_end = detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end)
        #         svs.append(i)
        #         print(i.line.strip())
        elif i.size > 500000 and not i.sv_type.startswith('inversion'):  # Other type of SV
            svs.append(i)
            print(i.line.strip())

    for sv in svs:  # integrate BPs and Segments
        find_bp_in_segment(sv.ref_c_id1, sv.ref_start, segments)  #
        find_bp_in_segment(sv.ref_c_id2, sv.ref_end, segments)
    # merging Bps for spiliting segments
    for s in segments:
        s.bp = merge_list(s.bp)

    # Graph Creation
    aa = 0
    g = Graph()
    counter = 0
    prev_chr = 0
    for index, s in enumerate(segments):  # Forach segment creates two vertices
        for inedx_bp, i in enumerate(s.bp):
            if inedx_bp == 0:  # if it is first BP in a segment
                aa += 1
                v = Vertices()  # Create Vetrices Object
                v.chromosome = s.chromosome
                v.id = counter
                v.pos = i
                v.cn = s.int_cn
                v.type = 'H'  # the start node is Head node
                g.vertices.append(v)
                if prev_chr == s.chromosome:
                    g.edges.append((counter - 1, counter, 0, 'R'))  # addign reference edge between two continues segments
                    g.return_node(counter - 1).append_edges(counter)  # update adjancency matrix
                    g.return_node(counter).append_edges(counter - 1)
                prev_chr = s.chromosome
                counter += 1
            elif inedx_bp == len(s.bp) - 1:  # if it is last BP in a segment
                aa += 1
                v = Vertices()
                v.chromosome = s.chromosome
                v.id = counter
                v.pos = i
                v.cn = s.int_cn
                v.type = 'T'  # last node should be Tail node
                g.vertices.append(v)
                g.edges.append((counter - 1, counter, s.int_cn, 'S'))
                g.return_node(counter).append_edges(counter - 1)
                g.return_node(counter - 1).append_edges(counter)
                # g.edges.append((counter,counter-1, s.int_cn,'S'))
                counter += 1
                prev_chr = s.chromosome
            else:  # Create both Tail and Head node
                aa += 1
                v = Vertices()
                v.chromosome = s.chromosome
                v.id = counter
                v.pos = i
                v.cn = s.int_cn
                v.type = 'T'
                g.vertices.append(v)
                g.edges.append((counter, counter - 1, s.int_cn, 'S'))
                g.return_node(counter - 1).append_edges(counter)
                g.return_node(counter).append_edges(counter - 1)
                # g.edges.append((counter-1,counter, s.int_cn,'S'))
                counter += 1
                ####
                aa += 1
                v = Vertices()
                v.chromosome = s.chromosome
                v.id = counter
                v.pos = i + 1
                v.type = 'H'
                v.cn = s.int_cn
                g.vertices.append(v)
                g.edges.append((counter - 1, counter, 0, 'R'))
                g.return_node(counter - 1).append_edges(counter)
                g.return_node(counter).append_edges(counter - 1)
                counter += 1
    for sv in svs:
        n_type1, n_type2 = detect_sv_directions(sv, xmap)
        # a and b are two nodes that ara connected by sv edge
        a = find_nodes(sv.ref_c_id1, sv.ref_start, g.vertices, n_type1)
        b = find_nodes(sv.ref_c_id2, sv.ref_end, g.vertices, n_type2)
        if int(sv.ref_c_id1) > int(sv.ref_c_id2) or (int(sv.ref_c_id1) == int(sv.ref_c_id2) and sv.ref_end < sv.ref_start):
            a = find_nodes(sv.ref_c_id1, sv.ref_start, g.vertices, n_type2)
            b = find_nodes(sv.ref_c_id2, sv.ref_end, g.vertices, n_type1)
        # if sv.sv_type == 'inversion' or sv.sv_type == 'inversion_paired':
        if sv.sv_type == 'inversion_paired':  # Lets complete the inversion
            if g.return_node(a).type == 'H':
                new_edge = (a - 1, b - 1, 0, 'SV')
            else:
                new_edge = (a + 1, b + 1, 0, 'SV')
            if new_edge not in g.edges:
                g.edges.append(new_edge)
                g.return_node(new_edge[0]).append_edges(new_edge[1])
                g.return_node(new_edge[1]).append_edges(new_edge[0])
        if b == a and (a, b, 0, 'SV') not in g.edges:  # it can be happend in duplication inversion
            g.edges.append((a, b, 0, 'SV'))
            g.return_node(a).append_edges(b)
        #     print(sv.line)
        #     print(a,b)
        elif sv.sv_type.startswith('delet'):  # telomere cite deletion prevent
            if (a, b, 0, 'SV') not in g.edges and abs(a - b) != 1:
                g.edges.append((a, b, 0, 'SV'))
                g.return_node(a).append_edges(b)
                g.return_node(b).append_edges(a)
        elif (a, b, 0, 'SV') not in g.edges:
            g.edges.append((a, b, 0, 'SV'))
            g.return_node(a).append_edges(b)
            g.return_node(b).append_edges(a)
        # g.edges.append((b,a,0,'SV'))
    g.print_node()
    g.output_node(f"{output_path}/{name}/{name}.preILP_nodes.txt")
    print('Siavash')
    print(g.edges)
    g.output_edges(f"{output_path}/{name}/{name}.preILP_edges.txt")
    Plot_graph(g, file, name)
    connected_components = find_connected_components(g)
    for component in connected_components:
        # if 102 in component:
        component_edges = estimating_edge_multiplicities_in_CC(component, g, xmap)
    connected_components = find_connected_components(g)
    Plot_graph(g, file2, name)
    paths = []
    edges_with_dummy = []

    component_counter = 0
    component_metadata = {}
    os.makedirs(f"{output_path}/{name}/postILP_components/", exist_ok=True)
    os.makedirs(f"{output_path}/{name}/all_edges_with_dummies/", exist_ok=True)
    print(g.vertices)
    for component in connected_components:
        component_metadata[component_counter] = component
        # if 102 in component:
        component_edges = return_all_edges_in_cc(component, g)
        print(component)
        print(component_edges)
        out_file = f"{output_path}/{name}/postILP_components/{name}.postILP_component_{component_counter}.txt"
        with open(out_file, 'w') as fp_write:
            fp_write.write(str(component) + '\n')
            for edge_itr in component_edges:
                fp_write.write(str(edge_itr) + '\n')

        out_file = f"{output_path}/{name}/all_edges_with_dummies/{name}.with_dummies_component_{component_counter}.txt"
        euler_tour, component_edges_with_dummies = printEulerTour(component, component_edges, g, out_file)
        paths.append(euler_tour)
        edges_with_dummy.append(component_edges_with_dummies)
        component_counter += 1
    metadata_file = f"{output_path}/{name}/postILP_components/{name}.postILP.metadata.txt"
    with open(metadata_file, 'w') as fp_write:
        for key, value in component_metadata.items():
            fp_write.write("{}\t{}\n".format(key, value))

    # print(paths)
    # write in the output
    with open(output3, 'w') as f2:
        with open(output, 'w') as f:
            f.write('Segment\tNumber\tChromosome\tStart\tEnd\tStartNode\tEndNode\n')
            f2.write('Segment\tNumber\tChromosome\tStart\tEnd\tStartNode\tEndNode\n')
            number = 1
            for i in range(0, len(g.vertices), 2):
                v = g.vertices[i]
                u = g.vertices[i + 1]
                f.write('Segment\t{id}\t{chrom}\t{start}\t{end}\t{snode}\t{enode}\n'.format(id=number, chrom=v.chromosome, start=v.pos, end=u.pos, snode=v.id,
                                                                                            enode=u.id))
                f2.write('Segment\t{id}\t{chrom}\t{start}\t{end}\t{snode}\t{enode}\n'.format(id=number, chrom=v.chromosome, start=v.pos, end=u.pos, snode=v.id,
                                                                                             enode=u.id))
                number += 1
            c = 1
            for path_idx, p in enumerate(paths):
                structures, scores = convert_path_to_segment(p, edges_with_dummy[path_idx], g)
                for jj in range(len(structures)):
                    structure = structures[jj]
                    if structure.endswith(' '):
                        structure = structure[:-1]
                    # f.write('Path'+str(c)+'='+','.join(str(z) for z in p)+'\n')
                    # print('path',p,check_non_centromeric_path(p,g, centro))
                    f.write('Path' + str(c) + ' = ' + structure + '\t score = ' + str(scores[jj]) + '\n')
                    paths_structure = structure.split(' ')
                    paths_structure2 = structure.split(' ')
                    for i in range(len(paths_structure)):
                        direction = paths_structure[i][-1]
                        seg = int(paths_structure[i][:-1])
                        if (direction == '-' and seg in segments_border_masked_regions_end) and i > 0:
                            prev_dir = paths_structure2[i - 1][-1]
                            prev_seg = int(paths_structure2[i - 1][:-1])
                            if prev_dir != direction or prev_seg - seg != 1:
                                paths_structure[i] = "*" + paths_structure[i]
                        if (direction == '-' and seg in segments_border_masked_regions_start) and i < len(paths_structure) - 1:
                            next_dir = paths_structure2[i + 1][-1]
                            next_seg = int(paths_structure2[i + 1][:-1])
                            if next_dir != direction or next_seg - seg != -1:
                                paths_structure[i] += "*"
                        if (direction == '+' and seg in segments_border_masked_regions_start) and i > 0:
                            prev_dir = paths_structure2[i - 1][-1]
                            prev_seg = int(paths_structure2[i - 1][:-1])
                            if prev_dir != direction or prev_seg - seg != -1:
                                paths_structure[i] = "*" + paths_structure[i]
                        if (direction == '+' and seg in segments_border_masked_regions_end) and i < len(paths_structure) - 1:
                            next_dir = paths_structure2[i + 1][-1]
                            next_seg = int(paths_structure2[i + 1][:-1])
                            if next_dir != direction or next_seg - seg != 1:
                                paths_structure[i] += "*"
                    structure = ' '.join(i for i in paths_structure)
                    f2.write('Path' + str(c) + ' = ' + structure + '\t score = ' + str(scores[jj]) + '\n')
                    c += 1

    ## post-processing
    post_process_input = output
    post_process_output = output4
    validation_code = validate_OMKar_output_format(output,
                                                   cont_allowance=200000,
                                                   span_allowance=200000,
                                                   forbidden_region_file=get_metadata_file_path('acrocentric_telo_cen.bed'))
    if validation_code:
        path_list, segment_dict = read_OMKar_output(post_process_input, return_segment_dict=True)
        processed_path_list, segment_obj_to_idx_dict = post_process_OMKar_output(path_list,
                                                                                 gap_merge_allowance=5,
                                                                                 isolate_centromere=True,
                                                                                 forbidden_region_file=get_metadata_file_path('acrocentric_telo_cen.bed'))
        write_MK_file(post_process_output, processed_path_list, segment_obj_to_idx_dict)
    else:
        ## if validation failed, keep original omkar output
        print(f"OMKar output format validation failed, the interpretation may include un-intended calls; check log for detailed reason; file: {post_process_output}", file=sys.stderr)
        shutil.copy2(post_process_input, post_process_output)

def find_input_file_paths(dir):
    """
    auto detect between Bionano Solve output structure, OR curated output structure
    :param dir: master dir of the sample containing all input files
    :return: the file paths for the four input files
    """
    cnv_filename = 'cnv_calls_exp.txt'
    rcmap_filename = 'cnv_rcmap_exp.txt'
    xmap_filename = 'exp_refineFinal1_merged.xmap'
    smap_filename = 'exp_refineFinal1_merged_filter_inversions.smap'

    # bionano default filepath
    cnv_filepath = f"{dir}/contigs/alignmolvref/copynumber/{cnv_filename}"
    rcmap_filepath = f"{dir}/contigs/alignmolvref/copynumber/{rcmap_filename}"
    xmap_filepath = f"{dir}/contigs/exp_refineFinal1_sv/merged_smaps/{xmap_filename}"
    smap_filepath = f"{dir}/contigs/exp_refineFinal1_sv/merged_smaps/{smap_filename}"
    if os.path.isfile(cnv_filepath) and os.path.isfile(rcmap_filepath) and os.path.isfile(xmap_filepath) and os.path.isfile(smap_filepath):
        return {'cnv': cnv_filepath,
                'rcmap': rcmap_filepath,
                'xmap': xmap_filepath,
                'smap': smap_filepath,
                'dir_structure': 'bionano_default'}

    # curated dir structure
    if os.path.isfile(f"{dir}/{cnv_filename}") and os.path.isfile(f"{dir}/{rcmap_filename}") and os.path.isfile(f"{dir}/{xmap_filename}") and os.path.isfile(f"{dir}/{smap_filename}"):
        cnv_filepath, rcmap_filepath, xmap_filepath, smap_filepath = f"{dir}/{cnv_filename}", f"{dir}/{rcmap_filename}", f"{dir}/{xmap_filename}", f"{dir}/{smap_filename}"
        return {'cnv': cnv_filepath,
                'rcmap': rcmap_filepath,
                'xmap': xmap_filepath,
                'smap': smap_filepath,
                'dir_structure': 'curated'}

    # debug use
    bionano_status1 = os.path.isfile(f'{dir}/contigs/alignmolvref/copynumber/{cnv_filename}')
    bionano_status2 = os.path.isfile(f'{dir}/contigs/alignmolvref/copynumber/{rcmap_filename}')
    bionano_status3 = os.path.isfile(f'{dir}/contigs/exp_refineFinal1_sv/merged_smaps/{xmap_filename}')
    bionano_status4 = os.path.isfile(f'{dir}/contigs/exp_refineFinal1_sv/merged_smaps/{smap_filename}')
    curated_status1 = os.path.isfile(f'{dir}/{cnv_filename}')
    curated_status2 = os.path.isfile(f'{dir}/{rcmap_filename}')
    curated_status3 = os.path.isfile(f'{dir}/{xmap_filename}')
    curated_status4 = os.path.isfile(f'{dir}/{smap_filename}')
    print(f"attempted {dir}/contigs/alignmolvref/copynumber/{cnv_filename}: {bionano_status1}")
    print(f"attempted {dir}/contigs/alignmolvref/copynumber/{rcmap_filename}: {bionano_status2}")
    print(f"attempted {dir}/contigs/exp_refineFinal1_sv/merged_smaps/{xmap_filename}: {bionano_status3}")
    print(f"attempted {dir}/contigs/exp_refineFinal1_sv/merged_smaps/{smap_filename}: {bionano_status4}")
    print(f"attempted {dir}/{cnv_filename}: {curated_status1}")
    print(f"attempted {dir}/{rcmap_filename}: {curated_status2}")
    print(f"attempted {dir}/{xmap_filename}: {curated_status3}")
    print(f"attempted {dir}/{smap_filename}: {curated_status4}")

    raise RuntimeError(f'file structure could not be parsed: {dir}')

def main():
    repo_dir = os.path.dirname(os.path.abspath(__file__))
    default_stdout = sys.stdout
    default_centro_path = f'{repo_dir}/hg38_centro.txt'

    parser = argparse.ArgumentParser()
    parser.add_argument("-dir", "--dir", help="path to directory containing input data", required=True)
    parser.add_argument("-centro", "--centro", default=default_centro_path, help="path to file contains centromere coordinates", required=False)
    parser.add_argument("-o", "--output", help="path to output dir", required=True)
    parser.add_argument("-single", "--single", help="flag to singal only run on one sample", action='store_true')
    parser.add_argument("-report", "--report", help="flag to generate html report", action='store_true')
    parser.add_argument("-noImage", "--noImage", help="flag to disable image generation", action='store_true')
    parser.add_argument("-reportDebug", "--reportDebug", help="generate report debug in html report", action='store_true')
    args = parser.parse_args()

    os.makedirs(f'{args.output}/omkar_output/', exist_ok=True)
    os.makedirs(f'{args.output}/logs/', exist_ok=True)
    if args.single:
        filepaths = find_input_file_paths(args.dir)
        sample_name = os.path.basename(os.path.normpath(args.dir))
        print(f'running sample: {sample_name}')
        sys.stdout = open(f"{args.output}/logs/{sample_name}.stdout.txt", 'w')
        sample_output_dir = f"{args.output}/omkar_output/"
        run_omkar(filepaths['cnv'], filepaths['smap'], filepaths['rcmap'], filepaths['xmap'], args.centro, sample_name, sample_output_dir)
        sys.stdout.close()
        sys.stdout = default_stdout
    else:
        for sample_dir in os.listdir(args.dir):
            filepaths = find_input_file_paths(f"{args.dir}/{sample_dir}/")
            sample_name = sample_dir
            print(f'running sample: {sample_name}')
            sys.stdout = open(f"{args.output}/logs/{sample_name}.stdout.txt", 'w')
            sample_output_dir = f"{args.output}/omkar_output/"
            run_omkar(filepaths['cnv'], filepaths['smap'], filepaths['rcmap'], filepaths['xmap'], args.centro, sample_name, sample_output_dir)
            sys.stdout.close()
            sys.stdout = default_stdout

    if args.report:
        from KarReporter import generate_html_report
        generate_image = not args.noImage
        sys.stdout = open(f"{args.output}/logs/__report.stdout.txt", 'w')
        os.makedirs(f'{args.output}/omkar_report/', exist_ok=True)
        omkar_input_dir_for_report = args.dir if not args.single else str(Path(args.dir).parent)
        generate_html_report(generate_image,
                             None,
                             '',
                             f'{args.output}/omkar_output/',
                             f'{args.output}/omkar_report/clustered_karyotype_images/',
                             f'{args.output}/omkar_report/',
                             omkar_input_dir_for_report,
                             debug=args.reportDebug)
        sys.stdout.close()
        sys.stdout = default_stdout
if __name__ == "__main__":
    main()
