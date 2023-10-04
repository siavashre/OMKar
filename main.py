from collections import defaultdict
from os.path import exists
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

class Vertices: #Class of vertices in out graph
    def __init__(self):
        self.type = '' #Type of vertices it can be Head or Tail 
        self.id = 0 # Id of vertices
        self.chromosome = 0 # which chromosome this vertices belong chrX = 23 and chrY = 24
        self.pos = 0 # in terms of bp the position of this vertices
        self.cn = 0 # integer CN of segment this vertices represent
        self.edges = [] # its like adjacenty matrix it contains all adjacent vertices

    def append_edges(self, i): #add new vertices to adjancty list
        self.edges.append(i)
        self.edges = list(set(self.edges))

    def remove_edges(self, i): # remove vertices from adjacenty list
        if i in self.edges:
            self.edges.remove(i)
        else:
            print('Not here baba')

    def print_v(self): #print vertices info
        print(self.id, self.chromosome, self.pos, self.cn, self.type, self.edges)


class Graph: #Class of graph 
    def __init__(self):
        self.vertices = [] #list of Object of vertices
        self.edges = [] #list of edges. each edges is tuple of (u,v,M,type) u and v are vertices connected to each other M represent edge multiplicities and type can be "S", "SV", "R", "D"

    def print_node(self): #print all nodes
        for v in self.vertices:
            print(v.id, v.chromosome, v.pos, v.cn, v.type, v.edges)

    def return_node(self, id): # return node by id
        for v in self.vertices:
            if id == v.id:
                return v
        return None

    def return_edges(self, a, b): #return all type of edge between two nodes. it return a list of edges or None
        ans = []
        for e in self.edges:
            if (a == e[0] and b == e[1]) or (a == e[1] and b == e[0]):
                ans.append(e)
        if len(ans)>0:
            return ans
        return None

    def add_dummy_edges(self, u, v):#add dummy edges with type "D" between two nodes if already edge exist between them increase the CN 
        e = self.return_edges(u, v)
        if max(u,v)%2 == 1 and abs(u-v) == 1:
           if len(e) == 1:
                e = e[0]
                cn = e[2] - 1
                self.edges.remove(e)
                self.edges.append((e[0], e[1], cn, e[3]))
        elif e == None:
            self.edges.append((u, v, 1, 'D'))
            self.return_node(u).append_edges(v)
            self.return_node(v).append_edges(u)
        else:
            if len(e) == 1:
                e = e[0]
                cn = e[2] + 1
                self.edges.remove(e)
                self.edges.append((e[0], e[1], cn, e[3]))
            else:
                e = e[0]
                cn = e[2] + 1
                self.edges.remove(e)
                self.edges.append((e[0], e[1], cn, e[3]))

    def update_edges(self, a, b, count,e_type): #update an edge between node a nad b with type = e_type to the new number
        e_list = self.return_edges(a, b)
        for e in e_list:
            if e[3] == e_type:
                self.edges.remove(e)
                if count > 0:
                    self.edges.append((e[0], e[1], count, e[3]))


def find_bp_in_segment(chromosome, point): #this function get chromosome and position as input(from sv call) and add this point to segment bp list. It helps us to then spliting segments to new one.
    for i in segments:
        if str(i.chromosome) == str(chromosome):
            if i.start <= point <= i.end:
                i.bp.append(point)



def merge_list(l): # this function merge all breakpoints in a segment within a window of 50Kbp
    l = list(set(l))
    l = sorted(l)
    if len(l) == 2:
        return l
    WINDOW = 50000
    ans = []
    prev = None
    group = []
    for item in l:
        if prev is None or abs(item - prev) <= WINDOW:
            group.append(item)
        else:
            ans.append(np.mean(group))
            group = [item]
        prev = item
    if len(group) > 0:
        ans.append(np.mean(group))
    return ans

######################### DOUBLE CHECK THIS FUNCTION ######################################
def detect_sv_directions(sv): #this function detect the direction of one smap like H/T to H/T. #need to be check Contain bug?
    if sv.sv_type == 'inversion_paired': # for inversion_paired always return T to T
        return 'T', 'T'
    dir1, dir2 = '', ''
    xmap_id1 = sv.xmap_id1 
    xmap_id2 = sv.xmap_id2
    alignment1 = xmap[str(xmap_id1)]
    alignment2 = xmap[str(xmap_id2)]
    swap = 0
    if min(alignment1['QryStartPos'], alignment1['QryEndPos']) > min(alignment2['QryStartPos'],
                                                                     alignment2['QryEndPos']):
        swap = 1
        alignment1, alignment2 = alignment2, alignment1 #sort alignment1 and alignment2 by query position. Always alignment1 has smaller position on the query
    if alignment1['Orientation'] == '+':
        dir1 = 'T'
    else:
        dir1 = 'H'
    if alignment2['Orientation'] == '+':
        dir2 = 'H'
    else:
        dir2 = 'T'
    dir1_a = dir1
    dir2_a = dir2
    dir1_b = dir1 
    dir2_b = dir2
    if swap == 1:
        dir1_a, dir2_a = dir2, dir1
    if int(alignment1['RefContigID']) > int(alignment2['RefContigID']) or (int(alignment1['RefContigID']) == int(alignment2['RefContigID'])
         and alignment1['RefStartPos']> alignment2['RefStartPos']):
        dir1_b, dir2_b = dir2, dir1
    if dir1_a !=dir1_b:
        print('Here we had BUG BUG BUG',sv)
    return dir1_b, dir2_b


def find_nodes(chromosome, pos, vertices, node_type): #find the closest node and return it id to the chromosome and position and type
    dist = 999999999 #as input take the chromosome number and position and type of node we are looking to it(H or T) and return it.
    ans = -1
    for v in vertices:
        if v.chromosome == chromosome and v.type == node_type:
            if abs(v.pos - pos) < dist:
                dist = abs(v.pos - pos)
                ans = v.id
    return ans


def find_start_end(prev_point, start, label_list):
    ans = []
    for i in label_list:
        if prev_point < i < start:
            ans.append(i)
    if len(ans) < 2:
        return 0, 0
    else:
        return min(ans), ans[-1]


def next_prev_label(chromosome, pos): #not used can be deleted
    for k in rcop.keys():
        label_list = list(rcop[k].keys())
        if k == chromosome:
            for index, i in enumerate(label_list):
                if i == pos:
                    return label_list[index - 1], label_list[(index + 1) % (len(label_list))]


def detect_overlap_map(chromosome, pos): # this function detect if there is a map(alignment that overlap this region)
    # maybe it is good idea instead of 1bp assume a window here.
    for xmapid in xmap.keys():
        x = xmap[xmapid]
        if int(chromosome) == int(x['RefContigID']):
            if  x['RefStartPos'] <=float(pos)<= x['RefEndPos']:
                return True
    return False

def Plot_graph(g, file, name): #this function plot the graph
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
                plt.annotate('(_,+)', xy=(node1.pos, node1.cn + 0.5+i))
            elif node1.type == 'H' and node2.type == 'T':
                plt.annotate('(_,_)', xy=(node1.pos, node1.cn + 0.5+i))
            elif node1.type == 'T' and node2.type == 'T':
                plt.annotate('(+,_)', xy=(node1.pos, node1.cn + 0.5+i))
            elif node1.type == 'T' and node2.type == 'H':
                plt.annotate('(+,+)', xy=(node1.pos, node1.cn + 0.5+i))
            exec(f"plt.subplot(grid{[int(node2.chromosome) - 1]})")
            j+=0.5
            if node1.type == 'H' and node2.type == 'H':
                plt.annotate('(_,+)', xy=(node2.pos, node2.cn + 0.5+j))
            elif node1.type == 'H' and node2.type == 'T':
                plt.annotate('(_,_)', xy=(node2.pos, node2.cn + 0.5+j))
            elif node1.type == 'T' and node2.type == 'T':
                plt.annotate('(+,_)', xy=(node2.pos, node2.cn + 0.5+j))
            elif node1.type == 'T' and node2.type == 'H':
                plt.annotate('(+,+)', xy=(node2.pos, node2.cn + 0.5+j))
    plt.savefig(file, dpi=200)


def find_in_smap(id): # return a smap call with the iD otherwise return None
    for s in smap:
        if s.smap_id == id:
            return s
    return None


def dfs(i, temp, g, visited): #run dfs alg on graph g on vertices i. visited is a list of seen vertices before visiting this node.
    visited.append(i)
    temp.append(i)
    for v in g.return_node(i).edges:
        if v not in visited:
            temp = dfs(v, temp, g, visited)
    return temp


def find_connected_components(g): # find connected components in a graph g
    cc = []
    visited = []
    for i in range(len(g.vertices)):
        if i not in visited:
            temp = []
            a = dfs(i, temp, g, visited)
            if len (a) > 1:
                cc.append(a)
    return cc


def return_all_edges_in_cc(c, g): #this function return all edges in graph g which is in connectec components C
    # c is a list of all vertices id in this connected component
    component_edges = []
    paired = itertools.combinations(c, 2) #all combination two of these vertices will check if the edges exist or not
    for p in paired:
        e_list = g.return_edges(p[0], p[1])
        if e_list != None:
            for e in e_list:
                component_edges.append(e)
    for i in c:
        e_list = g.return_edges(i, i)
        if e_list != None:
            for e in e_list:
                component_edges.append(e)
    return component_edges

def calculate_seg_length(e): #calculate segment length of segment edge e
    l = abs (g.return_node(e[0]).pos - g.return_node(e[1]).pos)
    return l , math.ceil(l / 5000000)


def estimating_edge_multiplicities_in_CC(component):
    Lp_prob = p.LpProblem('Problem', p.LpMinimize) # this line initiate the ILP problem
    objective = 0 # variable for some of our objective
    sv_sum = 0 # sum all variables for SV edges
    cn_tune = 0 # sum all variables for changing CN (tuning CN)
    component_edges = return_all_edges_in_cc(component, g)
    # print('Siavash', component_edges)
    for i in range(len(component_edges)):
        e = component_edges[i]
        if e[3] != 'S': #in it is not Segment edge
            component_edges[i] = [e, p.LpVariable('X' + str(i), lowBound=0, cat=p.LpInteger)] #create an ILP variable Xi for >= 0
            print('X' + str(i), e)
            if e[3] == 'SV':
                if g.return_node(e[0]).chromosome != g.return_node(e[1]).chromosome: # I want to give more weight to the translocation SVs
                    Lp_prob+= component_edges[i][1] >=0 #sum in the LP_probe
                    sv_sum += 8 * component_edges[i][1] #updating sv_sum with weight of 8 of this variable
                else:    
                    Lp_prob+= component_edges[i][1] >=0
                    sv_sum += component_edges[i][1] 
            elif e[3] == 'R': #if it is reference edge we want to have another constraint if a map overlap it it should be traverse at least one so updating the constraints
                node = g.return_node(e[0])
                if detect_overlap_map(node.chromosome, node.pos):
                    Lp_prob+= component_edges[i][1] >=1
        else: #if it is Segment edge
            component_edges[i] = [e, p.LpVariable('Y' + str(i), cat=p.LpInteger),p.LpVariable('Z' + str(i), cat=p.LpInteger)] # variable Y_i is for tuning CN of each segment if we need to change the CN
            print('Y' + str(i), e)
            # For having abselute value in ILP for we need to define variable Z as well. 
            Lp_prob += component_edges[i][2]>= -component_edges[i][1] #this is used for defining abselute value
            Lp_prob += component_edges[i][2]>= component_edges[i][1]
            # the following if is a condition for setting the limit on how much the CN can be changed for segment with length less than 1Mbp it is one
            if calculate_seg_length(e)[0] <= 1000000:
                Lp_prob+= component_edges[i][2] <=1
            else: # for rest it is 25% of their length
                Lp_prob += component_edges[i][2]<= math.ceil(e[2]/4)
    for v in component:
        v_edges = [] #all edges conntigin to node v
        for e in component_edges:
            if e[0][0] == v or e[0][1] == v:
                v_edges.append(e)
        if len(v_edges) > 1:
            cond = 0 # this is a CN balance condition 
            for e in v_edges:
                if e[0][3] == 'SV' and e[0][0] == e[0][1]:
                    cond = cond + 2 * e[1] # it there is a SV edge from node v to v
                    objective = objective - 2 * e[1]
                elif e[0][3] != 'S':
                    cond = cond + e[1] #sum all non Segment edge
                    objective = objective - e[1]
            for e in v_edges:
                if e[0][3] == 'S':
                    cond = cond + e[1] #summing tuning variable
                    Lp_prob += cond <= e[0][2] # this is important line wich we apply copy number balance condition
                    if calculate_seg_length(e[0])[0] > 50000: # for segment less than 50 Kbp no penalty applied for tuning segment CN
                        cn_tune += e[2] * calculate_seg_length(e[0])[1]
                    objective = objective + e[0][2] - e[1] #updating the Objective Function
        else:
            objective = objective + v_edges[0][0][2]
            if calculate_seg_length(v_edges[0][0])[0] > 50000:# for segment less than 50 Kbp no penalty applied for tuning segment CN
                cn_tune += v_edges[0][2] * calculate_seg_length(v_edges[0][0])[1]
    #Just for debug
    # print('obj', objective)
    # print('Sv_sum', sv_sum)
    # objective = 10 * objective  - 9 * sv_sum + 15 * cn_tune
    objective = 20 * objective  - 9 * sv_sum + 4 * cn_tune
    # print('obj', objective)
    Lp_prob += objective
    print(Lp_prob)
    status = Lp_prob.solve()
    print(p.LpStatus[status])  # The solution status
    print(Lp_prob.variables())
    for i in Lp_prob.variables():
        #Set edges multiplicity based on variiable values
        if i.name != '__dummy' and i.name.startswith('X'):
            index = int(i.name[1:])
            component_edges[index][0] = list(component_edges[index][0])
            component_edges[index][0][2] = int(i.varValue)
            component_edges[index]= tuple(component_edges[index][0])
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
        elif i[2] == 0: #if estimate of an edge is 0, update the graph  if two different type of edge connect same node it need to be check 
            #CHECK CHECK CHECK
            if g.return_edges(i[0], i[1]) == None:
                g.return_node(i[0]).remove_edges(i[1])
                if i[0]!=i[1]:
                    g.return_node(i[1]).remove_edges(i[0])
    return ans


def remove_edge(g2, e): #this function remove edge e from graph g2 and if edge e is not in g print is not in graph
    g = deepcopy(g2)
    component_edges = g.edges
    for i in range(len(component_edges)):
        if e == component_edges[i]:
            if component_edges[i][2] == 1:
                component_edges = component_edges[:i] + component_edges[i + 1:]
                if len(g.return_edges(e[0], e[1])) == 1: # if only one type of edge connect them to each other. 
                    g.return_node(e[0]).remove_edges(e[1])
                    if e[0]!=e[1]:
                        g.return_node(e[1]).remove_edges(e[0])
            else:
                component_edges[i] = list(component_edges[i])
                component_edges[i][2] -= 1
                component_edges[i] = tuple(component_edges[i])
            g.edges = component_edges
            return g
    print('Edge is not in graph', e)


def dfs_count(g, v, visited): #count number of available visired edge from v in graph g 
    #for more info look at here https://www.geeksforgeeks.org/fleurys-algorithm-for-printing-eulerian-path/
    count = 1
    visited.append(v)
    for i in g.return_node(v).edges:
        if i not in visited:
            count += dfs_count(g, i, visited)
    return count


def isValidNextEdge(g, u, v, e): # check this edge is bridge or not. Can we remove the edge?
    if len(g.return_node(u).edges) == 1 and g.return_node(u).edges[0] == v: #if there is no other way just traverse it
        return True
    else:
        #check if it is bridge or not
        visited = []
        count1 = dfs_count(g, u, visited)
        g2 = remove_edge(g, e)
        visited = []
        count2 = dfs_count(g2, u, visited)
        return False if count1 > count2 else True

def check_traverse_segment(prev,u, next): # this function check that after a any SV or R edge we traverse Segment edge #this one prevent none meaningfull path
    if prev//2 == u//2:
        return True
    elif u//2 == next//2:
        return True
    return False

def printEulerUtil(g, u, prev): #find Eulerian path or circuts in graph g starting with node u and previouse seen node is prev.
    valid = []
    next_find = False
    next_node = -1
    for v in g.return_node(u).edges:
        e_list = g.return_edges(u, v) #all edges between u and v
        for e in e_list:
            if isValidNextEdge(g, u, v, e) and check_traverse_segment(prev,u, v): #check traversing edge between u and v is valid and meaningfull
                valid.append((v,e))
                if abs(v- u) == 1 and v !=prev: #We want to force if next node is availanle traverse it at first between all options. 
                    next_find = True
                    # print('NextFind', u, v)
                    next_node = v
    valid = sorted(valid, key=lambda tup: (tup[0],tup[1][3])) # this is all valid options for traverse 
    # print('haa', u, prev, valid,next_find, next_node)
    if len(valid) == 1: #if length is equal to one Just traverse it. 
        g = remove_edge(g, valid[0][1])
        return [u] + printEulerUtil(g, valid[0][0], u)
    elif len(valid) > 1:
        if next_find: #if next is available just traverse it 
            for i in valid: 
                if i[0] == next_node:
                    g = remove_edge(g, i[1])    
                    # print('ali')
                    return [u] + printEulerUtil(g, i[0], u)
        find = False
        for i in valid: 
            if i[0]!= prev and i[0] != u: #first if next node is not previouse and current one jump to it
                g = remove_edge(g, i[1])
                find = True
                # print('reza')
                return [u] + printEulerUtil(g, i[0], u)
        if not find :
            for i in valid: 
                if i[0]!= prev:
                    g = remove_edge(g, i[1])
                    # print('mamad')
                    return [u] + printEulerUtil(g, i[0], u)
            g = remove_edge(g, valid[0][1])
            return [u] + printEulerUtil(g, valid[0][0], u)
    return [u]


def detect_segment_vertices(component, edges): # this function return those nodes which are telomere part (start and end of chromosome or no any edge like R or SV edge connect to them)
    v = set()
    for e in edges:
        if e[3] != 'S':
            v.add(e[0])
            v.add(e[1])
    return sorted(list(set(component) - v))

def detect_segment_odd_degree(component, component_edges): # detect vertices with odd degree
    ans = []
    d = {}
    for c in component:
        d[c] = 0
    for e in component_edges:
        d[e[0]] += e[2]
        d[e[1]] += e[2]
    for i in d.keys():
        if d[i]%2 != 0:
            ans.append(i)
    return ans

def scoring_paths(path_list, segment_vertices):
    best_score = 99999
    best_path = ''
    for p in path_list:
        ans = []
        temp = [p[0]]
        score = 0 
        for i in range(1,len(p)-1):
            if p[i] in segment_vertices and p[i-1]==p[i+1]:
                temp.append(p[i])
                score = score + 10 * abs(check_non_centromeric_path(temp) - 1)
                ans.append(temp)
                temp = [p[i]]
            else:
                temp.append(p[i])
        temp.append(p[-1])
        score = score + 10 * abs(check_non_centromeric_path(temp) - 1)
        ans.append(temp)
        print('paths_score', score, ans)
        if score < best_score:
            best_score = score
            best_path = p
    return best_path

def printEulerTour(component, component_edges, g): #Find Eulerian path/circuts in connected components in graph g
    g2 = Graph() #create new graph g2
    g2.edges = component_edges
    for v in component:
        g2.vertices.append(g.return_node(v))
    odd_vertices = []
    print('Component', component)
    print('Component edges', component_edges)

    segment_vertices = detect_segment_vertices(component, component_edges)
    odd_vertices = detect_segment_odd_degree(component, component_edges)
    print('ODD vertices', odd_vertices)
    # for i in range(0,len(segment_vertices), 2):
    #     g2.add_dummy_edges(segment_vertices[i], segment_vertices[i+1])
    # print(g2.edges)
    print('TOUR')
    print(segment_vertices)
    if len(odd_vertices) == 0: # Eulerian circuites exist 
        # a = printEulerUtil(g2, segment_vertices[0], -1)
        a = []
        for i in segment_vertices:
            a.append(printEulerUtil(g2, i, -1))
        a = scoring_paths(a,segment_vertices)
    elif len(odd_vertices) == 2: # Eulerian path exists
        if odd_vertices[0] in segment_vertices: #it is better to start path finding from telomere regions
            a = printEulerUtil(g2, odd_vertices[0], -1)
        elif odd_vertices[1] in segment_vertices:
            a = printEulerUtil(g2, odd_vertices[1], -1)
        else:#If not exist add a dummy edges to make it Eulerian and then find from segment vertices
            g2.add_dummy_edges(odd_vertices[0],odd_vertices[1])
            # a = printEulerUtil(g2, segment_vertices[0], -1)
            a = []
            for i in segment_vertices:
                a.append(printEulerUtil(g2, i, -1))
            a = scoring_paths(a,segment_vertices)
    else: # if more than two vertices with odd degree. Connect them to each other to make the graph Eulerian 
        if len(set(odd_vertices).intersection(set(segment_vertices))) == 0: #no telomere nodes with odd degree connect all of them to gether like previouse setp
            for i in range(0,len(odd_vertices),2):
                g2.add_dummy_edges(odd_vertices[i],odd_vertices[i+1])
            a = []
            for i in segment_vertices:
                a.append(printEulerUtil(g2, i, -1))
            a = scoring_paths(a,segment_vertices)
            # a = printEulerUtil(g2, segment_vertices[-1], -1)#Siavash in chaneed shode 
        else:
            count = 0
            save_index = 0
            for i in range(0,len(odd_vertices),2): #expect the segment vertices connect all to gether and then find Eulerian path from it
                if (odd_vertices[i] in segment_vertices or odd_vertices[i+1] in segment_vertices) and count ==0:
                    count +=1
                    save_index = i
                else:
                    g2.add_dummy_edges(odd_vertices[i],odd_vertices[i+1])
            if odd_vertices[save_index] in segment_vertices:
                a = printEulerUtil(g2, odd_vertices[save_index], -1)
            else:
                a = printEulerUtil(g2, odd_vertices[save_index+1], -1)   

    print(g2.edges)
    g2.print_node()
    print('Answer',a)
    return a 

def detect_del_dup_cn(chromosome, start, end): # this function detect that for a deletion or duplication, do we have CNV call as well or no
    #chromosome, start, end position of deletion call are input
    for i,s in enumerate(segments):#search in segments
        if int(s.chromosome) == int(chromosome):
            overlap = max(0,min(end, s.end) -max(start, s.start))
            #if these two windows are overlapping each other 0.9 and 1.2 of theirlength and the segment have CN different from one of it adjacent segments it return true.
            if overlap > 0.9 * (s.end - s.start):
                if (end - start)< 1.2 * (s.end - s.start):
                    if (s.int_cn != segments[(i+1)%len(segments)].int_cn and s.chromosome == segments[(i+1)%len(segments)].chromosome) or (s.int_cn != segments[i-1].int_cn and s.chromosome == segments[i-1].chromosome):
                        # print('Kir', chromosome, start, end, overlap)
                        # print('Kir', start , end, s.start,s.end)
                        return True, s.start , s.end
            if overlap > 0.9 * (end - start):
                if (s.end - s.start)< 1.2 * (end - start):
                    if (s.int_cn != segments[(i+1)%len(segments)].int_cn and s.chromosome == segments[(i+1)%len(segments)].chromosome) or (s.int_cn != segments[i-1].int_cn and s.chromosome == segments[i-1].chromosome):
                        return True, s.start , s.end
    return False, None, None

def detect_duplicatioon_inversion_cn(sv): # same as above for duplication calls. 
    window_lim = 50000
    node_dir1 , node_dir2 = detect_sv_directions(sv)
    if node_dir1 == 'T' and node_dir2 == 'T': #right fold back
        for i,s in enumerate(segments):
            if int(s.chromosome) == int(sv.ref_c_id1): #only compared with the next contigs
                if abs(s.end - sv.ref_end) < window_lim and (s.int_cn != segments[(i+1)%len(segments)].int_cn and s.chromosome == segments[(i+1)%len(segments)].chromosome): 
                    return True , s.end
    elif node_dir1 == 'H' and node_dir2 == 'H': # left foldback
        for i,s in enumerate(segments): #compared with prev contigs
            if int(s.chromosome) == int(sv.ref_c_id1):
                if abs(s.start - sv.ref_start) < window_lim and (s.int_cn != segments[i-1].int_cn and s.chromosome == segments[i-1].chromosome): 
                    return True , s.start
    return False, None


def detect_receprical_translocation(sv): #sometimes one of these reciprocal translocation have low confidence but this function we retrive it 
    window_lim = 200000
    for i in smap:
        if i.ref_c_id1 == sv.ref_c_id1 and i.ref_c_id2 == sv.ref_c_id2 and sv.smap_id != i.smap_id:
            if abs(i.ref_start - sv.ref_start) < window_lim and abs(i.ref_end - sv.ref_end) < window_lim:
                i_dir1, i_dir2 = detect_sv_directions(i)
                sv_dir1, sv_dir2 = detect_sv_directions(sv)
                if i_dir1 != sv_dir1 and i_dir2 != sv_dir2:
                    return True, i
    return False, None

def check_non_centromeric_path(p):
    count = 0
    for i in range(0,len(p)-1,2):
        u = g.return_node(p[i])
        v = g.return_node(p[i+1])
        if u.chromosome == v.chromosome:
            if (min(u.pos,v.pos)< min(centro['chr'+str(u.chromosome)]) and max(u.pos, v.pos)> min(centro['chr'+str(u.chromosome)])) or (min(u.pos,v.pos)< max(centro['chr'+str(u.chromosome)]) and max(u.pos, v.pos)> max(centro['chr'+str(u.chromosome)])):
                count += 1
    return count

def convert_path_to_segment(p,g): # this is important function that convert Eulerion path with vertices ID to segment path. 
    component = list(set(p))
    component_edges = return_all_edges_in_cc(component, g)
    segment_vertices = detect_segment_vertices(component, component_edges)
    ans = []
    temp = [p[0]]
    for i in range(1,len(p)-1):
        if p[i] in segment_vertices:# and p[i-1]==p[i+1]:
            temp.append(p[i])
            ans.append(temp)
            if p[i-1] == p[i+1]:
                temp = [p[i]]
            else:
                temp = []
        else:
            temp.append(p[i])
    temp.append(p[-1])
    ans.append(temp)
    ans2 = []
    for p in ans:
        temp = ''
        for i in range(0,len(p)-1,2):
            # print('as', i, len(p), p)
            seg_number = int((max(p[i],p[i+1])+1)/2)
            direction = '+'
            if p[i] > p[i+1]:
                direction = '-'
            temp = temp + str(seg_number)+direction+' '
        ans2.append(temp)

    return ans2

def check_exiest_call(chromosome , start , end , type): # if we have a call like gain or loss  but in CNV it is filtered retrive it
    for s in all_seg:
        if s.chromosome == chromosome and s.start >= start and s.end <= end and s.type[:4] == type[:4]:
            return True
    return False


def extend_segments_cn(segments):#this function check that if between two cnv call gap is less than 400Kbp and there is a call in CNV call but it marked or filtered we assume it is true and extend the segment length
    start = True
    for i in range(0 , len(segments)-1):
        s = segments[i]
        next_seg = segments[i+1]    
        if s.chromosome == next_seg.chromosome:
            if abs(next_seg.start - s.end) < 400000:
                if check_exiest_call(s.chromosome, s.end, next_seg.start, s.type):
                    s.end = next_seg.start -1
                    s.bp = [s.start , s.end]
                    segments[i] = s
    return segments
                
######################################################################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("-cnv", "--cnv", help="path to cnv call (cnv_call_exp.txt)", required=True)
parser.add_argument("-smap", "--smap", help="path to smap file", required=True)
parser.add_argument("-rcmap", "--rcmap", help="path to CNV rcmap file (cnv_rcmap_exp.txt)", required=True)
parser.add_argument("-xmap", "--xmap", help="path to contig alignments file xmap", required=True)
parser.add_argument("-centro", "--centro", help="path to file contains centromere coordinates", required=False)
parser.add_argument("-n", "--name", help="output name", required=True)
parser.add_argument("-o", "--output", help="path to output dir", required=True)
args = parser.parse_args()
segments, all_seg = parse_cnvcall(args.cnv)
smap = parse_smap(args.smap) 
rcov, rcop = parse_rcmap(args.rcmap)
chrY_cn = int(np.average(list(rcop['24'].values())) + 0.5)
chrX_cn = 2
if chrY_cn > 0:
    chrX_cn = 1
    # chrY_cn = 1
xmap = parse_xmap(args.xmap)
if args.centro is not None: # this will parse centromere region. It can be hard coded. 
    centro = parse_centro(args.centro)
else:
    centro = None
output = args.output+'/'+ args.name + '.txt'
file = args.output+'/'+ args.name + '.png'
name = args.name
svs = []
# for s in segments:
#         print(s.chromosome, s.start, s.end, s.int_cn)
segments = extend_segments_cn(segments) #fill the gap between calls. 
for k in rcop.keys():
    seg_list = []
    label_list = list(rcop[k].keys())
    for s in segments:
        if s.chromosome == k:
            if s.width > 200000: #if call has length greater than 200Kbp assume a segment
                seg_list.append(s)
    prev_point = list(rcop[k].keys())[0]
    if len(seg_list) == 0: #create segment for start of chromosme
        new_seg = Segments()
        new_seg.start = 0
        new_seg.end = list(rcop[k].keys())[-1]
        new_seg.width = list(rcop[k].keys())[-1]
        new_seg.chromosome = k
        new_seg.fractional_cn = 2
        new_seg.int_cn = 2 #assumption that sample is diploide. default CN = 2
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
            start, end = find_start_end(prev_point, s.start, label_list) # there are labels between two segments then create segment with CN =2 between them
            if start != 0 and end != 0: #
                new_seg = Segments()
                new_seg.start = start
                new_seg.end = end
                new_seg.width = end - start
                new_seg.chromosome = k
                new_seg.fractional_cn = 2
                new_seg.int_cn = 2 #assumption that sample is diploide. default CN = 2
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
            new_seg.fractional_cn = 2
            new_seg.int_cn = 2
            if int(k) == 23:
                new_seg.fractional_cn = chrX_cn
                new_seg.int_cn = chrX_cn
            if int(k) == 24:
                new_seg.fractional_cn = chrY_cn
                new_seg.int_cn = chrY_cn
            new_seg.bp = [start, end]
            segments.append(new_seg)

segments.sort(key=lambda x: (int(x.chromosome), x.start))
# for s in segments:
#     print('asli', s.chromosome, s.start, s.end, s.int_cn, sorted(s.bp))
for i in smap:
    # translocation applied filters. 
    if i.sv_type.startswith('trans') and i.confidence >= 0.05 and not i.sv_type.endswith(
            'segdupe') and not i.sv_type.endswith('common') and not i.sv_type.endswith('oveerlap'):
        svs.append(i)
        exist,s = detect_receprical_translocation(i)
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
                if detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end)[0]:
                    _, i.ref_start, i.ref_end = detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end)
                    svs.append(i)
                    print(i.line.strip())
            elif i.size > 500000 and abs(i.ref_end - i.ref_start) > 500000: # this would be for insertion length more than 500Kbp
                svs.append(i)
                print(i.line.strip())
    #if we have inversion SV
    elif i.sv_type == 'inversion' and i.confidence >= 0.7: # filter low confidance 
        start, end = 0, 0
        dir = ''
        dir1, dir2 = detect_sv_directions(i)
        s = find_in_smap(i.linkID) #inversion has two rwo in smap file. we find them with Link ID 
        if dir1 == 'H': #update inversion call. it is to complicated but baisically calculate inversion start and endpoint
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
        if abs(end - start) > 800000: #apply filter on size of inversion 
            svs.append(i)
            print(i.line.strip(), start, end,dir)
            print(s.line.strip())
    elif i.sv_type == 'inversion_paired' and i.confidence >= 0.7: #if it is full inversion
        s = find_in_smap(i.linkID)
        if abs(i.ref_start - s.ref_end) > 500000:
            i.ref_end = s.ref_end
            print(i.line.strip())
            svs.append(i)
    elif i.sv_type == 'duplication_inverted':
        if detect_duplicatioon_inversion_cn(i)[0]:
            _, fold_point = detect_duplicatioon_inversion_cn(i) #because CN is changed
            i.ref_start = fold_point
            i.ref_end = fold_point
            svs.append(i)
            print(i.line.strip())
    elif i.sv_type == 'duplication' or i.sv_type == 'duplication_split':
        if detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end)[0]:
            _, i.ref_start, i.ref_end = detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end)
            svs.append(i)
            print(i.line.strip())
    elif i.size > 500000 and not i.sv_type.startswith('inversion'): #Other type of SV
        svs.append(i)
        print(i.line.strip())

for sv in svs:#integrate BPs and Segments
    find_bp_in_segment(sv.ref_c_id1, sv.ref_start) #
    find_bp_in_segment(sv.ref_c_id2, sv.ref_end)
#merging Bps for spiliting segments 
for s in segments:
    s.bp = merge_list(s.bp)

#Graph Creation
aa = 0
g = Graph()
counter = 0
prev_chr = 0
for index, s in enumerate(segments): #Forach segment creates two vertices
    for inedx_bp, i in enumerate(s.bp):
        if inedx_bp == 0: #if it is first BP in a segment
            aa += 1
            v = Vertices() #Create Vetrices Object
            v.chromosome = s.chromosome
            v.id = counter
            v.pos = i
            v.cn = s.int_cn
            v.type = 'H' # the start node is Head node
            g.vertices.append(v)
            if prev_chr == s.chromosome:
                g.edges.append((counter - 1, counter, 0, 'R'))  #addign reference edge between two continues segments
                g.return_node(counter - 1).append_edges(counter) #update adjancency matrix
                g.return_node(counter).append_edges(counter - 1)
            prev_chr = s.chromosome
            counter += 1
        elif inedx_bp == len(s.bp) - 1: #if it is last BP in a segment
            aa += 1
            v = Vertices()
            v.chromosome = s.chromosome
            v.id = counter
            v.pos = i
            v.cn = s.int_cn
            v.type = 'T' #last node should be Tail node
            g.vertices.append(v)
            g.edges.append((counter - 1, counter, s.int_cn, 'S'))
            g.return_node(counter).append_edges(counter - 1)
            g.return_node(counter - 1).append_edges(counter)
            # g.edges.append((counter,counter-1, s.int_cn,'S'))
            counter += 1
            prev_chr = s.chromosome
        else: #Create both Tail and Head node
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
    n_type1, n_type2 = detect_sv_directions(sv)
    #a and b are two nodes that ara connected by sv edge
    a = find_nodes(sv.ref_c_id1, sv.ref_start, g.vertices, n_type1)
    b = find_nodes(sv.ref_c_id2, sv.ref_end, g.vertices, n_type2)
    # if sv.sv_type == 'inversion' or sv.sv_type == 'inversion_paired':
    if sv.sv_type == 'inversion_paired': #Lets complete the inversion
        if g.return_node(a).type == 'H':
            new_edge = (a-1,b-1,0,'SV')
        else:
            new_edge = (a+1,b+1,0,'SV')
        if new_edge not in g.edges:
                g.edges.append(new_edge)
                g.return_node(new_edge[0]).append_edges(new_edge[1])
                g.return_node(new_edge[1]).append_edges(new_edge[0])
    if b == a and (a, b, 0, 'SV') not in g.edges: #it can be happend in duplication inversion
        g.edges.append((a, b, 0, 'SV'))
        g.return_node(a).append_edges(b)
    #     print(sv.line)
    #     print(a,b)
    elif (a, b, 0, 'SV') not in g.edges:
        g.edges.append((a, b, 0, 'SV'))
        g.return_node(a).append_edges(b)
        g.return_node(b).append_edges(a)
    # g.edges.append((b,a,0,'SV'))
g.print_node()
print(g.edges)
Plot_graph(g,file,name)
connected_components = find_connected_components(g)
for component in connected_components:
    # if 148 in component:
        component_edges = estimating_edge_multiplicities_in_CC(component)
connected_components = find_connected_components(g)
paths = []
for component in connected_components:
    # if 148 in component:
        component_edges = return_all_edges_in_cc(component, g)
        print(component)
        print(component_edges)
        paths.append(printEulerTour(component, component_edges, g))
# print(paths)
#write in the output
with open(output , 'w') as f :
    f.write('Segment\tNumber\tChromosome\tStart\tEnd\tStartNode\tEndNode\n')
    number = 1
    for i in range(0,len(g.vertices),2):
        v = g.vertices[i]
        u = g.vertices[i+1]
        f.write('Segment\t{id}\t{chrom}\t{start}\t{end}\t{snode}\t{enode}\n'.format(id = number, chrom = v.chromosome, start= v.pos, end= u.pos, snode = v.id, enode = u.id))
        number += 1
    c = 1
    for p in paths:
        for structure in convert_path_to_segment(p,g):
            # f.write('Path'+str(c)+'='+','.join(str(z) for z in p)+'\n')
            f.write('Path'+str(c)+ ' = '+structure+'\n')
            c+=1
    
