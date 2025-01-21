import itertools
import pulp as p
from scripts.utill import *
from copy import deepcopy
from scripts.objects import *
def return_all_edges_in_cc(c, g):  # this function return all edges in graph g which is in connectec components C
    """
    Returns all edges in a connected component of a graph.

    Args:
        c (list): List of vertex IDs in the connected component.
        g (Graph): Graph object.

    Returns:
        list: List of edges present in the connected component.
    """
    # c is a list of all vertices id in this connected component
    component_edges = []
    paired = itertools.combinations(c, 2)  # all combination two of these vertices will check if the edges exist or not
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

def dfs(i, temp, g, visited):  # run dfs alg on graph g on vertices i. visited is a list of seen vertices before visiting this node.
    """
    Performs Depth-First Search (DFS) on a graph.

    Args:
        i (int): Current vertex ID.
        temp (list): List to store the vertices visited during this DFS call.
        g (Graph): Graph object.
        visited (list): List of visited vertex IDs.

    Returns:
        list: Updated list of vertices visited during DFS.
    """
    visited.append(i)
    temp.append(i)
    for v in g.return_node(i).edges:
        if v not in visited:
            temp = dfs(v, temp, g, visited)
    return temp


def find_connected_components(g):  # find connected components in a graph g
    """
    Identifies all connected components in a graph.

    Args:
        g (Graph): Graph object.

    Returns:
        list: List of connected components, where each component is represented
              as a list of vertex IDs.
    """
    cc = []
    visited = []
    for i in range(len(g.vertices)):
        if i not in visited:
            temp = []
            a = dfs(i, temp, g, visited)
            if len(a) > 1:
                cc.append(a)
    return cc



def estimating_edge_multiplicities_in_CC(component, g, xmap):
    """
    Estimates edge multiplicities for a connected component using Integer Linear Programming (ILP).

    Args:
        component (list): List of vertex IDs in the connected component.
        g (Graph): Graph object.
        xmap (dict): Dictionary containing alignment data.

    Returns:
        list: List of edges with updated multiplicities.
    """
    Lp_prob = p.LpProblem('Problem', p.LpMinimize)  # this line initiate the ILP problem
    objective = 0  # variable for some of our objective
    sv_sum = 0  # sum all variables for SV edges
    cn_tune = 0  # sum all variables for changing CN (tuning CN)
    odd_vertices_pen = 0
    component_edges = return_all_edges_in_cc(component, g)
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
        else:
            objective = objective + v_edges[0][0][2] - v_edges[0][1]
            aa = p.LpVariable('A' + str(v), lowBound=0, cat=p.LpInteger)
            print('A' + str(v), aa)
            bb = p.LpVariable('B' + str(v), lowBound=0, cat=p.LpBinary) #This variable shows if sum dgree this vertex is odd or no
            print('B' + str(v), bb)
            Lp_prob += v_edges[0][0][2] - v_edges[0][1] - 2 * aa == bb
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
    """
    Removes an edge from a graph.

    Args:
        g2 (Graph): Graph object.
        e (tuple): Edge to be removed.

    Returns:
        Graph: Updated graph object with the edge removed.
    """
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
    """
    Counts the number of reachable vertices from a given vertex using DFS.

    Args:
        g (Graph): Graph object.
        v (int): Starting vertex ID.
        visited (list): List of visited vertex IDs.

    Returns:
        int: Number of vertices reachable from the starting vertex.
    """
    count = 1
    visited.append(v)
    for i in g.return_node(v).edges:
        if i not in visited:
            count += dfs_count(g, i, visited)
    return count


def isValidNextEdge(g, u, v, e):  # check this edge is bridge or not. Can we remove the edge?
    """
    Determines if an edge is valid for traversal.

    Args:
        g (Graph): Graph object.
        u (int): Starting vertex ID.
        v (int): Ending vertex ID.
        e (tuple): Edge to be checked.

    Returns:
        bool: True if the edge is valid, False otherwise.
    """
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
    """
    Checks if traversing a segment edge is valid.

    Args:
        prev (int): Previous vertex ID.
        u (int): Current vertex ID.
        next (int): Next vertex ID.

    Returns:
        bool: True if traversal is valid, False otherwise.
    """
    if prev // 2 == u // 2:
        return True
    elif u // 2 == next // 2:
        return True
    return False

def check_traverse_two_consecutive_segment(prev, u , next,g, chrom):
    """
    Checks if two consecutive segment traversals are valid.

    Args:
        prev (int): Previous vertex ID.
        u (int): Current vertex ID.
        next (int): Next vertex ID.
        g (Graph): Graph object.
        chrom (int): Chromosome number.

    Returns:
        bool: True if traversal is valid, False otherwise.
    """
    if prev == next and g.return_edges(u, next) != None:
        if len( g.return_edges(u, next)) ==1 :
            if g.return_edges(u, next)[0][3]=='S':
                if g.return_node(max(u , next) + 1) != None  and g.return_node(max(min(u , next) - 1, 0)) != None:
                    return False
    return  True


def printEulerUtil(g, u, prev, chrom):  # find Eulerian path or circuts in graph g starting with node u and previouse seen node is prev.
    """
    Recursively finds an Eulerian path or circuit starting from a vertex.

    Args:
        g (Graph): Graph object.
        u (int): Starting vertex ID.
        prev (int): Previous vertex ID.
        chrom (int): Chromosome number.

    Returns:
        list: List of vertex IDs forming the Eulerian path or circuit.
    """
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
    if len(valid) == 1:  # if length is equal to one Just traverse it.
        g = remove_edge(g, valid[0][1])
        return [u] + printEulerUtil(g, valid[0][0], u, chrom)
    elif len(valid) > 1:
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
                return [u] + printEulerUtil(g, i[0], u, chrom)
        if not find:
            for i in valid:
                if i[0] != prev:
                    g = remove_edge(g, i[1])
                    return [u] + printEulerUtil(g, i[0], u, chrom)
            g = remove_edge(g, valid[0][1])
            return [u] + printEulerUtil(g, valid[0][0], u, chrom)
    return [u]





def detect_segment_odd_degree(component, component_edges):  # detect vertices with odd degree
    """
    Detects vertices with odd degrees in a connected component.

    Args:
        component (list): List of vertex IDs in the connected component.
        component_edges (list): List of edges in the connected component.

    Returns:
        list: List of vertex IDs with odd degrees.
    """
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
    """
    Detects residual degrees of vertices in a connected component.

    Args:
        component (list): List of vertex IDs in the connected component.
        component_edges (list): List of edges in the connected component.
        terminal_v_ids (list): List of terminal vertex IDs.

    Returns:
        list: List of tuples (vertex ID, residual degree).
    """
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
    """
    Scores Eulerian paths based on centromeric constraints and other factors.

    Args:
        path_list (list): List of potential Eulerian paths.
        segment_vertices (list): List of vertices representing genomic segments.
        g (Graph): Graph object.
        centro (dict): Dictionary of centromere positions.

    Returns:
        str: The best-scoring Eulerian path.
    """
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


def printEulerTour(component, component_edges, g, output_file, centro):  # Find Eulerian path/circuts in connected components in graph g
    """
    Finds Eulerian paths or circuits in a connected component of a graph.

    Args:
        component (list): List of vertex IDs in the connected component.
        component_edges (list): List of edges in the connected component.
        g (Graph): Graph object.
        output_file (str): Path to save the output.
        centro (dict): Dictionary of centromere positions.

    Returns:
        tuple: (Eulerian path, list of edges with dummy edges added).
    """
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
    print(segment_vertices)
    if len(odd_vertices) == 0:  # Eulerian circuites exist
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
