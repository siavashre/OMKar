import itertools
import math
######################### DOUBLE CHECK THIS FUNCTION ######################################
def detect_sv_directions(sv, xmap):  # this function detect the direction of one smap like H/T to H/T. #need to be check Contain bug?
    if sv.sv_type == 'inversion_paired':  # for inversion_paired always return T to T
        return 'T', 'T'
    elif sv.sv_type == 'duplication':  # for duplication always return T to T
        return 'H', 'T'
    dir1, dir2 = '', ''
    xmap_id1 = sv.xmap_id1
    xmap_id2 = sv.xmap_id2
    alignment1 = xmap[str(xmap_id1)]
    alignment2 = xmap[str(xmap_id2)]
    swap = 0
    if min(alignment1['QryStartPos'], alignment1['QryEndPos']) > min(alignment2['QryStartPos'],
                                                                     alignment2['QryEndPos']):
        swap = 1
        alignment1, alignment2 = alignment2, alignment1  # sort alignment1 and alignment2 by query position. Always alignment1 has smaller position on the query
    if alignment1['Orientation'] == '+':
        dir1 = 'T'
    else:
        dir1 = 'H'
    if alignment2['Orientation'] == '+':
        dir2 = 'H'
    else:
        dir2 = 'T'
    if int(alignment1['RefContigID']) > int(alignment2['RefContigID']) or (int(alignment1['RefContigID']) == int(alignment2['RefContigID'])
                                                                           and alignment1['RefStartPos'] > alignment2['RefStartPos']):
        if dir1 == 'T' and dir2 == 'H':
            return 'H', 'T'
        elif dir1 == 'H' and dir2 == 'T':
            return 'T', 'H'
    return dir1, dir2

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



def find_nodes(chromosome, pos, vertices, node_type):  # find the closest node and return it id to the chromosome and position and type
    dist = 999999999  # as input take the chromosome number and position and type of node we are looking to it(H or T) and return it.
    ans = -1
    for v in vertices:
        if v.chromosome == chromosome and v.type == node_type:
            if abs(v.pos - pos) < dist:
                dist = abs(v.pos - pos)
                ans = v.id
    return ans

def find_in_smap(id, smap):  # return a smap call with the iD otherwise return None
    for s in smap:
        if s.smap_id == id:
            return s
    return None

def is_overlapping(start1, end1, start2, end2):
    return start1 <= end2 and start2 <= end1
def is_overlapping_half_centro(start1, end1, start2, end2): # This function check if a segment overlap at least half of centromere region we count it as 1 otherwise 0
    #Centro start --> start2
    #centro end --> end2
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap_percent = (overlap_end - overlap_start) / (end2 - start2)
    return max(0, overlap_percent)
    # if overlap_percent >= 0.5:
    #     return True
    # return False
    # return start1 <= end2 and start2 <= end1

def check_non_centromeric_path(p, g, centro):
    count = 0
    for i in range(0, len(p) - 1, 2):
        u = g.return_node(p[i])
        v = g.return_node(p[i + 1])
        if u.chromosome == v.chromosome:
            # if (min(u.pos,v.pos)< min(centro['chr'+str(u.chromosome)]) and max(u.pos, v.pos)> min(centro['chr'+str(u.chromosome)])) or (min(u.pos,v.pos)< max(centro['chr'+str(u.chromosome)]) and max(u.pos, v.pos)> max(centro['chr'+str(u.chromosome)])):
            # if is_overlapping_half_centro(min(u.pos, v.pos), max(u.pos, v.pos), min(centro['chr' + str(u.chromosome)]), max(centro['chr' + str(u.chromosome)])):
            #     count += 1
            count += is_overlapping_half_centro(min(u.pos, v.pos), max(u.pos, v.pos), min(centro['chr' + str(u.chromosome)]), max(centro['chr' + str(u.chromosome)]))

    return math.ceil(count)


def detect_segment_vertices(component,
                            edges):  # this function return those nodes which are telomere part (start and end of chromosome or no any edge like R or SV edge connect to them)
    v = set()
    for e in edges:
        if e[3] == 'R':
            v.add(e[0])
            v.add(e[1])
    return sorted(list(set(component) - v))

def return_all_edges_in_cc(c, g):  # this function return all edges in graph g which is in connectec components C
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


def calculate_seg_length(e, g):  # calculate segment length of segment edge e
    l = abs(g.return_node(e[0]).pos - g.return_node(e[1]).pos)
    return l, 1 + math.ceil(l / 5000000)

##########################ZHAOYANG#################
import numpy as np


def sign(x):
    if x < 0:
        return -1
    elif x == 0:
        return 0
    else:
        return 1


def reverse_dict(input_dict):
    """
    requires mapping to be bijective
    :param input_dict:
    :return:
    """
    output_dict = {}
    for position, name in input_dict.items():
        output_dict[name] = position
    return output_dict


def generate_parabola(start_x, end_x, peak_x, peak_y, num_points=100):
    # Solve for the coefficients of the parabola
    A = np.array([[start_x ** 2, start_x, 1],
                  [peak_x ** 2, peak_x, 1],
                  [end_x ** 2, end_x, 1]])
    b = np.array([0.5, peak_y, 0.5])
    coef = np.linalg.solve(A, b)

    # Generate x values
    x = np.linspace(start_x, end_x, num_points)
    # Calculate y values for the parabola
    y = coef[0] * x ** 2 + coef[1] * x + coef[2]
    # Return x, y as tuples with precision to the sixth decimal point
    return np.array([(round(x_val, 6), round(y_val, 6)) for x_val, y_val in zip(x, y)])


def generate_circle(start_x, peak_y, circle_size_multiplier, num_points=100):
    diameter = abs(peak_y - 0.5) * circle_size_multiplier
    radius = diameter / 2
    is_upward_direction = (peak_y > 0.5)

    # Calculate the center of the circle
    center_x = start_x
    if is_upward_direction:
        center_y = 0.5 + radius  # Fixed y-coordinate
        # Generate angles for the circle (starting from the bottom and going clockwise)
        angles = np.linspace(3 * np.pi / 2, -np.pi / 2, num_points)
    else:
        center_y = 0.5 - radius
        # Generate angles for the circle (starting from the top and going counterclockwise)
        angles = np.linspace(np.pi / 2, -3 * np.pi / 2, num_points)

    # Generate x and y coordinates for the circle
    x = center_x + radius * np.cos(angles)
    y = center_y + radius * np.sin(angles)

    # Combine x and y coordinates into points
    points = np.array([[round(x_val, 6), round(y_val, 6)] for x_val, y_val in zip(x, y)])

    # Make sure the first and last points are equal to the starting point
    if points[0][0] != start_x or points[-1][0] != start_x:
        raise RuntimeError('circle coordinate formation error')

    return points
