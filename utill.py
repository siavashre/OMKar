import itertools
import math
import numpy as np

######################### DOUBLE CHECK THIS FUNCTION ######################################
def detect_sv_directions(sv, xmap):
    """Detect the direction of structural variations."""
    if sv.sv_type == 'inversion_paired':
        return 'T', 'T'
    elif sv.sv_type == 'duplication':
        return 'H', 'T'

    alignment1 = xmap[str(sv.xmap_id1)]
    alignment2 = xmap[str(sv.xmap_id2)]

    # Sort alignments by query position
    if min(alignment1['QryStartPos'], alignment1['QryEndPos']) > min(alignment2['QryStartPos'], alignment2['QryEndPos']):
        alignment1, alignment2 = alignment2, alignment1

    dir1 = 'T' if alignment1['Orientation'] == '+' else 'H'
    dir2 = 'H' if alignment2['Orientation'] == '+' else 'T'

    if int(alignment1['RefContigID']) > int(alignment2['RefContigID']) or (
            int(alignment1['RefContigID']) == int(alignment2['RefContigID']) and alignment1['RefStartPos'] > alignment2['RefStartPos']):
        if dir1 == 'T' and dir2 == 'H':
            return 'H', 'T'
        elif dir1 == 'H' and dir2 == 'T':
            return 'T', 'H'

    return dir1, dir2

def find_nodes(chromosome, pos, vertices, node_type):
    """Find the closest node by chromosome, position, and type."""
    dist = float('inf')
    ans = -1
    for v in vertices:
        if v.chromosome == chromosome and v.type == node_type and abs(v.pos - pos) < dist:
            dist = abs(v.pos - pos)
            ans = v.id
    return ans

def find_in_smap(id, smap):
    """Return a smap call with the given ID, otherwise return None."""
    return next((s for s in smap if s.smap_id == id), None)

def is_overlapping(start1, end1, start2, end2):
    """Check if two segments overlap."""
    return start1 <= end2 and start2 <= end1
    
def is_overlapping_half_centro(start1, end1, start2, end2):
    """Check if a segment overlaps at least half of the centromere region."""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap_percent = (overlap_end - overlap_start) / (end2 - start2)
    return max(0, overlap_percent)

def check_non_centromeric_path(p, g, centro):
    """Count the number of segments overlapping with centromere regions."""
    count = 0
    for i in range(0, len(p) - 1, 2):
        u = g.return_node(p[i])
        v = g.return_node(p[i + 1])
        if u.chromosome == v.chromosome:
            count += is_overlapping_half_centro(min(u.pos, v.pos), max(u.pos, v.pos), min(centro['chr' + str(u.chromosome)]), max(centro['chr' + str(u.chromosome)]))
    return math.ceil(count)


def detect_segment_vertices(component, edges):
    """Return nodes which are telomere parts."""
    telomere_nodes = {e[0] for e in edges if e[3] == 'R'} | {e[1] for e in edges if e[3] == 'R'}
    return sorted(set(component) - telomere_nodes)

def return_all_edges_in_cc(c, g):
    """Return all edges in a connected component."""
    component_edges = []
    paired = itertools.combinations(c, 2)
    for p in paired:
        e_list = g.return_edges(p[0], p[1])
        if e_list:
            component_edges.extend(e_list)
    for i in c:
        e_list = g.return_edges(i, i)
        if e_list:
            component_edges.extend(e_list)
    return component_edges


def sign(x):
    """Return the sign of x."""
    return -1 if x < 0 else 0 if x == 0 else 1

def reverse_dict(input_dict):
    """Reverse a bijective dictionary."""
    return {v: k for k, v in input_dict.items()}


def generate_parabola(start_x, end_x, peak_x, peak_y, num_points=100):
    """Generate points for a parabola."""
    A = np.array([[start_x ** 2, start_x, 1],
                  [peak_x ** 2, peak_x, 1],
                  [end_x ** 2, end_x, 1]])
    b = np.array([0.5, peak_y, 0.5])
    coef = np.linalg.solve(A, b)

    x = np.linspace(start_x, end_x, num_points)
    y = coef[0] * x ** 2 + coef[1] * x + coef[2]
    return np.array([(round(x_val, 6), round(y_val, 6)) for x_val, y_val in zip(x, y)])

def generate_circle(start_x, peak_y, circle_size_multiplier, num_points=100):
    """Generate points for a circle."""
    diameter = abs(peak_y - 0.5) * circle_size_multiplier
    radius = diameter / 2
    is_upward_direction = (peak_y > 0.5)

    center_x = start_x
    center_y = 0.5 + radius if is_upward_direction else 0.5 - radius
    angles = np.linspace(3 * np.pi / 2, -np.pi / 2, num_points) if is_upward_direction else np.linspace(np.pi / 2, -3 * np.pi / 2, num_points)

    x = center_x + radius * np.cos(angles)
    y = center_y + radius * np.sin(angles)

    points = np.array([[round(x_val, 6), round(y_val, 6)] for x_val, y_val in zip(x, y)])

    if points[0][0] != start_x or points[-1][0] != start_x:
        raise RuntimeError('circle coordinate formation error')

    return points
