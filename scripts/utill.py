import math
######################### DOUBLE CHECK THIS FUNCTION ######################################
def detect_sv_directions(sv, xmap):  # this function detect the direction of one smap like H/T to H/T. #need to be check Contain bug?
    """
    Detects the directions (Head/Tail) for structural variations (SVs).

    Args:
        sv: An object representing a structural variation.
        xmap (dict): Dictionary containing Xmap alignment data.

    Returns:
        tuple: Directions as ('H', 'T'), ('T', 'H'), or other combinations based on SV type and orientation.
    """
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
    """
    Merges all breakpoints within a 50Kbp window in a genomic segment.

    Args:
        l (list): List of breakpoint positions.

    Returns:
        list: List of merged breakpoint positions.
    """
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
    """
    Reverses a direction (Head/Tail).

    Args:
        a (str): Direction ('H' or 'T').

    Returns:
        str: Reversed direction ('T' or 'H').
    """
    if a == 'H':
        return 'T'
    return 'H'



def find_nodes(chromosome, pos, vertices, node_type):  # find the closest node and return it id to the chromosome and position and type
    """
    Finds the closest node by chromosome, position, and type.

    Args:
        chromosome (int): Chromosome number.
        pos (int): Position in base pairs.
        vertices (list): List of vertex objects.
        node_type (str): Type of the node ('H' or 'T').

    Returns:
        int: ID of the closest node, or -1 if none is found.
    """
    dist = 999999999  # as input take the chromosome number and position and type of node we are looking to it(H or T) and return it.
    ans = -1
    for v in vertices:
        if v.chromosome == chromosome and v.type == node_type:
            if abs(v.pos - pos) < dist:
                dist = abs(v.pos - pos)
                ans = v.id
    return ans

def find_in_smap(id, smap):  # return a smap call with the iD otherwise return None
    """
    Finds a structural variation in the Smap file by its ID.

    Args:
        id (int): Smap entry ID.
        smap (list): List of Smap entries.

    Returns:
        object or None: The matching Smap entry, or None if not found.
    """
    for s in smap:
        if s.smap_id == id:
            return s
    return None

def is_overlapping(start1, end1, start2, end2):
    """
    Checks if two genomic regions overlap.

    Args:
        start1 (int): Start position of the first region.
        end1 (int): End position of the first region.
        start2 (int): Start position of the second region.
        end2 (int): End position of the second region.

    Returns:
        bool: True if the regions overlap, False otherwise.
    """
    return start1 <= end2 and start2 <= end1
def is_overlapping_half_centro(start1, end1, start2, end2): # This function check if a segment overlap at least half of centromere region we count it as 1 otherwise 0
    """
    Checks if a segment overlaps at least half of the centromere region.

    Args:
        start1 (int): Start position of the segment.
        end1 (int): End position of the segment.
        start2 (int): Start position of the centromere.
        end2 (int): End position of the centromere.

    Returns:
        float: Percentage of overlap between the segment and the centromere.
    """
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap_percent = (overlap_end - overlap_start) / (end2 - start2)
    return max(0, overlap_percent)

def check_non_centromeric_path(p, g, centro):
    """
    Counts overlaps with centromeric regions in a genomic path.

    Args:
        p (list): List of vertex IDs representing the path.
        g (Graph): Graph object containing vertices and edges.
        centro (dict): Dictionary of centromere regions.

    Returns:
        int: Count of overlaps with centromeric regions.
    """
    count = 0
    for i in range(0, len(p) - 1, 2):
        u = g.return_node(p[i])
        v = g.return_node(p[i + 1])
        if u.chromosome == v.chromosome:
            count += is_overlapping_half_centro(min(u.pos, v.pos), max(u.pos, v.pos), min(centro['chr' + str(u.chromosome)]), max(centro['chr' + str(u.chromosome)]))

    return math.ceil(count)


def detect_segment_vertices(component,
                            edges):  # this function return those nodes which are telomere part (start and end of chromosome or no any edge like R or SV edge connect to them)
    """
    Detects vertices at the telomeric part of segments.

    Args:
        component (set): Set of vertex IDs.
        edges (list): List of edges.

    Returns:
        list: List of vertices that are telomeric parts.
    """
    v = set()
    for e in edges:
        if e[3] == 'R':
            v.add(e[0])
            v.add(e[1])
    return sorted(list(set(component) - v))


def calculate_seg_length(e, g):  # calculate segment length of segment edge e
    """
    Calculates the length of a genomic segment based on its edges.

    Args:
        e (tuple): Edge representing the segment.
        g (Graph): Graph object containing vertices and edges.

    Returns:
        tuple: Length of the segment and its scaled length.
    """
    l = abs(g.return_node(e[0]).pos - g.return_node(e[1]).pos)
    return l, 1 + math.ceil(l / 5000000)

def find_bp_in_segment(chromosome, point,
                       segments):  # this function get chromosome and position as input(from sv call) and add this point to segment bp list. It helps us to then spliting segments to new one.
    """
    Adds a breakpoint to the appropriate genomic segment.

    Args:
        chromosome (int): Chromosome number.
        point (int): Position of the breakpoint.
        segments (list): List of segment objects.
    """
    for i in segments:
        if str(i.chromosome) == str(chromosome):
            if i.start <= point <= i.end:
                i.bp.append(point)

def find_start_end(prev_point, start, label_list):
    """
    Finds the start and end of a segment given a list of labels.

    Args:
        prev_point (int): Previous point position.
        start (int): Start position of the current segment.
        label_list (list): List of positions.

    Returns:
        tuple: Start and end positions for the segment.
    """
    ans = []
    for i in label_list:
        if prev_point <= i < start:
            ans.append(i)
    if len(ans) < 2:
        return 0, 0
    else:
        return min(min(ans)+1, prev_point+1), start-1

def check_contains_masked_region(l_m, start, end):
    """
    Checks if a genomic region contains or overlaps with a masked region.

    Args:
        l_m (list): List of masked regions.
        start (int): Start position of the region.
        end (int): End position of the region.

    Returns:
        tuple: Flags indicating if the region contains or borders a masked region.
    """
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

def detect_overlap_map(chromosome, pos, xmap):  # this function detect if there is a map(alignment that overlap this region)
    # maybe it is good idea instead of 1bp assume a window here.
    """
    Detects whether an alignment overlaps with a genomic position.

    Args:
        chromosome (int): Chromosome number.
        pos (int): Genomic position.
        xmap (dict): Dictionary of alignment data.

    Returns:
        bool: True if an overlap exists, False otherwise.
    """
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

def average_values_greater_than_a(input_dict, a):
    """
    Calculates the average of values greater than a threshold.

    Args:
        input_dict (dict): Dictionary of values.
        a (int): Threshold value.

    Returns:
        float or None: Average of values greater than the threshold, or None if no such values exist.
    """
    total = 0
    count = 0
    for key, value in input_dict.items():
        if key > a:
            total += value
            count += 1
    return total / count if count > 0 else None


def cn_in_mask_N_region(chromosome, start, end, cop, centro):
    """
    Calculates the copy number in masked N regions of specific chromosomes.

    Args:
        chromosome (str): Chromosome number.
        start (int): Start position of the region.
        end (int): End position of the region.
        cop (dict): Dictionary of copy numbers.
        centro (dict): Dictionary of centromere positions.

    Returns:
        int: Copy number in the masked N region.
    """
    if chromosome not in ['22', '21', '15', '14', '13']:  # Andy shared this chromosome with me that Bionano CNV call in P arm is not reliable
        return round(np.average(list(cop.values()))) #instead of assumption that it is always diploid return the average cn of that chromosome
    else:
        if chromosome == '22' and start < max(centro['chr22']) and end < max(centro['chr22']):
            return round(average_values_greater_than_a(cop, max(centro['chr22'])))
        if chromosome == '21' and start < max(centro['chr21']) and end < max(centro['chr21']):
            return round(average_values_greater_than_a(cop, max(centro['chr21'])))
        if chromosome == '15' and start < max(centro['chr15']) and end < max(centro['chr15']):
            return round(average_values_greater_than_a(cop, max(centro['chr15'])))
        if chromosome == '14' and start < max(centro['chr14']) and end < max(centro['chr14']):
            return round(average_values_greater_than_a(cop, max(centro['chr14'])))
        if chromosome == '13' and start < max(centro['chr13']) and end < max(centro['chr13']):
            return round(average_values_greater_than_a(cop, max(centro['chr13'])))
    return 2

##########################ZHAOYANG#################
import numpy as np


def sign(x):
    """
    Determines the sign of a number.

    Args:
        x (float): Input number.

    Returns:
        int: -1 for negative, 1 for positive, and 0 for zero.
    """
    if x < 0:
        return -1
    elif x == 0:
        return 0
    else:
        return 1


def reverse_dict(input_dict):
    """
    Reverses a bijective dictionary.

    Args:
        input_dict (dict): Input dictionary.

    Returns:
        dict: Reversed dictionary.
    """
    output_dict = {}
    for position, name in input_dict.items():
        output_dict[name] = position
    return output_dict


def generate_parabola(start_x, end_x, peak_x, peak_y, num_points=100):
    """
    Generates a parabola between specified points.

    Args:
        start_x (float): Starting x-coordinate.
        end_x (float): Ending x-coordinate.
        peak_x (float): X-coordinate of the peak.
        peak_y (float): Y-coordinate of the peak.
        num_points (int): Number of points to generate along the parabola.

    Returns:
        np.ndarray: Array of (x, y) points forming the parabola.
    """
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
    """
    Generates a circle based on given parameters.

    Args:
        start_x (float): X-coordinate of the circle center.
        peak_y (float): Y-coordinate of the peak.
        circle_size_multiplier (float): Multiplier for the circle's size.
        num_points (int): Number of points to generate for the circle.

    Returns:
        np.ndarray: Array of (x, y) points forming the circle.
    """
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
