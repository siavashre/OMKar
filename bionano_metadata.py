import pandas as pd
from scripts.utill import *
def map_cyto_coords(coords, cytoband_filtered, orientation):
    """
    Map cytogenetic coordinates based on the provided cytoband data.

    Args:
        coords (str): Coordinate string in the format 'chrom:start-end'.
        cytoband_filtered (pd.DataFrame): Filtered cytoband data.
        orientation (str): Orientation of the coordinates.

    Returns:
        str: A formatted string with chromosome, band, start, end, and orientation info.
    """
    chrom = coords.split(':')[0].replace('chr','')
    start,end = map(float,coords.split(':')[1].split('-'))
    start = int(start)
    end = int(end)
    band = cytoband_filtered[(cytoband_filtered['chr'] == chrom) & ((cytoband_filtered['start'] >= start) | (cytoband_filtered['end'] >= start)) & ((cytoband_filtered['end'] <= end)|(cytoband_filtered['start'] <= end))]['band'].values
    if len(band) > 1:
        band = '{}{}'.format(band[0],band[-1])
    if len(band) == 1:
        band = band[0]
    out_handle = "{chrom}{band}({start}-{end})x1({orientation})".format(chrom=chrom,band=band,start=start,end=end,orientation=orientation)
    return out_handle

def read_in_cyto(cytobands = 'resources/cytoBand.txt'):
    """
    Read and process cytoband data from a given file.

    Args:
        cytobands (str): Path to the cytoband file. Default is 'resources/cytoBand.txt'.

    Returns:
        pd.DataFrame: Filtered and processed cytoband data.
    """
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX', 'chrY', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']
    cytoband = pd.read_table(cytobands,header=None)
    cytoband.columns = ['chr','start','end','band','stain']
    cytoband_filtered = cytoband[cytoband['chr'].isin(chroms)]
    cytoband_filtered['chr'] = cytoband_filtered['chr'].str.replace("chr", '').replace('X','23').replace('Y','24')
    cytoband_filtered['band_base'] = cytoband_filtered['band'].str.split('.',expand=True)[0]
    cytoband_filtered['start'] = cytoband_filtered['start'].astype(int)
    cytoband_filtered['end'] = cytoband_filtered['end'].astype(int)
    return cytoband_filtered
def node_to_map(svs, xmap, g):
    """
    Create mapping dictionaries for nodes based on SVs.

    Args:
        svs (list): List of structural variations.
        xmap (xmap.parser): Information about the xmap.
        g (Graph): Graph structure.

    Returns:
        tuple: Two dictionaries - mapping of nodes to map and nodes to smap.
    """
    node_to_map_dict = {}
    node_to_smap_dict = {}
    for sv in svs:
        q_id = sv.q_id
        smap_id = sv.smap_id
        n_type1, n_type2 = detect_sv_directions(sv, xmap)
        #a and b are two nodes that ara connected by sv edge
        a = find_nodes(sv.ref_c_id1, sv.ref_start, g.vertices, n_type1)
        b = find_nodes(sv.ref_c_id2, sv.ref_end, g.vertices, n_type2)
        a = str(a)
        b = str(b)
        if a not in node_to_map_dict:
            node_to_map_dict[a] = []
            node_to_smap_dict[a] = []
        node_to_map_dict[a].append(q_id)
        node_to_smap_dict[a].append(smap_id)
        if b not in node_to_map_dict:
            node_to_map_dict[b] = []
            node_to_smap_dict[b] = []
        node_to_map_dict[b].append(q_id)
        node_to_smap_dict[b].append(smap_id)
    return node_to_map_dict, node_to_smap_dict

def return_mapids(v,u,node_to_map_dict):
    """
    Return map IDs for the given nodes v and u.

    Args:
        v (int): Node v.
        u (int): Node u.
        node_to_map_dict (dict): Mapping dictionary of nodes to map IDs.

    Returns:
        str: String representation of intersecting nodes and differences between nodes v and u.
    """
    v = str(v)
    u = str(u)
    if v in node_to_map_dict:
        start_node_maps = node_to_map_dict[v]
    else:
        start_node_maps = ''
    if u in node_to_map_dict:
        end_node_maps = node_to_map_dict[u]
    else:
        end_node_maps = ''
    intersect_nodes = set(start_node_maps).union(set(end_node_maps))
    start_dif = set(start_node_maps) - set(end_node_maps)
    end_dif = set(end_node_maps) - set(start_node_maps)
    map_str = '{}|{}|{}'.format(','.join(map(str,intersect_nodes)),','.join(map(str,start_dif)),','.join(map(str,end_dif)))
    return map_str

def convert_path(structure, path_map, cytoband_filtered):
    """
    Convert the provided structure into merged coordinates and ISCN coordinates.

    Args:
        structure (str): String representation of the structure.
        path_map (dict): Mapping of paths to coordinates.
        cytoband_filtered (pd.DataFrame): Filtered cytoband data.

    Returns:
        tuple: Merged coordinates and ISCN coordinates as strings.
    """
    structure_split = structure.split()
    coord_list = []
    iscn_list = []
    for s in structure_split:
        orientation = s[-1]
        p = s[:-1]
        coords = path_map[p]
        update_coords = '{}({})'.format(coords,orientation)
        coord_list.append(update_coords)
        iscn_mapped = map_cyto_coords(coords, cytoband_filtered, orientation)
        iscn_list.append(iscn_mapped)
    merged_coords = ','.join(coord_list)
    iscn_coords = ','.join(iscn_list)
    return merged_coords, iscn_coords

def smap_to_segment(pathnumber, smapids, smap):
    """ Associates smap ids with pathnumber, pulls smap entries

    Args:
        pathnumber (int): path
        smapids (_type_): _description_
        smap (list): list of parsed smap entries from parsers.SmapEntry
    """
    all_ids = smapids.split('|')[0].split(',')
    frame_list = []
    if (len(all_ids) >=1) & (all_ids[0] != ''):
        for smapid in all_ids:
            smap_entry = find_in_smap(int(smapid),smap)
            smap_frame = pd.DataFrame(columns=['smap_id','q_id','ref_c_id1','ref_c_id2','ref_start','ref_end','confidence','sv_type','size','VAF'],data=[[smapid, smap_entry.q_id, smap_entry.ref_c_id1, smap_entry.ref_c_id2, smap_entry.ref_start, smap_entry.ref_end, smap_entry.confidence, smap_entry.sv_type, smap_entry.size, smap_entry.VAF]],index=['Segment {}'.format(pathnumber)])
            frame_list.append(smap_frame)
    if len(frame_list)>0:
        joined_frame_list = pd.concat(frame_list)
    else:
        joined_frame_list = None
    return joined_frame_list

def find_sv_node_edges(svs, xmap, g):
    """
    Identify nodes in a graph that are connected by a structural variation (SV) edge.

    Parameters:
    - svs (list): A list of structural variations (SVs) where each SV has attributes:
        * q_id: Query identifier of the SV.
        * smap_id: Smap identifier of the SV.
        * ref_c_id1: Reference chromosome ID for the start position of the SV.
        * ref_start: Start position of the SV on the reference chromosome.
        * ref_c_id2: Reference chromosome ID for the end position of the SV.
        * ref_end: End position of the SV on the reference chromosome.
    - xmap (object): An object representing the cross-map for the SVs, used by the `detect_sv_directions` function.
    - g (object): A graph object with an attribute `vertices` representing all the vertices (or nodes) in the graph.

    Returns:
    - set: A set of tuple pairs, where each tuple (a, b) represents two nodes (a and b) in the graph
           that are connected by an SV edge. The set includes both (a, b) and (b, a) for bidirectionality.

    Note:
    The function makes use of two additional functions:
    1. `detect_sv_directions(sv, xmap)` - Determines the directions/types of the nodes based on the SV and xmap.
    2. `find_nodes(ref_c_id, position, vertices, n_type)` - Finds nodes in the graph based on chromosome ID,
       position, list of vertices, and node type.

    Example:
    Consider svs as a list where each element has attributes like q_id, smap_id, ref_c_id1, etc.
    Given svs, xmap, and g, the function will return a set of node pairs connected by SV edges.
    """
    tuples_list = []
    for sv in svs:
        q_id = sv.q_id
        smap_id = sv.smap_id
        n_type1, n_type2 = detect_sv_directions(sv, xmap)
        #a and b are two nodes that ara connected by sv edge
        a = find_nodes(sv.ref_c_id1, sv.ref_start, g.vertices, n_type1)
        b = find_nodes(sv.ref_c_id2, sv.ref_end, g.vertices, n_type2)
        tuples_list.append((a,b))
        tuples_list.append((b,a))
    sv_tuples_set = set(tuples_list)
    return sv_tuples_set

def associate_segments_to_svs(paths,g,sv_tuples_set,node_to_smap_dict,centro):
    """
    Associate segments to structural variations (SVs) by converting Eulerian paths with vertex IDs to segment paths.

    Parameters:
    - paths (list of lists): A list of Eulerian paths where each path is represented as a list of vertices.
    - g (object): A graph object representing the underlying graph.
    - sv_tuples_set (set): A set of tuple pairs, where each tuple (a, b) represents two nodes in the graph
                           connected by an SV edge.
    - node_to_smap_dict (dict): A dictionary mapping nodes to their corresponding SMAP IDs.

    Returns:
    - list: A list of segments associated to SVs, where each segment is represented as a list
            [Segment_Name, Path_Name, SMAP_IDs].

    Internal Functions and Steps:
    1. return_all_edges_in_cc: Returns all edges in a connected component of the graph.
    2. detect_segment_vertices: Identifies the vertices that belong to segments in the graph.
    3. check_non_centromeric_path: Checks if a path is non-centromeric.
    4. return_mapids: Returns the SMAP IDs for a pair of nodes.

    Example:
    Given paths, g, sv_tuples_set, and node_to_smap_dict, the function will return a list of segments associated to SVs.

    Notes:
    - This function prints intermediate paths.
    - Assumes some functions (like return_all_edges_in_cc, detect_segment_vertices, check_non_centromeric_path,
      return_mapids) exist in the surrounding context.
    """
    c = 1
    segs_list = []
    for p in paths:
        component = list(set(p))
        component_edges = return_all_edges_in_cc(component, g)
        segment_vertices = detect_segment_vertices(component, component_edges)
        ans = []
        temp = [p[0]]
        for i in range(1,len(p)-1):
            if p[i] in segment_vertices and p[i-1]==p[i+1]:
                temp.append(p[i])
                if check_non_centromeric_path(temp,g,centro):
                    ans.append(temp)
                temp = [p[i]]
            else:
                temp.append(p[i])
        temp.append(p[-1])
        if check_non_centromeric_path(temp,g,centro):
            ans.append(temp)
        ans2 = []
        for subpath in ans:
            temp = ''
            path = f"Path {c}"
            print(path)
            for i in range(0,len(subpath)-1,2):
                seg_number = int((max(subpath[i],subpath[i+1])+1)/2)
                try:
                    print((subpath[i+1],subpath[i+2]))
                    if (subpath[i],subpath[i+1]) in sv_tuples_set:
                        smap_ids = return_mapids(str(subpath[i]),str(subpath[i+1]),node_to_smap_dict)
                        segs_list.append([f"Segment {seg_number}", path, smap_ids])
                    if (subpath[i+1],subpath[i+2]) in sv_tuples_set:
                        smap_ids = return_mapids(str(subpath[i+1]),str(subpath[i+2]),node_to_smap_dict)
                        segs_list.append([f"Segment {seg_number}", path, smap_ids])
                except:
                    continue
            c += 1
    return segs_list
