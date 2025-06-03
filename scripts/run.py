from scripts.parsers import *
import csv
import shutil
from KarReporter.KarUtils.read_OMKar_output import *
import sys
from scripts.visualization import Plot_graph
from scripts.eulerian_path_finder import *
from scripts.path_decomposition import *
def parse_inputs(cnv_path, smap_path, rcmap_path, xmap_path, centro_path):
    """
    Parses input files to extract genomic segments, structural variation (SV) data, centromere regions, and other required data.

    Args:
        cnv_path (str): Path to the CNV call file.
        smap_path (str): Path to the Smap file.
        rcmap_path (str): Path to the RCmap file.
        xmap_path (str): Path to the Xmap file.
        centro_path (str): Path to the centromere file.

    Returns:
        tuple: Parsed data including:
            - segments: List of genomic segments.
            - all_seg: List of all segments.
            - smap: Parsed structural variation data.
            - centro: Centromere information.
            - rcov: Coverage data.
            - rcop: Copy number data.
            - xmap: Alignment data.
    """
    segments, all_seg = parse_cnvcall(cnv_path)
    smap = parse_smap(smap_path)
    segments, all_seg = fix_coordinate(segments, all_seg, smap)
    segments.sort(key=lambda x: (int(x.chromosome), x.start, x.end))
    all_seg.sort(key=lambda x: (int(x.chromosome), x.start, x.end))
    centro = parse_centro(centro_path)
    segments = merge_segments_all_seg_smap(segments, all_seg, smap, centro)  # Need to debug this function
    segments.sort(key=lambda x: (int(x.chromosome), x.start))
    rcov, rcop = parse_rcmap(rcmap_path)
    xmap = parse_xmap(xmap_path)
    return segments, all_seg, smap, centro, rcov, rcop, xmap

def detect_and_filter_svs(smap, rcop, segments, xmap, centro):
    """
    Detects and filters structural variations (SVs) based on size, confidence, and copy number criteria.
    Also integrates SVs into genomic segments.

    Args:
        smap (list): List of Smap entries representing SVs.
        rcop (dict): Copy number data for chromosomes.
        segments (list): List of genomic segments.
        xmap (dict): Xmap alignment data.
        centro (dict): Centromere regions for chromosomes.

    Returns:
        tuple: Updated segments and filtered SVs.
    """
    chrY_cn = int(np.average(list(rcop['24'].values())) + 0.5)
    chrX_cn = round(np.average(list(rcop['23'].values())))
    if chrY_cn > 0:
        chrX_cn = max(1, chrX_cn)
    svs = []
    segments = extend_segments_cn(segments, all_seg)  # fill the gap between calls.
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
            new_seg.int_cn = cn_in_mask_N_region(k, new_seg.start, new_seg.end, rcop[k],
                                                 centro)  # assumption that sample is diploide. default CN = 2
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
            new_seg.int_cn = cn_in_mask_N_region(k, new_seg.start, new_seg.end, rcop[k],
                                                 centro)  # assumption that sample is diploide. default CN = 2
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
                    new_seg.int_cn = cn_in_mask_N_region(k, new_seg.start, new_seg.end, rcop[k],
                                                         centro)  # assumption that sample is diploide. default CN = 2
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

    with open(output2.replace('txt', 'bedpe'), 'w', newline='') as bed:
        bed_writer = csv.writer(bed, delimiter='\t')
        bed_writer.writerow(['chromosome1', 'bp1', 'chromosome2', 'bp2', 'SV_type', 'smap_id', 'confidence'])
        for i in smap:
            if not i.sv_type.startswith('inser'):
                bed_writer.writerow([i.ref_c_id1, i.ref_start,i.ref_c_id2, i.ref_end, i.sv_type, i.smap_id, i.confidence])
    for i in smap:
        # translocation applied filters.
        if i.sv_type.startswith('trans') and i.confidence >= 0.05 and not i.sv_type.endswith('common') and not i.sv_type.endswith('segdupe') and (
                i.ref_c_id1 != i.ref_c_id2 or abs(i.ref_end - i.ref_start) > 300000):
            svs.append(i)  # fixed
            exist, s = detect_receprical_translocation(i, xmap, smap)
            if exist:
                svs.append(s)
        # indels
        elif i.sv_type.startswith('inse') or i.sv_type.startswith('delet'):
            if not i.sv_type.endswith('nbase') and not i.sv_type.endswith('tiny') and i.confidence >= 0:
                if i.sv_type.startswith('delet') and i.size > 200000:
                    if detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end, segments)[0]:
                        _, i.ref_start, i.ref_end = detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end, segments)
                        svs.append(i)
                    elif i.confidence >= 0.98 and i.size > 2000000:  # This is a limit for detection thos misslabeld translocation as deletion in Bionano Pipeline
                        svs.append(i)
                elif i.size > 500000 and abs(i.ref_end - i.ref_start) > 500000:  # this would be for insertion length more than 500Kbp
                    svs.append(i)
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
        elif i.sv_type == 'inversion_paired' and i.confidence >= 0.7:  # if it is full inversion
            s = find_in_smap(i.linkID, smap)
            if abs(i.ref_start - s.ref_end) > 500000:
                i.ref_end = s.ref_end
                svs.append(i)
        elif i.sv_type == 'duplication_inverted':
            if detect_duplicatioon_inversion_cn(i, xmap, segments)[0]:
                _, fold_point = detect_duplicatioon_inversion_cn(i, xmap, segments)  # because CN is changed
                i.ref_start = fold_point
                i.ref_end = fold_point
                svs.append(i)
        elif i.sv_type == 'duplication' or i.sv_type == 'duplication_split':
            if detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end, segments)[0]:
                _, i.ref_start, i.ref_end = detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end, segments)
                svs.append(i)
        elif i.size > 500000 and not i.sv_type.startswith('inversion'):  # Other type of SV
            svs.append(i)

    for sv in svs:  # integrate BPs and Segments
        find_bp_in_segment(sv.ref_c_id1, sv.ref_start, segments)  #
        find_bp_in_segment(sv.ref_c_id2, sv.ref_end, segments)
    # merging Bps for spiliting segments
    for s in segments:
        s.bp = merge_list(s.bp)
    return segments, svs

def build_graph(segments, svs, output_path, name):
    # Graph Creation
    """
    Constructs a graph from genomic segments and structural variations.

    Args:
        segments (list): List of genomic segments.
        svs (list): List of structural variations.
        output_path (str): Path to save the graph and related files.
        name (str): Name of the dataset for labeling outputs.

    Returns:
        Graph: A graph object representing the genomic structure.
    """
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
                    g.edges.append(
                        (counter - 1, counter, 0, 'R'))  # addign reference edge between two continues segments
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
        if int(sv.ref_c_id1) > int(sv.ref_c_id2) or (
                int(sv.ref_c_id1) == int(sv.ref_c_id2) and sv.ref_end < sv.ref_start):
            a = find_nodes(sv.ref_c_id1, sv.ref_start, g.vertices, n_type2)
            b = find_nodes(sv.ref_c_id2, sv.ref_end, g.vertices, n_type1)
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
        elif sv.sv_type.startswith('delet'):  # telomere cite deletion prevent
            if (a, b, 0, 'SV') not in g.edges and abs(a - b) != 1:
                g.edges.append((a, b, 0, 'SV'))
                g.return_node(a).append_edges(b)
                g.return_node(b).append_edges(a)
        elif (a, b, 0, 'SV') not in g.edges:
            g.edges.append((a, b, 0, 'SV'))
            g.return_node(a).append_edges(b)
            g.return_node(b).append_edges(a)
    g.print_node()
    g.output_node(f"{output_path}/{name}/{name}.preILP_nodes.txt")
    print(g.edges)
    g.output_edges(f"{output_path}/{name}/{name}.preILP_edges.txt")
    file = f"{output_path}/{name}/{name}.pdf"
    Plot_graph(g, file, name, centro)
    return g
def process_connected_components(g, output_path, name):
    """
    Processes connected components in the graph to estimate edge multiplicities
    and identify Eulerian paths.

    Args:
        g (Graph): The graph object representing the genomic structure.
        output_path (str): Path to save the processed components.
        name (str): Name of the dataset for labeling outputs.

    Returns:
        tuple: Processed data including:
            - connected_components: List of connected components in the graph.
            - paths: List of Eulerian paths in the graph.
            - edges_with_dummy: List of edges with dummy edges added.
    """
    file2 = f"{output_path}/{name}/{name}_2.png"
    connected_components = find_connected_components(g)
    for component in connected_components:
        # if 102 in component:
        component_edges = estimating_edge_multiplicities_in_CC(component, g, xmap)
    connected_components = find_connected_components(g)
    Plot_graph(g, file2, name, centro)
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
        euler_tour, component_edges_with_dummies = printEulerTour(component, component_edges, g, out_file, centro)
        paths.append(euler_tour)
        edges_with_dummy.append(component_edges_with_dummies)
        component_counter += 1
    metadata_file = f"{output_path}/{name}/postILP_components/{name}.postILP.metadata.txt"
    with open(metadata_file, 'w') as fp_write:
        for key, value in component_metadata.items():
            fp_write.write("{}\t{}\n".format(key, value))
    return connected_components, paths, edges_with_dummy

def generate_outputs(output_path, name, g, connected_components, paths, edges_with_dummy):
    # write in the output
    """
    Generates and writes outputs, including raw and processed genomic segment data.

    Args:
        output_path (str): Path to save the output files.
        name (str): Name of the dataset for labeling outputs.
        g (Graph): The graph object representing the genomic structure.
        connected_components (list): List of connected components in the graph.
        paths (list): List of Eulerian paths in the graph.
        edges_with_dummy (list): List of edges with dummy edges added.

    Returns:
        None: Outputs are written to specified files.
    """
    segments_border_masked_regions_start = []
    segments_border_masked_regions_end = []
    output3 = f"{output_path}/{name}/{name}_flagged.txt"
    output4 = f"{output_path}/{name}/{name}.txt"
    output = f"{output_path}/{name}/{name}_raw.txt"
    with open(output3, 'w') as f2:
        with open(output, 'w') as f:
            f.write('Segment\tNumber\tChromosome\tStart\tEnd\tStartNode\tEndNode\n')
            f2.write('Segment\tNumber\tChromosome\tStart\tEnd\tStartNode\tEndNode\n')
            number = 1
            for i in range(0, len(g.vertices), 2):
                v = g.vertices[i]
                u = g.vertices[i + 1]
                f.write(
                    'Segment\t{id}\t{chrom}\t{start}\t{end}\t{snode}\t{enode}\n'.format(id=number, chrom=v.chromosome,
                                                                                        start=v.pos, end=u.pos,
                                                                                        snode=v.id,
                                                                                        enode=u.id))
                f2.write(
                    'Segment\t{id}\t{chrom}\t{start}\t{end}\t{snode}\t{enode}\n'.format(id=number, chrom=v.chromosome,
                                                                                        start=v.pos, end=u.pos,
                                                                                        snode=v.id,
                                                                                        enode=u.id))
                number += 1
            c = 1
            for path_idx, p in enumerate(paths):
                structures, scores = convert_path_to_segment(p, edges_with_dummy[path_idx], g, centro)
                for jj in range(len(structures)):
                    structure = structures[jj]
                    if structure.endswith(' '):
                        structure = structure[:-1]
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
                        if (direction == '-' and seg in segments_border_masked_regions_start) and i < len(
                                paths_structure) - 1:
                            next_dir = paths_structure2[i + 1][-1]
                            next_seg = int(paths_structure2[i + 1][:-1])
                            if next_dir != direction or next_seg - seg != -1:
                                paths_structure[i] += "*"
                        if (direction == '+' and seg in segments_border_masked_regions_start) and i > 0:
                            prev_dir = paths_structure2[i - 1][-1]
                            prev_seg = int(paths_structure2[i - 1][:-1])
                            if prev_dir != direction or prev_seg - seg != -1:
                                paths_structure[i] = "*" + paths_structure[i]
                        if (direction == '+' and seg in segments_border_masked_regions_end) and i < len(
                                paths_structure) - 1:
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
                                                   forbidden_region_file=get_metadata_file_path(
                                                       'acrocentric_telo_cen.bed'))
    if validation_code:
        path_list, segment_dict = read_OMKar_output(post_process_input, return_segment_dict=True)
        processed_path_list, segment_obj_to_idx_dict = post_process_OMKar_output(path_list,
                                                                                 gap_merge_allowance=5,
                                                                                 isolate_centromere=True,
                                                                                 forbidden_region_file=get_metadata_file_path(
                                                                                     'acrocentric_telo_cen.bed'))
        write_MK_file(post_process_output, processed_path_list, segment_obj_to_idx_dict)
    else:
        ## if validation failed, keep original omkar output
        print(
            f"OMKar output format validation failed, the interpretation may include un-intended calls; check log for detailed reason; file: {post_process_output}",
            file=sys.stderr)
        shutil.copy2(post_process_input, post_process_output)


def run_omkar(cnv_path, smap_path, rcmap_path, xmap_path, centro_path, name, output_path):
    """
    Main function to execute the OMKar pipeline. It parses inputs, detects SVs,
    builds a graph, processes connected components, and generates outputs.

    Args:
        cnv_path (str): Path to the CNV call file.
        smap_path (str): Path to the Smap file.
        rcmap_path (str): Path to the RCmap file.
        xmap_path (str): Path to the Xmap file.
        centro_path (str): Path to the centromere file.
        name (str): Name of the dataset for labeling outputs.
        output_path (str): Path to save all output files.

    Returns:
        None: Outputs are saved to the specified directory.
    """
    global centro, all_seg, output2, xmap
    segments, all_seg, smap, centro, rcov, rcop, xmap = parse_inputs(cnv_path, smap_path, rcmap_path, xmap_path, centro_path)
    os.makedirs(f'{output_path}/{name}/', exist_ok=True)
    output2 = f"{output_path}/{name}/{name}_SV.txt"
    segments, svs = detect_and_filter_svs(smap, rcop, segments, xmap, centro)
    g = build_graph(segments, svs, output_path, name)
    connected_components, paths, edges_with_dummy = process_connected_components(g, output_path, name)
    generate_outputs(output_path, name, g, connected_components, paths, edges_with_dummy)
