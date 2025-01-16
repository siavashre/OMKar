from scripts.utill import *

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

def convert_path_to_segment(p, component_edges,g, centro):  # this is important function that convert Eulerion path with vertices ID to segment path.
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

    return  adjust_and_remove_overlapping_segments(segments), adjust_and_remove_overlapping_segments(all_seg)
