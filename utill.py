import itertools

######################### DOUBLE CHECK THIS FUNCTION ######################################
def detect_sv_directions(sv, xmap):  # this function detect the direction of one smap like H/T to H/T. #need to be check Contain bug?
    if sv.sv_type == 'inversion_paired':  # for inversion_paired always return T to T
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
    # dir1_a = dir1
    # dir2_a = dir2
    # dir1_b = dir1
    # dir2_b = dir2
    # if swap == 1:
    #     dir1_a, dir2_a = dir2, dir1
    # if int(alignment1['RefContigID']) > int(alignment2['RefContigID']) or (int(alignment1['RefContigID']) == int(alignment2['RefContigID'])
    #      and alignment1['RefStartPos']> alignment2['RefStartPos']):
    #     dir1_b, dir2_b = dir2, dir1
    # if dir1_a !=dir1_b:
    #     print('Here we had BUG BUG BUG',sv)
    # return dir1_b, dir2_b


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

def check_non_centromeric_path(p, g, centro):
    count = 0
    for i in range(0, len(p) - 1, 2):
        u = g.return_node(p[i])
        v = g.return_node(p[i + 1])
        if u.chromosome == v.chromosome:
            # if (min(u.pos,v.pos)< min(centro['chr'+str(u.chromosome)]) and max(u.pos, v.pos)> min(centro['chr'+str(u.chromosome)])) or (min(u.pos,v.pos)< max(centro['chr'+str(u.chromosome)]) and max(u.pos, v.pos)> max(centro['chr'+str(u.chromosome)])):
            if is_overlapping(min(u.pos, v.pos), max(u.pos, v.pos), min(centro['chr' + str(u.chromosome)]), max(centro['chr' + str(u.chromosome)])):
                count += 1
    return count


def detect_segment_vertices(component,
                            edges):  # this function return those nodes which are telomere part (start and end of chromosome or no any edge like R or SV edge connect to them)
    v = set()
    for e in edges:
        if e[3] != 'S':
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
