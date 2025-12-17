class Vertices:  # Class of vertices in out graph
    """
    Represents a vertex in a graph, corresponding to genomic information.

    Attributes:
        type (str): Type of the vertex, e.g., "Head" or "Tail".
        id (int): Unique identifier for the vertex.
        chromosome (int): Chromosome number the vertex belongs to (e.g., 23 for chrX, 24 for chrY).
        pos (int): Genomic position of the vertex in base pairs.
        cn (int): Integer copy number for the segment represented by this vertex.
        edges (list): List of adjacent vertices (adjacency matrix).
    """
    def __init__(self):
        self.type = ''  # Type of vertices it can be Head or Tail
        self.id = 0  # Id of vertices
        self.chromosome = 0  # which chromosome this vertices belong chrX = 23 and chrY = 24
        self.pos = 0  # in terms of bp the position of this vertices
        self.cn = 0  # integer CN of segment this vertices represent
        self.edges = []  # its like adjacenty matrix it contains all adjacent vertices

    def append_edges(self, i):  # add new vertices to adjancty list
        """
        Adds a new vertex to the adjacency list of the current vertex.

        Args:
            i (int): Vertex ID to add to the adjacency list.
        """
        self.edges.append(i)
        self.edges = list(set(self.edges))

    def remove_edges(self, i):  # remove vertices from adjacenty list
        """
        Removes a vertex from the adjacency list of the current vertex.

        Args:
            i (int): Vertex ID to remove from the adjacency list.
        """
        if i in self.edges:
            self.edges.remove(i)
        else:
            print('No vertices here')

    def print_v(self):  # print vertices info
        """
        Prints the vertex information, including ID, chromosome, position, copy number, type, and adjacent vertices.
        """
        print(self.id, self.chromosome, self.pos, self.cn, self.type, self.edges)


class Graph:  # Class of graph
    """
    Represents a graph structure for genomic information.

    Attributes:
        vertices (list): List of `Vertices` objects in the graph.
        edges (list): List of edges, each represented as a tuple (u, v, M, type), where:
                      u and v are vertex IDs, M is edge multiplicity, and type indicates the edge type.
        terminal_vertices_ids (list): List of vertex IDs that mark the start or end of chromosomes.
    """
    def __init__(self):
        self.vertices = []  # list of Object of vertices
        self.edges = []  # list of edges. each edges is tuple of (u,v,M,type) u and v are vertices connected to each other M represent edge multiplicities and type can be "S", "SV", "R", "D"
        self.terminal_vertices_ids = []  # subset of vertices that begin/end a chromosome

    def print_node(self):  # print all nodes
        """
        Prints information for all nodes in the graph.
        """
        for v in self.vertices:
            print(v.id, v.chromosome, v.pos, v.cn, v.type, v.edges)

    def output_node(self, out_file):
        """
        Outputs node information to a file.

        Args:
            out_file (str): Path to the output file.
        """
        out_str = ""
        for v in self.vertices:
            out_str += "{}\t{}\t{}\t{}\t{}\t{}\n".format(v.id, v.chromosome, v.pos, v.cn, v.type, v.edges)
        with open(out_file, 'w') as fp_write:
            fp_write.write(out_str)

    def output_edges(self, out_file):
        """
        Outputs edge information to a file.

        Args:
            out_file (str): Path to the output file.
        """
        out_str = ""
        for e in self.edges:
            out_str += str(e) + '\n'
        with open(out_file, 'w') as fp_write:
            fp_write.write(out_str)

    def return_node(self, id):  # return node by id
        """
        Returns a node object by its ID.

        Args:
            id (int): Vertex ID.

        Returns:
            Vertices or None: The vertex with the given ID, or None if not found.
        """
        for v in self.vertices:
            if id == v.id:
                return v
        return None

    def return_edges(self, a, b):  # return all type of edge between two nodes. it return a list of edges or None
        """
        Returns all edges between two nodes.

        Args:
            a (int): First node ID.
            b (int): Second node ID.

        Returns:
            list or None: List of edges between the nodes, or None if no edges exist.
        """
        ans = []
        for e in self.edges:
            if (a == e[0] and b == e[1]) or (a == e[1] and b == e[0]):
                ans.append(e)
        if len(ans) > 0:
            return ans
        else:
            return None

    def add_dummy_edges(self, u, v, cnn):  # add dummy edges with type "D" between two nodes if already edge exist between them increase the CN
        """
        Adds dummy edges between two nodes, or increments copy number if an edge already exists.

        Args:
            u (int): ID of the first node.
            v (int): ID of the second node.
            cnn (int): Copy number to add for the edge.
        """
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
        """
        Updates an edge's multiplicity between two nodes for a specific edge type.

        Args:
            a (int): First node ID.
            b (int): Second node ID.
            count (int): New multiplicity for the edge.
            e_type (str): Type of the edge to update.
        """
        e_list = self.return_edges(a, b)
        for e in e_list:
            if e[3] == e_type:
                self.edges.remove(e)
                if count > 0:
                    self.edges.append((e[0], e[1], count, e[3]))

    def label_terminal_vertices(self):
        ## mark the first and the last vertices on this graph, for each chromosome
        """
        Identifies and labels the terminal vertices (start and end points) of chromosomes in the graph.
        """
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

class Segments: #Class represent each genomic segment
    """
    Represents a genomic segment.

    Attributes:
        id (str): Unique identifier for the segment.
        chromosome (int): Chromosome number.
        start (int): Start position of the segment.
        end (int): End position of the segment.
        width (int): Width of the segment.
        type (str): Segment type.
        fractional_cn (float): Fractional copy number.
        int_cn (int): Integer copy number.
        conf (float): Confidence score.
        line (str): Original line from the input file.
        bp (list): List of breakpoints within the segment.
    """
    id = ''
    chromosome = 0
    start = 0
    end = 0
    width = 0
    type = ''
    fractional_cn = 0
    int_cn = 0
    conf = 0
    line = ''
    bp = [] # this list contains all SVs intersect middle of this genomic segment including start and end position of segment, we will use this then for spiliting this to new genomic segments
class BP: #Class of SV breakpoints
    """
    Represents a structural variation (SV) breakpoint.

    Attributes:
        contig_id (str): Contig ID where the breakpoint is located.
        direction1 (str): Direction of the first breakpoint ("+" or "-").
        direction2 (str): Direction of the second breakpoint ("+" or "-").
        pos1 (str): Lower genomic coordinate of the breakpoint.
        pos2 (str): Upper genomic coordinate of the breakpoint.
        chrom1 (str): Chromosome of the first breakpoint.
        chrom2 (str): Chromosome of the second breakpoint.
        line (str): Original line from the input file.
        type (str): Type of the breakpoint.
    """
    contig_id = '' # contig id the breakpoint is called
    direction1 = '' # Forward direction is + / Reverse direction is -
    direction2 = '' # Forward direction is + / Reverse direction is -
    pos1 = '' # Always lower genomic coordinate and lower chromosomes are sorted
    pos2 = ''
    chrom1 = ''
    chrom2 = ''
    line = ''
    type = ''

class SmapEntry: #Class of each entry in Smap file
    """
    Represents an entry in an Smap file.

    Attributes:
        smap_id (str): Smap entry ID.
        q_id (str): Query contig ID.
        ref_c_id1 (str): Reference contig ID 1.
        ref_c_id2 (str): Reference contig ID 2.
        ref_start (float): Start position on the reference.
        ref_end (float): End position on the reference.
        query_start (float): Start position on the query.
        query_end (float): End position on the query.
        confidence (float): Confidence score.
        xmap_id1 (str): Xmap ID 1.
        xmap_id2 (str): Xmap ID 2.
        sv_type (str): Type of structural variation.
        line (str): Original line from the input file.
        size (float): Size of the structural variation.
        linkID (str): Link ID.
        VAF (float): Variant allele frequency.
    """
    smap_id = ''
    q_id = ''
    ref_c_id1 = ''
    ref_c_id2 = ''
    ref_start = 0
    ref_end = 0
    query_start = 0
    query_end = 0
    confidence = 0
    xmap_id1 = ''
    xmap_id2 = ''
    sv_type = ''
    line = ''
    size = 0
    linkID = ''
    VAF = 0
