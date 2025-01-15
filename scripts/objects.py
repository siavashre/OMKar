class Vertices:  # Class of vertices in out graph
    def __init__(self):
        self.type = ''  # Type of vertices it can be Head or Tail
        self.id = 0  # Id of vertices
        self.chromosome = 0  # which chromosome this vertices belong chrX = 23 and chrY = 24
        self.pos = 0  # in terms of bp the position of this vertices
        self.cn = 0  # integer CN of segment this vertices represent
        self.edges = []  # its like adjacenty matrix it contains all adjacent vertices

    def append_edges(self, i):  # add new vertices to adjancty list
        self.edges.append(i)
        self.edges = list(set(self.edges))

    def remove_edges(self, i):  # remove vertices from adjacenty list
        if i in self.edges:
            self.edges.remove(i)
        else:
            print('Not here baba')

    def print_v(self):  # print vertices info
        print(self.id, self.chromosome, self.pos, self.cn, self.type, self.edges)


class Graph:  # Class of graph
    def __init__(self):
        self.vertices = []  # list of Object of vertices
        self.edges = []  # list of edges. each edges is tuple of (u,v,M,type) u and v are vertices connected to each other M represent edge multiplicities and type can be "S", "SV", "R", "D"
        self.terminal_vertices_ids = []  # subset of vertices that begin/end a chromosome

    def print_node(self):  # print all nodes
        for v in self.vertices:
            print(v.id, v.chromosome, v.pos, v.cn, v.type, v.edges)

    def output_node(self, out_file):
        out_str = ""
        for v in self.vertices:
            out_str += "{}\t{}\t{}\t{}\t{}\t{}\n".format(v.id, v.chromosome, v.pos, v.cn, v.type, v.edges)
        with open(out_file, 'w') as fp_write:
            fp_write.write(out_str)

    def output_edges(self, out_file):
        out_str = ""
        for e in self.edges:
            out_str += str(e) + '\n'
        with open(out_file, 'w') as fp_write:
            fp_write.write(out_str)

    def return_node(self, id):  # return node by id
        for v in self.vertices:
            if id == v.id:
                return v
        return None

    def return_edges(self, a, b):  # return all type of edge between two nodes. it return a list of edges or None
        ans = []
        for e in self.edges:
            if (a == e[0] and b == e[1]) or (a == e[1] and b == e[0]):
                ans.append(e)
        if len(ans) > 0:
            return ans
        return None

    def add_dummy_edges(self, u, v, cnn):  # add dummy edges with type "D" between two nodes if already edge exist between them increase the CN
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
        e_list = self.return_edges(a, b)
        for e in e_list:
            if e[3] == e_type:
                self.edges.remove(e)
                if count > 0:
                    self.edges.append((e[0], e[1], count, e[3]))

    def label_terminal_vertices(self):
        ## mark the first and the last vertices on this graph, for each chromosome
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

