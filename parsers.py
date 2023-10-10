from collections import defaultdict
from os.path import exists

#################################### parsing RefAligner Copy Number file ##########################################
def parse_rcmap(cmap_dir):
    cov = {}
    cov = defaultdict(lambda: {}, cov)
    cop = {}
    cop = defaultdict(lambda: {}, cop)
    with open(cmap_dir, 'r') as f:
        for line in f:
            if line.startswith("#"):
                head = line[1:].rstrip().rsplit()
            if not line.startswith('#'):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                chrom = fD['CMapId']
                pos = int(float(fD['Position']))
                cover = int(fD['Coverage'])
                copynumber = int(fD['CopyNumber'])
                cov[chrom][pos] = cover
                cop[chrom][pos] = copynumber
    return cov, cop # return two dictionary of dictionary which keys are chromosome number and keys are label position and return coverage/ copynumber
    #cov return coverage / #cop return copynumber

class BP: #Class of SV breakpoints 
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

#################################### parsing RefAligner smap file ##########################################
def parse_smap(smap_dir):
    with open(smap_dir, 'r') as f:
        bfb_count = {}
        bfb_count = defaultdict(lambda: [], bfb_count)
        breakpoints = []
        translocations = []
        for line in f:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]
            if not line.startswith('#'):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                smap_entry = SmapEntry()
                smap_entry.ref_start = float(fD['RefStartPos'])
                smap_entry.ref_end = float(fD['RefEndPos'])
                smap_entry.xmap_id1 = int(fD['XmapID1'])
                smap_entry.xmap_id2 = int(fD['XmapID2'])
                smap_entry.q_id = int(fD['QryContigID'])
                smap_entry.ref_c_id1 = fD['RefcontigID1']
                smap_entry.ref_c_id2 = fD['RefcontigID2']
                smap_entry.smap_id = int(fD['SmapEntryID'])
                smap_entry.confidence = float(fD['Confidence'])
                smap_entry.query_start = float(fD['QryStartPos'])
                smap_entry.query_end = float(fD['QryEndPos'])
                smap_entry.sv_type = fD['Type']
                smap_entry.size =float(fD['SVsize'])
                smap_entry.line = line
                smap_entry.linkID = int(fD['LinkID'])
                # smap_entry.VAF = float(fD['VAF'])
                breakpoints.append(smap_entry)
    return breakpoints # return list of smap entries

class Segments: #Class represent each genomic segment
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

######################## parse CNV call ########################
def parse_cnvcall(cnvcall):
    segment_list = []
    all_seg = []
    with open(cnvcall, 'r') as f:
        for line in f:
            if line.startswith("#Id"):
                head = line.rstrip().rsplit()
            if not line.startswith('#'):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                segment = Segments()
                segment.id = int(fD['#Id'])
                segment.chromosome = fD['Chromosome']
                segment.start = float(fD['Start'])
                segment.end = float(fD['End'])
                segment.width = float(fD['Width'])
                segment.type = fD['Type']
                segment.fractional_cn = float(fD['fractionalCopyNumber'])
                segment.int_cn = int(fD['CopyNumber'])
                segment.conf = float(fD['Confidence'])
                segment.line = line
                segment.bp = [segment.start, segment.end]
                if segment.width > 200000  and not segment.type.endswith('masked') and segment.conf>= 0.95: #Apply filters on CNV call masked region and segments legnth 200000bp is a limit of filtering segments
                    if segment.width < 500000:
                        if len(segment_list) == 0:
                            segment_list.append(segment)
                        elif segment.int_cn >= segment_list[-1].int_cn: #If segment length is between 200 and 500 Kbp. if the CN is great than previouse segment length will add it( this prevent having small deletions between 200 and 500 Kbp)
                            segment_list.append(segment)
                    else:# if segment length is greater than 500Kbp add it. 
                        segment_list.append(segment)
                all_seg.append(segment)
    return segment_list, all_seg #return two lists of segment class. one the filtered one and one all of them. 

#################################### parsing Alignment xmap file ##########################################
def parse_xmap(xmapf):
    detailFields = ["XmapEntryID", "QryContigID", "RefContigID", "Orientation", "QryLen", "RefLen",
                    "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", "Alignment","Confidence"]
    numeric = ["QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos","Confidence"]
    xmapPair = {}
    with open(xmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                alnstring = ")" + fD["Alignment"] + "("
                # xmapAln[fD["XmapEntryID"]] = alnstring.rsplit(")(")[1:-1]

                xmapPair[fD["XmapEntryID"]] = {x: fD[x] for x in detailFields}
                for keyword in numeric:
                    xmapPair[fD["XmapEntryID"]][keyword] = float(xmapPair[fD["XmapEntryID"]][keyword])

    return xmapPair #return  dict of dict which key is Xmapentry 

################################ Parse bed file ###########################
def parse_bed(bed_dir):
    l = []
    with open(bed_dir, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            chrom = int(line[0][3:])
            start = int(float(line[1]))
            end = int(float(line[2]))
            l.append([chrom, start, end])
    return l 
############################################# 
#This function parse the centromere region
def parse_centro(centro):
    r = {}
    r = defaultdict(lambda: [], r)
    with open(centro, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            key = line[0]
            pos1 = int(line[1])
            pos2 = int(line[2])
            r[key].append(pos1)
            r[key].append(pos2)
    return r

