from collections import defaultdict
from os.path import exists
from scripts.objects import Segments, SmapEntry, BP

#################################### parsing RefAligner Copy Number file ##########################################
def parse_rcmap(cmap_dir):
    """
    Parses a RefAligner Copy Number file to extract coverage and copy number information.

    Args:
        cmap_dir (str): Path to the RefAligner Copy Number file.

    Returns:
        tuple: Two nested dictionaries:
            - `cov`: Coverage data where keys are chromosome IDs and values are dictionaries of positions to coverage.
            - `cop`: Copy number data where keys are chromosome IDs and values are dictionaries of positions to copy numbers.
    """
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
                cover = float(fD['Coverage'])
                copynumber = int(fD['CopyNumber'])
                cov[chrom][pos] = cover
                cop[chrom][pos] = copynumber
    return cov, cop # return two dictionary of dictionary which keys are chromosome number and keys are label position and return coverage/ copynumber
    #cov return coverage / #cop return copynumber

#################################### parsing RefAligner smap file ##########################################
def parse_smap(smap_dir):
    """
    Parses a RefAligner Smap file to extract structural variation breakpoints.

    Args:
        smap_dir (str): Path to the Smap file.

    Returns:
        list: A list of `SmapEntry` objects representing structural variation breakpoints.
    """
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

######################## parse CNV call ########################
def parse_cnvcall(cnvcall):
    """
    Parses a CNV (Copy Number Variation) call file to extract genomic segments.

    Args:
        cnvcall (str): Path to the CNV call file.

    Returns:
        tuple: Two lists of `Segments` objects:
            - Filtered segments based on length, confidence, and type.
            - All segments, including those filtered out.
    """
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
                if 'Width' in fD.keys():
                    segment.width = float(fD['Width'])
                else:
                    segment.width = float(fD['Size'])
                segment.type = fD['Type']
                segment.fractional_cn = float(fD['fractionalCopyNumber'])
                segment.int_cn = int(fD['CopyNumber'])
                segment.conf = float(fD['Confidence'])
                segment.line = line
                segment.bp = [segment.start, segment.end]
                if segment.width > 200000  and not segment.type.endswith('masked') and segment.conf>= 0.98: #Apply filters on CNV call masked region and segments legnth 200000bp is a limit of filtering segments
                        segment_list.append(segment)
                all_seg.append(segment)
    return segment_list, all_seg #return two lists of segment class. one the filtered one and one all of them. 

#################################### parsing Alignment xmap file ##########################################
def parse_xmap(xmapf):
    """
    Parses an Xmap file to extract alignment information.

    Args:
        xmapf (str): Path to the Xmap file.

    Returns:
        dict: A dictionary where keys are Xmap entry IDs and values are dictionaries of alignment details.
    """
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
    """
    Parses a BED file to extract genomic regions.

    Args:
        bed_dir (str): Path to the BED file.

    Returns:
        list: A list of lists, each containing chromosome ID, start position, and end position of a region.
    """
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
    """
    Parses a centromere file to extract centromere regions.

    Args:
        centro (str): Path to the centromere file.

    Returns:
        dict: A dictionary where keys are chromosome IDs and values are lists of start and end positions of centromeres.
    """
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


def parse_forbiden_region(forbiden):
    """
    Parses a forbidden regions file to extract masked genomic regions.

    Args:
        forbiden (str): Path to the forbidden regions file.

    Returns:
        dict: A dictionary where keys are chromosome IDs and values are lists of masked regions.
              Each region is represented as a list containing start position, end position, and region type.
    """
    masked_regions = {}
    with open(forbiden, 'r') as file:
        lines = file.readlines()
    for line in lines[1:]:
        columns = line.strip().split('\t')
        if columns[0] == 'ChrX':
            columns[0] = 23
        elif columns[0] == 'ChrY':
            columns[0] = 24
        else:
            columns[0] = int(columns[0][3:])
        region = {
            'Chr': columns[0],
            'StartPos': int(columns[1]),
            'EndPos': int(columns[2]),
            'Type': columns[3]
        }
        if columns[0] not in masked_regions:
            masked_regions[columns[0]] = []
        masked_regions[columns[0]].append([int(columns[1]), int(columns[2]), columns[3]])
    return masked_regions
