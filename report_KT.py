import argparse
from collections import defaultdict
import pandas as pd
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
    for k in r:
        r[k] = [min(r[k]), max(r[k])]
    return r
def parse_decipher(dec_file):
    a = {}
    a = defaultdict(lambda: [], a)
    b = {}
    b = defaultdict(lambda: '', b)
    with open(dec_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0]!= 'gene symbol':
                a[line[0]].append(line[3]+':'+line[2])
                b[line[0]] = line[1]+':'+ line[0]
    return a, b


def parse_genes(genes_file):
    r = {}
    r = defaultdict(lambda: [], r)
    d = {}
    with open(genes_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            gene = line[3]
            # print(line[0],gene)
            chrom = line[0][3:]
            if len(chrom )< 4 :
                if chrom=='X':
                    chrom = '23'
                if chrom == 'Y':
                    chrom = '24'
                pos1 = int(line[1])
                pos2 = int(line[2])
                d[gene] = chrom
                r[gene] = [chrom,min(pos1, pos2), max(pos1,pos2),gene]
            # print(chrom,gene)
    # for k in r:
    #     r[k] = [d[k],min(r[k]), max(r[k]),k]
    return r

def find_duplicates(input_list):
    seen = set()
    duplicates = set()
    for item in input_list:
        if item in seen:
            duplicates.add(item)
        else:
            seen.add(item)
    return list(duplicates)
def intervals_overlap(interval1, interval2):
    start1, end1 = interval1
    start2, end2 = interval2

    # Check for non-overlapping conditions
    if end1 < start2 or end2 < start1:
        return 0
    elif(start2<=start1<=end2 and end2<=end1) or (start1<=start2<=end1 and end1<=end2):
        return 2
    else:
        return 1
def return_gene_list(genes,sv_segment,integration_point, segments):
    ans_genes = set()
    a = {}
    a = defaultdict(lambda: set(), a)
    for i in genes.values():
        for j in integration_point:
            # print(j, i)
            if j[0] == i[0]:
                if intervals_overlap((i[1],i[2]),(j[1],j[2]))==1:
                    ans_genes.add(i[-1])
                    a['chr'+j[0]+':'+str(j[1])+'-'+str(j[2])].add(i[-1])
                elif intervals_overlap((i[1],i[2]),(j[1],j[2]))==2:
                    ans_genes.add(i[-1]+'*')
                    a['chr'+j[0]+':'+str(j[1])+'-'+str(j[2])].add(i[-1]+'*')
        for j in sv_segment:
            if segments[j][0] == i[0]:
                if intervals_overlap((i[1],i[2]),(segments[j][1],segments[j][2]))==1:
                    ans_genes.add(i[-1])
                    a[j].add(i[-1])
                elif intervals_overlap((i[1],i[2]),(segments[j][1],segments[j][2]))==2:
                    ans_genes.add(i[-1]+'*')
                    a[j].add(i[-1]+'*')
    # print(ans_genes)
    return ans_genes ,a
def parse_file(file_path, centro,decipher):
    segments = {}
    with open(file_path, 'r') as file:
        for line in file:
            direction = []
            numbers = []
            if line.startswith('Path'):
                p = False
                path = line.strip().split('\t')[0].split('=')[1][1:]
                sv_segment = set()
                integration_point = []
                for i in path.split(' '):
                    direction.append(i[-1])
                    numbers.append(int(i[:-1]))
                if len(set(direction)) != 1:
                    p = True
                else:
                    if len(numbers) ==2 :
                        if abs(numbers[0]-numbers[1]) !=1:
                            p = True
                    else:
                        for i in range(1,len(numbers)-1):
                            if numbers[i-1] - numbers[i] != numbers[i] - numbers[i+1]:
                                p = True
                                break
                if p :
                    reverse_path = False
                    path_number = line.strip().split('\t')[0].split('=')[0][:-1]
                    chrom = ''
                    path = path.split(' ')
                    for i in range(len(path)):
                        n = path[i][:-1]
                        if segments[n][3]:
                            chrom = chrom + 'chr' + segments[n][0]+ ' '
                            if path[i][-1] == '-':
                                reverse_path = True
                            path[i] = path[i] + '*'

                    if reverse_path:
                        path.reverse()
                        for i in range(len(path)):
                            if '+'  in path[i]:
                                path[i] = path[i].replace('+','-')
                            else:
                                path[i] = path[i].replace('-','+')
                    for i in path:
                        if i[-1]=='*':
                            i = i[:-1]
                        if i[-1] == '-':
                            sv_segment.add(i[:-1])
                    a = find_duplicates(numbers)
                    for j in a:
                        sv_segment.add(str(j))
                    f = False
                    for i in range(len(numbers)-1):
                        if abs(numbers[i+1] - numbers[i]) >1:
                            a1 = segments[str(numbers[i])]
                            a2 = segments[str(numbers[i+1])]
                            if a1[0] == a2[0]:
                                a = []
                                a.append(a1[1])
                                a.append(a1[2])
                                a.append(a2[1])
                                a.append(a2[2])
                                a = sorted(a)
                                if a[2] - a[1] < 3000000:
                                    for j in range(min(numbers[i],numbers[i+1])+1,max(numbers[i],numbers[i+1])):
                                        sv_segment.add(str(j))
                            if a1[0] != a2[0] or (a1[0] == a2[0] and a[2] - a[1] > 3000000):
                                if direction[i] == '+':
                                    integration_point.append((a1[0],a1[2]-30000,a1[2]))
                                else:
                                    integration_point.append((a1[0],a1[1],a1[1]+30000))
                                if direction[i+1] == '+':
                                    integration_point.append((a2[0],a2[2]-30000,a2[2]))
                                else:
                                    integration_point.append((a2[0],a2[1],a2[1]+30000))


                    # print(integration_point)
                    # li = path_number + '='+ ' '.join(i for i in path) + '\tChromosome='+chrom + '\trev='+str(reverse_path)+'\tSV_segment='+','.join(i for i in sv_segment)+'\tIntegrationPoint='+','.join('chr'+str(i[0])+':'+str(int(i[1]))+'-'+str(int(i[2])) for i in integration_point)
                    name =args.input[-5:]
                    li = name[:1] + '\t'+path_number + '='+ ' '.join(i for i in path) + '\t'+chrom + '\t'+str(reverse_path)+'\t'+','.join(i for i in sv_segment)+'\t'+','.join('chr'+str(i[0])+':'+str(int(i[1]))+'-'+str(int(i[2])) for i in integration_point)
                    ans_gene, ans_gene2 = return_gene_list(genes,sv_segment,integration_point, segments)
                    dec_list = []
                    genes_dec_list  = []
                    for g in ans_gene:
                        g = g.replace('*','')
                        for dd in decipher[g]:
                            dec_list.append(dd)
                        if g in genes_mim.keys():
                            genes_dec_list.append(genes_mim[g])


                    # li = li + '\tgenes='+','.join(s for s in ans_gene)
                    # li = li + '\tdecipher='+','.join(s for s in dec_list)
                    # li = li + '\t'+';'.join(s for s in ans_gene)
                    li = li + '\t'+' '.join(s+':('+(';'.join(ss for ss in ans_gene2[s]))+')' for s in ans_gene2.keys())
                    li = li + '\t'+';'.join(s for s in genes_dec_list)
                    li = li + '\t'+'; '.join(s for s in dec_list)
                    print(li)
                    # print(line.strip().split('\t')[0])
            elif not line.startswith('Segment\tNumber'):
                line = line.strip().split('\t')
                segment_number = line[1]
                segment_chr = line[2]
                segment_start = float(line[3])
                segment_end = float(line[4])
                contains_centro = not(segment_end < centro['chr'+segment_chr][0] or centro['chr'+segment_chr][1] < segment_start)
                # print(segment_number, segment_start,segment_end, centro['chr'+segment_chr][0],centro['chr'+segment_chr][1],contains_centro)
                segments[segment_number] = (segment_chr,segment_start,segment_end, contains_centro)
    # print(segments)


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input_file", required=True)
parser.add_argument("-centro", "--centro", help="path to file contains centromere coordinates", required=False)
parser.add_argument("-genes", "--genes", help="path to file contains genes coordinates", required=False)
parser.add_argument("-dec", "--dec", help="path to decipher file", required=False)
args  = parser.parse_args()
centro = parse_centro(args.centro)
genes = parse_genes(args.genes)
decipher , genes_mim  = parse_decipher(args.dec)
# print(genes)
parse_file(args.input, centro,decipher)
