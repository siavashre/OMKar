import pandas as pd


def get_genes_in_region(chrom, start, end, gene_file='KT_Reports/Metadata/gtf_protein_coding.bed'):
    """
    report all genes with intersection with the region (inclusive)
    :param chrom: chr{1-22, X, Y}
    :param start:
    :param end:
    :param gene_file: bed file containing genes and location
    :return:
    """
    gene_df = pd.read_csv(gene_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'gene'])
    chrom = chrom[:3].lower() + chrom[3:]

    intersection_gene_df = gene_df[gene_df['chrom'] == chrom]
    intersection_gene_df = intersection_gene_df[intersection_gene_df['start'] <= end]
    intersection_gene_df = intersection_gene_df[intersection_gene_df['end'] >= start]
    intersection_gene_df = intersection_gene_df.sort_values(by='start')
    intersection_gene_df = intersection_gene_df.drop_duplicates(subset=['gene'], keep='first')
    return intersection_gene_df['gene'].tolist()


def get_DDG_overlapped_genes(input_gene_list, DDG_file='KT_Reports/Metadata/DDG2P_14_11_2023.csv'):
    DDG_df = pd.read_csv(DDG_file, sep='\t')
    overlapped_gene_df = DDG_df[DDG_df['gene symbol'].isin(input_gene_list)]
    return overlapped_gene_df


def tostring_gene_disease_omim(filtered_DDG_df):
    # for the same gene (under same gene OMIM), we can have multiple diseases (different disease OMIM)
    gene_list = []  # [(str, int)]; (gene, OMIM)
    disease_list = []  # [[(str, int)]]; [(disease name, corresponding OMIM)]

    genes = filtered_DDG_df['gene symbol'].unique()
    for gene in genes:
        gene_df = filtered_DDG_df[filtered_DDG_df['gene symbol'] == gene]
        diseases = []
        first_row = True
        for index, row in gene_df.iterrows():
            if first_row:
                gene_list.append((gene, row['gene mim']))
                first_row = False
            disease_name = row['disease name'].replace('&', '\&')  # prevent latex errors
            disease_mim = row['disease mim']
            diseases.append((disease_name, disease_mim))
        disease_list.append(diseases)

    return gene_list, disease_list


def test_if_DDG_has_duplicate():
    DDG_file = 'Metadata/DDG2P_14_11_2023.csv'
    DDG_df = pd.read_csv(DDG_file, sep='\t')
    print(DDG_df.shape)
    filtered_DDG_df = DDG_df.drop_duplicates(subset=['gene symbol'], keep=False)
    print(filtered_DDG_df.shape)
    overlapping_rows = DDG_df.groupby('gene symbol').filter(lambda x: len(x) > 1)
    print(overlapping_rows)


def get_band_location(chrom, nt_idx, cyto_file='KT_Reports/Metadata/cytoBand.txt'):
    cyto_df = pd.read_csv(cyto_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'band_name', 'stain'])
    chrom = chrom[:3].lower() + chrom[3:]
    for index, row in cyto_df.iterrows():
        if row['chrom'] != chrom:
            continue
        if row['start'] <= nt_idx <= row['end']:
            return row['band_name']
    raise RuntimeError('no band found')


def test_get_genes():
    genes = get_genes_in_region('Chr22', 25200725, 25560371)
    print(genes)
    df = get_DDG_overlapped_genes(genes)
    a, b = tostring_gene_disease_omim(df)
    for idx, a_itr in enumerate(a):
        print(a_itr, b[idx])


def test_get_band():
    x = get_band_location('Chr22', 25200725)
    print(x)


if __name__ == "__main__":
    test_get_band()
