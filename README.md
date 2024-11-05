
# OMKar

OMKar is a computational tool for automated karyotyping that utilizes Structural Variation (SV) and Copy Number Variation (CNV) calls derived from Optical Genome Mapping (OGM) data. OMKar generates a virtual karyotype by analyzing SV and CNV information, providing insights into chromosomal abnormalities associated with constitutional genetic disorders.

## Table of Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [File Structure](#file-structure)
- [Output and Interpretation](#output-and-interpretation)
- [Known Issues and Limitations](#known-issues-and-limitations)
- [License](#license)
- [Contact](#contact)

---
### Requirements
- **Python packages**: All dependencies are Python-based for ease of installation.
  - Python 3.6 or later.
  - `numpy` >=1.25.0
  - `matplotlib` >=3.8.0
  - `pandas` >=2.1.0
  - `scipy` >=1.11.0
  - `pulp` >=2.7.0
  - `argparse` >=1.4.0
- **Optional PDF Reports**: Requires TeXWork for generating PDF reports (HTML report generation only requires Python packages).
  `
### Installation
Clone OMKar from GitHub and navigate to the working directory:
```shell
git clone https://github.com/siavashre/OMKar.git
```
### File Structure
The input data directory must follow this structure:
```
<data_dir>/
    <sample_name1>/
        cnv_calls_exp.txt
        cnv_rcmap_exp.txt
        exp_refineFinal1_merged.xmap
        exp_refineFinal1_merged_filter_inversions.smap
    <sample_name2>/
        cnv_calls_exp.txt
        cnv_rcmap_exp.txt
        exp_refineFinal1_merged.xmap
        exp_refineFinal1_merged_filter_inversions.smap
    ...
```

### Usage
Run OMKar in the `OMKar/` directory using the batch script:
```shell
bash batch_run.sh data_dir [--out_dir out_dir] [--report report] [--debug debug]
```

**Arguments**:
| Argument                         | Type | Description                                                                                        |
|----------------------------------|------|----------------------------------------------------------------------------------------------------|
| `data_dir`                       | DIR  | Contains the data in the correct hierarchy                                                         |
| (optional) `--out_dir <outdir>`  | DIR  | (default: OMKar/outputs/) Output location                                                          |
| (optional) `--report <html/pdf>` | STR  | (default: None) Generating an HTML/PDF report for all cases ran                                    |
| (optional) `--debug true`        | STR  | (default: false) Only useful when generating a report; include debug information in the report output |

---
**Usage**:
For running OMkar(Only the karyotype part) on a single case you can use:
```shell
python3 main.py -cnv cnv_calls_exp.txt -smap exp_refineFinal1_merged_filter_inversions.smap -rcmap cnv_rcmap_exp.txt -xmap exp_refineFinal1_merged.xmap -n test -o output/ -centro hg38_centro.txt -cyto cytoBand.txt
```
- `-smap`: Path to the structural variations map file (`exp_refineFinal1_merged_filter_inversions.smap`).
- `-rcmap`: Path to the reference copy number map (`cnv_rcmap_exp.txt`).
- `-xmap`: Path to the XMAP alignment file (`exp_refineFinal1_merged.xmap`).
- `-n`: Name of the output file (`test`).
- `-o`: Output directory path for storing OMKar results.
- `-centro`: Path to the centromere file (`hg38_centro.txt`).
- `-cyto`: Path to the cytoband file (`cytoBand.txt`).

As an example, you can use the test files in the `test_files` directory. If you have installed the prerequisites correctly, you can run the following command to see the outputs of OMKar for the test case which contains Balanced translocation between
Chr2 and Chr14, and duplication in Chr2:

```shell
python3 main.py \
  -cnv test_files/cnv_calls_exp.txt \
  -smap test_files/exp_refineFinal1_merged_filter_inversions.smap \
  -rcmap test_files/cnv_rcmap_exp.txt \
  -xmap test_files/exp_refineFinal1_merged.xmap \
  -n test \
  -o test_files/ \
  -centro hg38_centro.txt \
  -cyto  cytoBand.txt
```
### Output
After running OMKar, the following output files are generated, providing various analyses and visualizations of the karyotype and structural variations:

1. **Chromosomal Graph PDF (`{name}.pdf`)**:
   - A PDF file with a graphical representation of each chromosome. This document shows detected structural variations and other chromosomal features across all chromosomes for easy visualization and interpretation.
2. **Segment Pathway File (`{name}.txt`)**:
   - This file provides a detailed representation of the karyotype in terms of chromosomal segments. OMKar outputs a "Molecular Karyotype" in this text file format, which includes:
      - **Segment List**: Each defined segment across the reference genome is listed with its segment number, chromosome, start and end coordinates, and the graph nodes representing the segment. Each segment has two nodes, connected by a segment edge, and is forward-oriented (end coordinate ≥ start coordinate). Segments are sorted by chromosome and coordinate, are non-overlapping, and exclude telomeric regions where applicable.
      - **Reconstructed Paths**: These paths represent the karyotype by listing segments in the format "Path number = segment number followed by direction." The direction (`+` or `-`) indicates traversal direction, with `+` being forward (start to end) and `-` being reverse (end to start). For instance, "Path1 = 1+ 2+ 3-" means that the path traverses segment 1 forward, segment 2 forward, and segment 3 in reverse.
      - **Centromere Count**: Each path includes the number of centromeres present, which should ideally be one to indicate a valid chromosome structure.
3. **Structural Variation BED File (`{name}_SV.bed`)**:
   - A BED format file detailing the locations of structural variations (SVs) identified in the genome. This file is useful for visualization in genome browsers and compatibility with other bioinformatics tools (e.g., UCSC Genome Browser).

4. **Structural Variation Summary File (`{name}_SV.txt`)**:
   - A summary of structural variations identified, including types (e.g., deletion, inversion, duplication), positions, confidence scores, and additional details like zygosity and allele frequency in the smap file format. 


### License
MIT License.

### Contact
For issues, questions, or contributions, visit [GitHub Issues](https://github.com/siavashre/OMKar/issues) or contact the maintainers directly.
