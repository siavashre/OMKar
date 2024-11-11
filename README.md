
# OMKar

OMKar is a computational tool for automated karyotyping that utilizes Structural Variation (SV) and Copy Number Variation (CNV) calls derived from Optical Genome Mapping (OGM) data. OMKar generates a virtual karyotype by analyzing SV and CNV information, providing insights into chromosomal abnormalities associated with constitutional genetic disorders.

## Table of Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Input File Structure](#input-file-structure)
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
  `
### Installation
Clone OMKar from GitHub and navigate to the working directory:
```shell
git clone --recurse-submodules https://github.com/siavashre/OMKar.git
```
*The submodules are needed for generating the HTML reports.*

### Input File Structure
The input data directory needs to be either the Bionano Solve's default output structure, or a curated directory structure:

Bionano default output DIR structure:
```
<data_dir>/
    <sample_name1>/
        contigs/
            alignmolvref/
                copynumber/
                    cnv_calls_exp.txt
                    cnv_rcmap_exp.txt
            exp_refineFinal1_sv/
                merged_smaps/
                    exp_refineFinal1_merged.xmap
                    exp_refineFinal1_merged_filter_inversions.smap
    <sample_name2>/
        contigs/
            alignmolvref/
                copynumber/
                    cnv_calls_exp.txt
                    cnv_rcmap_exp.txt
            exp_refineFinal1_sv/
                merged_smaps/
                    exp_refineFinal1_merged.xmap
                    exp_refineFinal1_merged_filter_inversions.smap
    ...
```
Curated DIR structure:
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
The default usage of OMKar is batch run all samples within an input data directory (detailed in [Input File Structure](#input-file-structure)). An optional flag `-single` can be used for a single-sample run. In a single sample run, 
please pass the sample-directory instead of the master data-directory as the input directory.

```shell
python3 main.py -dir input_dir -o output_dir [-centro custome_centromere_file] [-single] [-report] [-reportDebug]
```
- `-dir`: Path to the input directory (detailed in [Input File Structure](#input-file-structure)). If using the default batch run, this should be the master data-directory containing each sample-directories. If using the single run flag, this should be the individual sample-directory.
- `-o`: Path to the output directory (detailed in [Output Files](#output-files)).
- `-centro`: (default: hg38 centromere coordinates, 'hg38_centro.txt') a custom centromere coordinate can be used.
- `-single`: Flag for single-sample run. If used together with `-report`, the report will also be generated for the single sample.
- `-report`: Flag to output the HTML report (this takes much longer than the standalone OMKar).
- `-reportDebug`: Flag to output debugging information for the report, including the interpretation information.

As an example, you can use the test files in the `test_files` directory. If you have installed the prerequisites correctly, you can run the following command to see the outputs of OMKar for the test case which contains Balanced translocation between
Chr2 and Chr14, and duplication in Chr2 (intended output in `test_intended_output/`):

```shell
python3 main.py \
  -dir test_input/ \
  -o tets_output/ \
  -report
```
### Output
After running OMKar, the following output files are generated, providing various analyses and visualizations of the karyotype and structural variations. They are organized into subdirectories `OMKar_output/`, `OMKar_report/` (optional), and `logs/` (for debugging).

In `OMKar_output/`:
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

In `OMKar_report/`:
1. **Report Summary (`dashboard.html`)**:
    - Open in a browser, this gives the summary HTML report for all samples ran. The summary includes the tally of SVs on each sample, prediction of disruption on DDG2P genes, and a summary visualization of the karyotype. It also links to the individual sample's HTML report.
2. **Sample Report (`{name}.html`)**:
    - Open in a browser, this gives the HTML report for a particular sample. This report is separated by chromosome clusters. Each cluster includes the ISCN SV list, predicted disrupted DDG2P genes, Molecular Karyotype output, tabulated BED file of the SV calls not incorporated, and a visualization of the chromosome cluster.

### License
MIT License.

### Contact
For issues, questions, or contributions, visit [GitHub Issues](https://github.com/siavashre/OMKar/issues) or contact the maintainers directly.
