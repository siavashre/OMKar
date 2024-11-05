
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

### Output and Interpretation
OMKar provides:
1. **Molecular Karyotype**: A detailed text output listing chromosomal segments.
2. **Graphical Karyotype**: A visual representation of chromosomes with event annotations.
3. **ISCN Event Interpretation**: Events described in ISCN notation, detailing structural variations and chromosomal abnormalities.
4. **Gene List Report**: Important genes near breakpoints are highlighted, aiding in genotype-to-phenotype interpretation.

**Example**:
- A report may identify balanced translocations, aneuploidies, and SV types such as inversions and duplications with ISCN notation and provide corresponding affected gene information.

### Known Issues and Limitations
- **False Positives**: The latest update includes more inversion and duplication-inversion calls, which may introduce some false positives.
- **Acrocentric Region of Chr21**: Recurrent structural variations have been observed in the p-arm region of Chr21. This issue is under investigation.

### Additional Notes
OMKar is actively updated based on simulations and real sample performance. If needed, reset to an earlier version:
```shell
git reset --hard 61a1c62
```
For live updates and issue tracking, please check the GitHub repository.

---

### License
MIT License.

### Contact
For issues, questions, or contributions, visit [GitHub Issues](https://github.com/siavashre/OMKar/issues) or contact the maintainers directly.
