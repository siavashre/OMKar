
# OMKar

OMKar is a computational tool for automated karyotyping using Optical Genome Mapping (OGM) data, focusing on detecting structural variations (SV) and copy number variations (CNV) with high resolution. OMKar provides a virtual karyotype and detailed analysis of chromosomal abnormalities, making it ideal for identifying constitutional genetic disorders.

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
