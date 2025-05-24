
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
All dependencies are Python-based for ease of installation.
  - Python 3.7.3 or later.
  - `numpy` >=1.25.0
  - `matplotlib` >=3.8.0
  - `pandas` >=2.1.0
  - `scipy` >=1.11.0
  - `pulp` >=2.7.0
  - `argparse` >=1.4.0
  - `jinja2` >= 3.1.4

We have prepared a dependency file that list all packages and versions.
A convenient installation guide is included in the next section.

### Installation
**Clone OMKar from GitHub and navigate to the working directory:**
```shell
git clone --recurse-submodules https://github.com/siavashre/OMKar.git
```
*The submodules are needed for generating the HTML reports.*

If download speed is limited, you may add `--depth 1` and/or
`--shallow-submodules` to prevents files removed in previous versions from being downloaded (they are for development use only).


**Update OMKar from GitHub (if a previous clone exists):**
```shell
git pull
git submodule update --init --recursive
```

**Set up dependencies:**

You may optionally create a python virtual environment for running OMKar. This may help
with cleaner dependency installation. `python -m venv omkar_env` to create the environment.
`source omkar_env/bin/activate` to launch the environment. To deactivate, run `deactivate`.

Full list of dependencies and version are prepared in `dependencies.txt`. 
Installation of all dependencies can be done by running the following command.
If a python virtual environment is used, first activate the environment, then
run this command:

```shell
pip install -r dependencies.txt
```

**Validate environment:**

To test if installation was correct, we included test files in the `test_files` directory. 
You can run the following command to see the outputs of OMKar for the test case which contains Balanced translocation between
Chr2 and Chr14, and duplication in Chr2:

```shell
python3 main.py -dir test_input/ -o test_output/ -single -report
```

*When running for the first time, it may take about a minute for matplotlib to adjust to the font environment.*

A validation code is also provided to compare the output against the intended output:

```shell
python3 validate_installation.py
```

### Input File Structure
The four required files from Bionano Solve output are the following:
1) cnv_calls_exp.txt
2) cnv_rcmap_exp.txt
3) exp_refineFinal1_merged_filter_inversions_orig.smap
4) exp_refineFinal1_merged.xmap

The input data directory needs to be either the Bionano Solve's default output structure, 
or a curated directory structure. They are illustrated below (with the other un-used files hidden).

Bionano default output DIR structure:
```
<master_data_dir>/
├── <sample1_dir>/
│   └── contigs/
│       ├── alignmolvref/
│       │   └── copynumber/
│       │       ├── cnv_calls_exp.txt
│       │       └── cnv_rcmap_exp.txt
│       └── exp_refineFinal1_sv/
│           └── merged_smaps/
│               ├── exp_refineFinal1_merged_filter_inversions_orig.smap
│               └── exp_refineFinal1_merged.xmap
├── <sample2_dir>/
│   └── contigs/
│       ├── alignmolvref/
│       │   └── copynumber/
│       │       ├── cnv_calls_exp.txt
│       │       └── cnv_rcmap_exp.txt
│       └── exp_refineFinal1_sv/
│           └── merged_smaps/
│               ├── exp_refineFinal1_merged_filter_inversions_orig.smap
│               └── exp_refineFinal1_merged.xmap
...
```
Curated DIR structure:
```
<master_data_dir>/
├── <sample1_dir>/
│   ├── cnv_calls_exp.txt
│   ├── cnv_rcmap_exp.txt
│   ├── exp_refineFinal1_merged_filter_inversions.smap
│   └── exp_refineFinal1_merged.xmap
├── <sample2_dir>/
│   ├── cnv_calls_exp.txt
│   ├── cnv_rcmap_exp.txt
│   ├── exp_refineFinal1_merged_filter_inversions.smap
│   └── exp_refineFinal1_merged.xmap
...
```
## 🛠 Running OMKar with conda
### 🚀 Installation
You can install **OMKar** via `conda` from the [Anaconda Cloud](https://anaconda.org/siavashrei/omkar):
### ✅ Requirements

- `conda` (Miniconda or Anaconda)
- Python 3.10+

### 🔧 Installation (Recommended)

```bash
conda install -c siavashrei omkar
```
### 🧪 Test Installation
```bash
omkar --help
```
You should see the OMKar command-line usage instructions.
## 🛠 Running OMKar with Docker  

OMKar is available as a Docker container, making it easy to run without installing dependencies.  
### 📥 Pull the Docker Image  
To get the latest version of OMKar from Docker Hub, run:  

```bash
docker pull siavashre/omkar
```
### 🚀 Running OMKar in a Container
To run OMKar with an input directory and output file:
```bash
docker run --rm omkar -dir test_input -o test_output -single -report
```
If your "test_input" directory and "test_output" file exist on your local machine, you need to mount them inside the container.

### Usage
The default usage of OMKar is batch run all samples within an input master data directory (detailed in [Input File Structure](#input-file-structure)). 

```shell
python3 main.py -dir input_dir -o output_dir [-centro custome_centromere_file] [-single] [-report] [-noImage] [-reportDebug] 
```
- `-dir`: Path to the input directory (detailed in [Input File Structure](#input-file-structure)). If using the default batch run, this should be the master data-directory containing each sample-directories. If using the single run flag, this should be the individual sample-directory.
- `-o`: Path to the output directory (detailed in [Output](#output)).
- `-centro`: (default: hg38 centromere coordinates, 'hg38_centro.txt') a custom centromere coordinate can be used.
- `-single`: Flag for single-sample run. If used together with `-report`, the report will also be generated for the single sample.
- `-report`: Flag to output the HTML report (this takes much longer than the standalone OMKar).
- `-noImage`: Flag to output the HTML report without any image/visualizations. This saves a lot of time if only the ISCN interpretations are needed.
- `-reportDebug`: Flag to output debugging information for the report, including the interpretation information.

### Output
The output of OMKar is organized into subdirectories `OMKar_output/`, `OMKar_report/` (optional), and `logs/` (for debugging). This is consistent for both batch and single runs.
```
<output_dir>/
├── omkar_output/
│   ├── <sample1>/
│   │   ├── <sample1>.pdf
│   │   ├── <sample1>.txt
│   │   ├── <sample1>_SV.bed
│   │   └── <sample1>_SV.txt
│   ├── <sample2>/
│   │   ├── <sample2>.pdf
│   │   ├── <sample2>.txt
│   │   ├── <sample2>_SV.bed
│   │   └── <sample2>_SV.txt
│   └── ...
└── omkar_report/
│   ├── report_summary.html
│   ├── <sample1>.html
│   ├── <sample2>.html
│   └── ...
└── logs/
    ├── <sample1>.stdout.txt
    ├── <sample2>.stdout.txt
    ├── ...
    └── __report.stdout.txt

```

Of all the output, `{name}.txt` in each `omkar_output/{sample}/` contains the
reconstructed karyotype. To view the report and ISCN interpretation,
`omkar_report/report_summary.html` should be viewed in your browser first, 
where the individual sample's report can be accessed within this HTML page.

Below are the detailed information of each output files, providing various analyses and visualizations of the karyotype and structural variations.

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
1. **Report Summary (`report_summary.html`)**:
    - Open in a browser, this gives the summary HTML report for all samples ran. The summary includes the tally of SVs on each sample, prediction of disruption on DDG2P genes, and a summary visualization of the karyotype. It also links to the individual sample's HTML report.
2. **Sample Report (`{name}.html`)**:
    - Open in a browser, this gives the HTML report for a particular sample. This report is separated by chromosome clusters. Each cluster includes the ISCN SV list, predicted disrupted DDG2P genes, Molecular Karyotype output, tabulated BED file of the SV calls not incorporated, and a visualization of the chromosome cluster.

### Known Bug in Solve v3.7

If you are using the Solve pipeline version 3.7, please be aware of a reported bug. This version may introduce artifacts that report:

1. **Deletion on chromosome 22** in the USP18 region (hg38: 18,181,891–18,754,038).
2. **Inversion on chromosome 7** in the PMS2 region (hg38: 5,986,390–6,519,332).

These artifacts are not observed in Solve versions 3.6. If OMKar reports events in these regions while using Solve v3.7, it is highly likely that these are artifacts. Please interpret such results with caution.

### License
BSD-3-Clause license.

### Contact
For issues, questions, or contributions, visit [GitHub Issues](https://github.com/siavashre/OMKar/issues) or contact the maintainers directly.
