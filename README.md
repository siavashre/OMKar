# OMKar

#### Requirements

Use the following command in your working directory to clone OMKar:
```shell
git clone https://github.com/siavashre/OMKar.git
```

All dependencies of OMKar are python packages, so the installation should be relatively simple.

Optionally, there are two methods for generating downstream report, either in PDF or in HTML. 
Generating the HTML reports only requires python packages, whereas the PDF version requires a locally installed
TeXWork package. The HTML version currently is less refined looking, but contains the same information.

The input of the batch code takes in a data directory path and output OMKar for all data entries in the directory.
The data directory must be organized in the below hierarchy:
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

#### Usage/Parameters

Run OMKar in the `OMKar/` folder:
```shell
bash batch_run.sh data_dir [--out_dir out_dir] [--report report] [--debug debug]
```

| Argument                         | Type | Description                                                                                            |
|----------------------------------|------|--------------------------------------------------------------------------------------------------------|
| `data_dir`                       | DIR  | Contains the data in the correct hierarchy                                                             |
| (optional) `--out_dir <outdir>`  | DIR  | (default: OMKar/outputs/) Output location                                                              |
| (optional) `--report <html/pdf>` | STR  | (default: None) Generating an HTML/PDF report for all cases ran                                        |
| (optional) `--debug true`        | STR  | (default: false) Only useful when generating a report; include debug information in the report output, |

#### Additional Notes:
We are actively updating OMKar based on Simulation and Real Sample Performance, 
so in general it should be getting better. However, in some situations, we may introduce
errors/bugs in newer release. In the case that you may want to run the initial OMKar release, run the following
in the `OMkar\` folder to reset:
```shell
git reset --hard 61a1c62
```
Check back to this live update thread for issues and notices:
- We are calling more inversions and duplication inversions in the newer version. 
Although this may introduce a few false positives, but based on simulations, we are introducing a lot
more true positives.
- There is an issue where a set of recurrent SVs are appearing in the p-arm (acrocentric) region of Chr21; 
we are investigating
