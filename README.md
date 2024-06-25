# OMKar

#### Requirements

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

```shell
bash batch_run.sh data_dir [--out_dir out_dir] [--report report] [--debug debug]
```

| Argument                         | Type | Description                                                                                            |
|----------------------------------|------|--------------------------------------------------------------------------------------------------------|
| `data_dir`                       | DIR  | Contains the data in the correct hierarchy                                                             |
| (optional) `--out_dir <outdir>`  | DIR  | (default: OMKar/outputs/) Output location                                                              |
| (optional) `--report <html/pdf>` | STR  | (default: None) Generating an HTML/PDF report for all cases ran                                        |
| (optional) `--debug true`        | STR  | (default: false) Only useful when generating a report; include debug information in the report output, |
