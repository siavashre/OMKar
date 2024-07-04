centromere_file=hg38_centro.txt

data_dir=$1; shift
output_param="default"
report_param="default"
debug_param="default"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --out_dir) output_param="$2"; shift ;;
        --report) report_param="$2"; shift ;;
        --debug) debug_param="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ "$output_param" == "default" ]]; then
    output_param=outputs/
fi

if [[ "$debug_param" == "default" ]]; then
    debug_param=0
else
    debug_param=1
fi

output_dir=${output_param}/omkar_output/
log_dir=${output_param}/omkar_logs/
path_dir=${output_param}/omkar_paths/

export PYTHONPATH=./:$PYTHONPATH
export PYTHONPATH=html_utils/:$PYTHONPATH

mkdir -p ${output_dir}
mkdir -p ${log_dir}
mkdir -p ${path_dir}

for folder in ${data_dir}/*
do
  folder_name=$(basename "$folder")
  mkdir -p ${output_dir}/${folder_name}
  echo $folder_name
  python main.py \
		-cnv ${folder}/cnv_calls_exp.txt \
		-smap ${folder}/exp_refineFinal1_merged_filter_inversions.smap \
		-rcmap ${folder}/cnv_rcmap_exp.txt \
		-xmap ${folder}/exp_refineFinal1_merged.xmap \
		-n $folder_name \
		-o ${output_dir}/${folder_name} \
		-centro ${centromere_file} > ${log_dir}/${folder_name}.log.txt
done

for folder in ${output_dir}/*
do
	for txt_file in ${folder}/*.txt; do
        # Check if the file is not named *.SV.txt
        if [[ ! "$txt_file" =~ SV\.txt$ && ! "$txt_file" =~ preILP_nodes\.txt$ && ! "$txt_file" =~ preILP_edges\.txt$ ]]; then
            # Copy the file to the destination folder
            cp "$txt_file" "$path_dir"
            # echo "Copied $txt_file to $output_dir"
        fi
    done
done

if [[ "$report_param" == "html" ]]; then
    title=$(basename "$data_dir")
    image_dir=${output_param}/report_images/
    python KT_Reports/KT_interpreter_html_report.py \
      --title ${title} \
      --data_dir ${path_dir} \
      --image_dir ${image_dir} \
      --output_dir ${output_param} \
      --debug ${debug_param}
fi
