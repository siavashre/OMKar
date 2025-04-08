import sys
import argparse
from matplotlib import rcParams
from KarReporter.KarUtils.read_OMKar_output import *
from pathlib import Path
rcParams['pdf.fonttype'] = 42
from scripts.run import run_omkar
######################################################################################################################################

def find_input_file_paths(dir):
    """
    Automatically detects input file paths for Bionano Solve or curated directory structures.

    Args:
        dir (str): Path to the master directory containing input files.

    Returns:
        dict: Dictionary with paths to the four required input files:
            - 'cnv': CNV call file path.
            - 'rcmap': RCmap file path.
            - 'xmap': Xmap file path.
            - 'smap': Smap file path.
            - 'dir_structure': Directory structure type ('bionano_default' or 'curated').

    Raises:
        RuntimeError: If the input file structure could not be parsed.
    """
    """
    auto detect between Bionano Solve output structure, OR curated output structure
    :param dir: master dir of the sample containing all input files
    :return: the file paths for the four input files
    """
    cnv_filename = 'cnv_calls_exp.txt'
    rcmap_filename = 'cnv_rcmap_exp.txt'
    xmap_filename = 'exp_refineFinal1_merged.xmap'
    smap_filename = 'exp_refineFinal1_merged_filter_inversions.smap'

    # bionano default filepath
    cnv_filepath = f"{dir}/contigs/alignmolvref/copynumber/{cnv_filename}"
    rcmap_filepath = f"{dir}/contigs/alignmolvref/copynumber/{rcmap_filename}"
    xmap_filepath = f"{dir}/contigs/exp_refineFinal1_sv/merged_smaps/{xmap_filename}"
    smap_filepath = f"{dir}/contigs/exp_refineFinal1_sv/merged_smaps/{smap_filename}"
    if os.path.isfile(cnv_filepath) and os.path.isfile(rcmap_filepath) and os.path.isfile(xmap_filepath) and os.path.isfile(smap_filepath):
        return {'cnv': cnv_filepath,
                'rcmap': rcmap_filepath,
                'xmap': xmap_filepath,
                'smap': smap_filepath,
                'dir_structure': 'bionano_default'}

    # curated dir structure
    if os.path.isfile(f"{dir}/{cnv_filename}") and os.path.isfile(f"{dir}/{rcmap_filename}") and os.path.isfile(f"{dir}/{xmap_filename}") and os.path.isfile(f"{dir}/{smap_filename}"):
        cnv_filepath, rcmap_filepath, xmap_filepath, smap_filepath = f"{dir}/{cnv_filename}", f"{dir}/{rcmap_filename}", f"{dir}/{xmap_filename}", f"{dir}/{smap_filename}"
        return {'cnv': cnv_filepath,
                'rcmap': rcmap_filepath,
                'xmap': xmap_filepath,
                'smap': smap_filepath,
                'dir_structure': 'curated'}

    # debug use
    bionano_status1 = os.path.isfile(f'{dir}/contigs/alignmolvref/copynumber/{cnv_filename}')
    bionano_status2 = os.path.isfile(f'{dir}/contigs/alignmolvref/copynumber/{rcmap_filename}')
    bionano_status3 = os.path.isfile(f'{dir}/contigs/exp_refineFinal1_sv/merged_smaps/{xmap_filename}')
    bionano_status4 = os.path.isfile(f'{dir}/contigs/exp_refineFinal1_sv/merged_smaps/{smap_filename}')
    curated_status1 = os.path.isfile(f'{dir}/{cnv_filename}')
    curated_status2 = os.path.isfile(f'{dir}/{rcmap_filename}')
    curated_status3 = os.path.isfile(f'{dir}/{xmap_filename}')
    curated_status4 = os.path.isfile(f'{dir}/{smap_filename}')
    print(f"attempted {dir}/contigs/alignmolvref/copynumber/{cnv_filename}: {bionano_status1}")
    print(f"attempted {dir}/contigs/alignmolvref/copynumber/{rcmap_filename}: {bionano_status2}")
    print(f"attempted {dir}/contigs/exp_refineFinal1_sv/merged_smaps/{xmap_filename}: {bionano_status3}")
    print(f"attempted {dir}/contigs/exp_refineFinal1_sv/merged_smaps/{smap_filename}: {bionano_status4}")
    print(f"attempted {dir}/{cnv_filename}: {curated_status1}")
    print(f"attempted {dir}/{rcmap_filename}: {curated_status2}")
    print(f"attempted {dir}/{xmap_filename}: {curated_status3}")
    print(f"attempted {dir}/{smap_filename}: {curated_status4}")

    raise RuntimeError(f'file structure could not be parsed: {dir}')

def main():
    """
    Main function for running the OMKar pipeline.

    Parses command-line arguments to determine the input directory, centromere coordinates,
    output directory, and additional flags for single-sample processing, report generation,
    image generation, and debug mode.

    Workflow:
        - Detects input file paths.
        - Executes the OMKar pipeline for either a single sample or multiple samples.
        - Optionally generates an HTML report summarizing the results.

    Command-line Arguments:
        -dir (str): Path to the directory containing input data (required).
        -centro (str): Path to the file containing centromere coordinates (default: 'hg38_centro.txt').
        -o (str): Path to the output directory (required).
        -single (flag): Process only a single sample.
        -report (flag): Generate an HTML report.
        -noImage (flag): Disable image generation in the report.
        -reportDebug (flag): Enable debug information in the HTML report.

    Returns:
        None: Outputs are saved to the specified output directory.
    """
    repo_dir = os.path.dirname(os.path.abspath(__file__))
    default_stdout = sys.stdout
    default_centro_path = f'{repo_dir}/hg38_centro.txt'

    parser = argparse.ArgumentParser()
    parser.add_argument("-dir", "--dir", help="path to directory containing input data", required=True)
    parser.add_argument("-centro", "--centro", default=default_centro_path, help="path to file contains centromere coordinates", required=False)
    parser.add_argument("-o", "--output", help="path to output dir", required=True)
    parser.add_argument("-single", "--single", help="flag to singal only run on one sample", action='store_true')
    parser.add_argument("-report", "--report", help="flag to generate html report", action='store_true')
    parser.add_argument("-noImage", "--noImage", help="flag to disable image generation", action='store_true')
    parser.add_argument("-reportDebug", "--reportDebug", help="generate report debug in html report", action='store_true')
    args = parser.parse_args()

    os.makedirs(f'{args.output}/omkar_output/', exist_ok=True)
    os.makedirs(f'{args.output}/logs/', exist_ok=True)
    if args.single:
        filepaths = find_input_file_paths(args.dir)
        sample_name = os.path.basename(os.path.normpath(args.dir))
        print(f'running sample: {sample_name}')
        sys.stdout = open(f"{args.output}/logs/{sample_name}.stdout.txt", 'w')
        sample_output_dir = f"{args.output}/omkar_output/"
        run_omkar(filepaths['cnv'], filepaths['smap'], filepaths['rcmap'], filepaths['xmap'], args.centro, sample_name, sample_output_dir)
        sys.stdout.close()
        sys.stdout = default_stdout
    else:
        for sample_dir in os.listdir(args.dir):
            if not os.path.isdir(os.path.join(args.dir, sample_dir)):
                continue  # this ignores non-dir files such as .DS_Store
            filepaths = find_input_file_paths(f"{args.dir}/{sample_dir}/")
            sample_name = sample_dir
            print(f'running sample: {sample_name}')
            sys.stdout = open(f"{args.output}/logs/{sample_name}.stdout.txt", 'w')
            sample_output_dir = f"{args.output}/omkar_output/"
            run_omkar(filepaths['cnv'], filepaths['smap'], filepaths['rcmap'], filepaths['xmap'], args.centro, sample_name, sample_output_dir)
            sys.stdout.close()
            sys.stdout = default_stdout

    if args.report:
        from KarReporter import generate_html_report
        generate_image = not args.noImage
        sys.stdout = open(f"{args.output}/logs/__report.stdout.txt", 'w')
        os.makedirs(f'{args.output}/omkar_report/', exist_ok=True)
        omkar_input_dir_for_report = args.dir if not args.single else str(Path(args.dir).parent)
        generate_html_report(generate_image,
                             None,
                             '',
                             f'{args.output}/omkar_output/',
                             f'{args.output}/omkar_report/clustered_karyotype_images/',
                             f'{args.output}/omkar_report/',
                             omkar_input_dir_for_report,
                             debug=args.reportDebug)
        sys.stdout.close()
        sys.stdout = default_stdout
if __name__ == "__main__":
    main()
