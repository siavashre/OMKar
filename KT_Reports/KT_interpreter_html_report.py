from KT_visualizer import *
from jinja2 import Environment, FileSystemLoader
from read_OMKar_output import *
import os
import base64
import re

def image_to_base64(image_path):
    try:
        with open(image_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    except FileNotFoundError:
        print(f"Error: File {image_path} not found.")
        return ""


def batch_populate_contents(omkar_output_dir, image_dir, file_of_interest=None, compile_image=False, debug=False):
    headers = []
    cases_with_events = []
    image_paths = []
    iscn_reports = []
    genes_reports = []
    debug_outputs = []  # list of dicts [{'segs': [], 'mt_haps': [], 'wt_haps': []}]
    files = [file for file in os.listdir(omkar_output_dir)]
    for file in files:
        if file_of_interest is not None:
            if file not in file_of_interest:
                continue
        filename = file.split('.')[0]
        file_path = omkar_output_dir + file
        print(file)
        mt_indexed_lists, mt_path_chrs, segment_to_index_dict, segment_size_dict = read_OMKar_to_indexed_list(file_path, forbidden_region_file)
        index_to_segment_dict = reverse_dict(segment_to_index_dict)
        mt_path_chrs = [info.split(': ')[-1] for info in mt_path_chrs]
        wt_path_dict = generate_wt_from_OMKar_output(segment_to_index_dict)
        wt_indexed_lists = populate_wt_indexed_lists(mt_path_chrs, wt_path_dict)
        events, aligned_haplotypes = interpret_haplotypes(mt_indexed_lists, wt_indexed_lists, mt_path_chrs, segment_size_dict)
        if len(events) == 0:
            continue
        else:
            cases_with_events.append(filename)
        dependent_clusters, cluster_events = form_dependent_clusters(events, aligned_haplotypes, index_to_segment_dict)
        print(dependent_clusters)
        ## iterate over all clusters
        n_clusters = len(dependent_clusters)
        for image_cluster_idx, (c_cluster, c_events) in enumerate(zip(dependent_clusters, cluster_events)):
            # to remove all later file names, check cluster_idx != 0
            headers.append('{}: cluster {} (out of {})'.format(filename, image_cluster_idx + 1, n_clusters))
            ## include all homologues
            event_chr = set()
            for cluster_idx in c_cluster:
                event_chr.add(aligned_haplotypes[cluster_idx].chrom)
            hap_idx_to_plot = []
            for hap_idx, hap in enumerate(aligned_haplotypes):
                if hap.chrom in event_chr:
                    hap_idx_to_plot.append(hap_idx)

            c_aligned_haplotypes = [aligned_haplotypes[i] for i in hap_idx_to_plot]

            ## generate report text
            c_events = sort_events(c_events)
            iscn_events, genes_report = format_report(c_events, aligned_haplotypes, index_to_segment_dict, debug=debug)
            ## generate image
            c_vis_input = generate_visualizer_input(c_events, c_aligned_haplotypes, segment_to_index_dict)

            def vis_key(input_vis):
                chr_val = input_vis['chr'][3:]
                if chr_val == "X":
                    return_val = 23.0
                elif chr_val == "Y":
                    return_val = 24.0
                else:
                    return_val = float(chr_val)
                if input_vis['highlight']:
                    return_val += 0.5  # highlight always later
                return return_val

            c_vis_input = sorted(c_vis_input, key=vis_key)
            image_prefix = "{}/{}_imagecluster{}".format(image_dir, filename, image_cluster_idx)
            image_path = image_prefix + '_rotated.png'
            relative_image_path = image_dir.replace('latex_reports/', '') + image_path.split('/')[-1]
            if compile_image:
                if len(c_vis_input) <= 4:
                    make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix, IMG_LENGTH_SCALE_VERTICAL_SPLIT)
                else:
                    make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix, IMG_LENGTH_SCALE_HORIZONTAL_SPLIT)

            image_paths.append(relative_image_path)
            iscn_reports.append(iscn_events)
            genes_reports.append(genes_report)

            ## generate debug output
            debug_segs = set()
            debug_mt_haps = []
            debug_wt_haps = []
            debug_hap_ids = []
            debug_mt_aligned = []
            debug_wt_aligned = []
            for aligned_haplotype in c_aligned_haplotypes:
                unique_segs = aligned_haplotype.unique_segment_indices()
                for seg in unique_segs:
                    seg_object = index_to_segment_dict[int(seg)]
                    debug_segs.add((seg, seg_object.chr_name, f"{seg_object.start:,}", f"{seg_object.end:,}", f"{len(seg_object):,}"))
            debug_segs = list(debug_segs)
            debug_segs = sorted(debug_segs, key=lambda x: int(x[0]))

            for c_vis in c_vis_input:
                debug_hap_ids.append(c_vis['hap_id'])
                # find the correct aligned_haplotype for mt/wt
                hap_found = False
                for aligned_haplotype in aligned_haplotypes:
                    if aligned_haplotype.id == c_vis['hap_id']:
                        debug_mt_haps.append(aligned_haplotype.mt_hap)
                        debug_wt_haps.append(aligned_haplotype.wt_hap)
                        debug_mt_aligned.append(aligned_haplotype.mt_aligned)
                        debug_wt_aligned.append(aligned_haplotype.wt_aligned)
                        hap_found = True
                        break
                if not hap_found:
                    raise RuntimeError('hap not found')
            debug_outputs.append({'segs': debug_segs, 'mt_haps': debug_mt_haps, 'wt_haps': debug_wt_haps, 'IDs': debug_hap_ids,
                                  'mt_aligned': debug_mt_aligned, 'wt_aligned': debug_wt_aligned})
    return headers, cases_with_events, image_paths, iscn_reports, genes_reports, debug_outputs


def html_hyperlink_coordinates(input_str, proximity=50000):
    return_dict = {}  # {replacement_string: hyperlinked_string}

    pattern = r'Chr(\d+|X|Y): (\d{1,3}(?:,\d{3})*)-(\d{1,3}(?:,\d{3})*) \(.*?\)'
    matches_itr = re.finditer(pattern, input_str)
    for match in matches_itr:
        replacement_str = input_str[match.start(): match.end()]
        start_pos = int(match.group(2).replace(',', ''))
        end_pos = int(match.group(3).replace(',', ''))
        ucsc_url = get_ucsc_url('chr' + match.group(1), start_pos, end_pos)
        hyperlinked_str = f'<a href="{ucsc_url}">{replacement_str}</a>'
        return_dict[replacement_str] = hyperlinked_str

    pattern = r'Chr(\d+): (\d{1,3}(?:,\d{3})*) \(.*?\)'
    matches_itr = re.finditer(pattern, input_str)
    for match in matches_itr:
        replacement_str = input_str[match.start(): match.end()]
        chrom = 'chr' + match.group(1)
        pos = int(match.group(2).replace(',', ''))
        c_chr_length = get_chr_length_from_forbidden_file(chrom)
        ucsc_url = get_ucsc_url(chrom, max(0, pos - proximity), min(c_chr_length, pos + proximity))
        hyperlinked_str = f'<a href="{ucsc_url}">{replacement_str}</a>'
        return_dict[replacement_str] = hyperlinked_str

    return return_dict


def hyperlink_iscn_interpretation(input_str):
    hyperlinked_mapping = html_hyperlink_coordinates(input_str)
    for replacement_str, hyperlinked_str in hyperlinked_mapping.items():
        input_str = input_str.replace(replacement_str, hyperlinked_str)
    return input_str


def generate_report(compile_image, cases_of_interest, title, omkar_output_dir, image_output_dir, output_file, debug=False):
    os.makedirs(image_output_dir, exist_ok=True)

    # one tuple per cluster (event)
    headers, cases_with_events, image_paths, iscn_reports, genes_reports, debug_outputs = batch_populate_contents(omkar_output_dir, image_output_dir,
                                                                                          file_of_interest=cases_of_interest, compile_image=compile_image, debug=debug)
    images_base64 = [image_to_base64(img) for img in image_paths]

    formatted_genes_reports = [format_genes_report(genes_report) for genes_report in genes_reports]
    columns_order = ['SV', 'rationale', 'gene name', 'gene omim']

    ## hyperlinking
    for iscn_report in iscn_reports:
        for iscn_report_idx, (iscn, sv_interpretation) in enumerate(iscn_report):
            hyperlinked_sv_interpretation = hyperlink_iscn_interpretation(sv_interpretation)
            iscn_report[iscn_report_idx][1] = hyperlinked_sv_interpretation

    content = [(header, text, image, table_content, debug) for header, text, image, table_content, debug in zip(headers, iscn_reports, images_base64, formatted_genes_reports, debug_outputs)]

    env = Environment(loader=FileSystemLoader('KT_Reports/html_utils/'))
    template = env.get_template('template.html')
    rendered_html = template.render(title=title, content=content, columns_order=columns_order, debug=debug)

    with open(output_file, 'w') as f:
        f.write(rendered_html)
    print(f"HTML file generated: {os.path.abspath(output_file)}")


def manual_test():
    # Define the data
    title = "My Text and Images"
    texts = ['text1',
            'text2',
            'text3']

    ## ZJ: image paths need to be relative path
    images = ['/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_plots_new/3_imagecluster0_rotated.png',
              '/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_plots_new/3_imagecluster1_rotated.png',
              '/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_plots_new/3_imagecluster2_rotated.png']
    images_base64 = []
    for img in images:
        images_base64.append(image_to_base64(img))


    table_contents = [
        [["Row1-Col1", "Row1-Col2"], ["Row2-Col1", "Row2-Col2"], ["Row3-Col1", "Row3-Col2"], ["Row4-Col1", "Row4-Col2"],
         ["Row5-Col1", "Row5-Col2"], ["Row6-Col1", "Row6-Col2"], ["Row7-Col1", "Row7-Col2"], ["Row8-Col1", "Row8-Col2"],
         ["Row9-Col1", "Row9-Col2"], ["Row10-Col1", "Row10-Col2"], ["Row11-Col1", "Row11-Col2"], ["Row12-Col1", "Row12-Col2"]],

        [["Row1-Col1", "Row1-Col2"], ["Row2-Col1", "Row2-Col2"], ["Row3-Col1", "Row3-Col2"], ["Row4-Col1", "Row4-Col2"],
         ["Row5-Col1", "Row5-Col2"], ["Row6-Col1", "Row6-Col2"], ["Row7-Col1", "Row7-Col2"], ["Row8-Col1", "Row8-Col2"],
         ["Row9-Col1", "Row9-Col2"], ["Row10-Col1", "Row10-Col2"], ["Row11-Col1", "Row11-Col2"], ["Row12-Col1", "Row12-Col2"]],

        [["Row1-Col1", "Row1-Col2"], ["Row2-Col1", "Row2-Col2"], ["Row3-Col1", "Row3-Col2"], ["Row4-Col1", "Row4-Col2"],
         ["Row5-Col1", "Row5-Col2"], ["Row6-Col1", "Row6-Col2"], ["Row7-Col1", "Row7-Col2"], ["Row8-Col1", "Row8-Col2"],
         ["Row9-Col1", "Row9-Col2"], ["Row10-Col1", "Row10-Col2"], ["Row11-Col1", "Row11-Col2"], ["Row12-Col1", "Row12-Col2"]]
    ]

    content = [(text, image, table_content) for text, image, table_content in zip(texts, images_base64, table_contents)]

    # Create an environment for Jinja2
    env = Environment(loader=FileSystemLoader('html_utils/'))
    template = env.get_template('template.html')

    # Render the template with the data
    rendered_html = template.render(title=title, content=content)

    # Write the rendered HTML to a file
    output_file = 'html_utils/test.html'
    with open(output_file, 'w') as f:
        f.write(rendered_html)

    print(f"HTML file generated: {os.path.abspath(output_file)}")


if __name__ == "__main__":
    forbidden_region_file = "KT_Reports/Metadata/acrocentric_telo_cen.bed"
    # i_title, i_omkar_output_dir, i_image_output_dir, i_output_file = 'keyhole', 'real_case_data/keyhole_OMKar_output_paths/', 'html_utils/keyhole_plots/', 'html_utils/keyhole.html'
    # i_title, i_omkar_output_dir, i_image_output_dir, i_output_file = 'sunnyside', 'real_case_data/sunnyside_OMKar_output_paths/', 'html_utils/sunnyside_plots/', 'html_utils/sunnyside.html'
    # i_title, i_omkar_output_dir, i_image_output_dir, i_output_file = 'dremsek', 'real_case_data/dremsek_OMKar_output_paths/', 'html_utils/dremsek_plots/', 'html_utils/dremsek.html'
    # i_title, i_omkar_output_dir, i_image_output_dir, i_output_file = 'karsim', 'omkar_analyses_pipeline/builds/b14/omkar_paths/', 'html_utils/karsim_plots/', 'html_utils/karsim.html'
    # i_title, i_omkar_output_dir, i_image_output_dir, i_output_file = 'test', 'omkar_analyses_pipeline/builds/b14/omkar_paths/', 'html_utils/test_plots/', 'html_utils/test.html'
    # test(True, ['23Y_CMT1A_r2.1.txt'], i_title, i_omkar_output_dir, i_image_output_dir, i_output_file, debug=True)

    # i_title, i_omkar_output_dir, i_image_output_dir, i_output_file = 'sunnyside', 'real_case_data/sunnyside_OMKar_output_paths/', 'html_utils/sunnyside_plots/', 'html_utils/sunnyside.html'
    # test(True, None, i_title, i_omkar_output_dir, i_image_output_dir, i_output_file, debug=True)

    #################################OMKAR INTEGRATION####################
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--title', default='Title', type=str, help='title of report')
    parser.add_argument('--data_dir', default='outputs/omkar_paths/', type=str, help='path to DIR of OMKar paths output')
    parser.add_argument('--output_dir', default='outputs/', type=str, help='path to DIR of report output')
    parser.add_argument('--image_dir', default='outputs/report_images/', type=str, help='path to DIR for temporary image storage')
    parser.add_argument('--debug', default=0, type=int, help='output debug info in report')
    args = parser.parse_args()
    i_title, i_omkar_output_dir, i_report_output_dir, i_image_output_dir, debug = args.title, args.data_dir, args.output_dir, args.image_dir, args.debug
    i_output_file = '{}/{}.html'.format(i_report_output_dir, i_title)
    generate_report(True, None, i_title, i_omkar_output_dir, i_image_output_dir, i_output_file, debug=debug)
