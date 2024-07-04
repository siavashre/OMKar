import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import colorsys
from PIL import Image

from KT_interpreter import *
from utill import *
from Structures import *


def reduce_saturation(color, factor):
    """Reduce the saturation of a color by a given factor (0-1)."""
    rgb = mcolors.to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(*rgb)
    s *= factor
    return colorsys.hls_to_rgb(h, l, s)


def get_text_color(bg_color):
    """Determine if black or white text would be more readable on the given background color."""
    rgb = mcolors.to_rgb(bg_color)
    luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return 'black' if luminance > 0.5 else 'white'


# constant parameters, x and y based on the unrotated orientation
IMAGE_DPI = 500
IMG_LENGTH_SCALE_VERTICAL_SPLIT = 0.92
IMG_LENGTH_SCALE_HORIZONTAL_SPLIT = 0.55
MAX_CHR_LEN_IF_NO_SCALE = 250
SV_LABEL_MIN_DISTANCE = 5
MAX_LABEL_EACH_LINE = 2

MIN_LENGTH_ARROW_WITHOUT_SCALE = 1.5  # smallest orientation contig to have arrow labels
MIN_LENGTH_SHOW_ORIGIN_NAME = 2.5
ORIENTATION_ARROW_X_OFFSET = 0.15
ORIENTATION_COLOR = 'black'
ORIENTATION_ALPHA = 0.7
ORIENTATION_BAR_WEIGHT = 0.5

WHOLE_CHR_Y_OFFSET = 2
CHR_HEADER_Y_OFFSET = -0.2 + WHOLE_CHR_Y_OFFSET
CHR_BAND_Y_OFFSET = WHOLE_CHR_Y_OFFSET
CHR_BAND_MARK_Y_OFFSET = 0.05 + CHR_BAND_Y_OFFSET
CHR_HEADER_HIGHLIGHT_Y_OFFSET = -0.85 + WHOLE_CHR_Y_OFFSET
ORIGIN_Y_OFFSET = 1.3 + WHOLE_CHR_Y_OFFSET
TICK_Y_OFFSET = -0.05 + WHOLE_CHR_Y_OFFSET
TICK_MARKING_Y_OFFSET = -0.3 + WHOLE_CHR_Y_OFFSET
LABEL_BAR_Y_OFFSET = -0.05 + WHOLE_CHR_Y_OFFSET
LABEL_MARK_Y_OFFSET = -2.05 + WHOLE_CHR_Y_OFFSET

CHR_HEADER_X_OFFSET = 10
CHR_HEADER_HIGHLIGHT_COLOR = 'red'
CHR_HEADER_HIGHLIGHT_X_OFFSET = -0.03
BAND_WIDTH = 1
ORIGIN_WIDTH = 0.5
ORIGIN_MARK_Y_OFFSET = 0.05

BAND_SATURATION = 1
BAND_ALPHA = 0.85
BAND_TEXT_WEIGHT = 'normal'
ORIGIN_ALPHA = 0.7

TICK_MARKING_X_OFFSET = 0.25  # helps to center the tickmarkings
TICK_LEN = 0.2
TICK_THICKNESS = 0.9
TICK_ALPHA = 0.65
TICK_MARKING_ALPHA = 0.85
SUBTICK_LEN = 0.065
SUBTICK_THICKNESS = 0.9
SUBTICK_ALPHA = 0.5

LABEL_BAR_LEN1 = 0.95
LABEL_BAR_LEN2 = 0.7
LABEL_BAR_THICKNESS = 1.25
LABEL_BAR_ALPHA = 1
LABEL_MARK_ALPHA = 1

BAND_FONTSIZE = 3
LABEL_MARK_FONTSIZE = 6
CHR_HEADER_HIGHLIGHT_FONTSIZE = 12
CHR_HEADER_FONTSIZE = 7
SCALE_MARKING_FONTSIZE = 5
ORIGIN_FONTSIZE = 5
LABEL_MARK_COLOR = 'red'

BAND_RECT_LINEWIDTH = 0.5
ORIGIN_RECT_LINEWIDTH = 0.5

# constant parameters: do not adjust, dependent to above
TICK_END_Y_OFFSET = TICK_Y_OFFSET - TICK_LEN
SUBTICK_END_Y_OFFSET = TICK_Y_OFFSET - SUBTICK_LEN
LABEL_BAR_END_Y_OFFSET1 = LABEL_BAR_Y_OFFSET - LABEL_BAR_LEN1
LABEL_BAR_END_Y_OFFSET2 = LABEL_BAR_Y_OFFSET - LABEL_BAR_LEN2

color_mapping = {
    'gneg': 'white',
    'gpos25': '#C0C0C0',  # light grey
    'gpos50': '#808080',  # medium grey
    'gpos75': '#404040',  # dark grey
    'gpos100': 'black',  # full black
    'acen': 'red',  # centromere
    'gvar': 'blue',  # variable region
    'stalk': '#87CEEB'  # light blue (skyblue)
}
chr_color_mapping = {
    '1': '#73a9a3',
    '2': '#ffffd9',
    '3': '#a5a4c6',
    '4': '#fb8072',
    '5': '#73a9c3',
    '6': '#d6a780',
    '7': '#a8d8b4',
    '8': '#fcebf5',
    '9': '#d9d9d9',
    '10': '#b780bd',
    '11': '#c5e7c5',
    '12': '#ffed9f',
    '13': '#ff6a6a',
    '14': '#73a8d3',
    '15': '#78af78',
    '16': '#a94ea3',
    '17': '#ff7f40',
    '18': '#ffff99',
    '19': '#a56e40',
    '20': '#f7a1bf',
    '21': '#999999',
    '22': '#66c2a5',
    'X': '#fc8d62',
    'Y': '#8da0cb'
}

# Reduce saturation of the colors in the color mapping
reduced_saturation_mapping = {k: reduce_saturation(v, BAND_SATURATION) for k, v in color_mapping.items()}


def plot_chromosome(ax, chromosome_data, y_offset, x_offset, len_scaling):
    ## Chrom header
    ax.text(x_offset, y_offset + CHR_HEADER_Y_OFFSET, chromosome_data['name'],
            va='bottom', fontsize=CHR_HEADER_FONTSIZE, rotation=90, weight='bold')

    ## Bands
    for band in chromosome_data['bands']:
        start = band['start']
        end = band['end']
        name = band['band']
        stain = band['stain']
        color = reduced_saturation_mapping[stain]
        text_color = get_text_color(color)

        chrom_bands = patches.Rectangle((x_offset + start + CHR_HEADER_X_OFFSET, y_offset + CHR_BAND_Y_OFFSET), end - start, BAND_WIDTH,
                                        linewidth=1, edgecolor='black', facecolor=color, alpha=BAND_ALPHA, lw=BAND_RECT_LINEWIDTH)
        ax.add_patch(chrom_bands)
        if end - start > (1 / len_scaling):
            # do not label band that are too narrow
            ax.text(x_offset + (start + end) / 2 + CHR_HEADER_X_OFFSET, y_offset + BAND_WIDTH / 2 + CHR_BAND_MARK_Y_OFFSET, name,
                    ha='center', va='center', fontsize=BAND_FONTSIZE, color=text_color, rotation=90, weight=BAND_TEXT_WEIGHT)

    # ## Origins
    # for origin in chromosome_data['origins']:
    #     start = origin['start']
    #     end = origin['end']
    #     name = origin['name']
    #     color = chr_color_mapping[name]
    #     text_color = get_text_color(color)
    #     origin_bands = patches.Rectangle((x_offset + start + CHR_HEADER_X_OFFSET, y_offset + ORIGIN_Y_OFFSET), end - start, ORIGIN_WIDTH,
    #                                      linewidth=1, edgecolor='black', facecolor=color, alpha=ORIGIN_ALPHA, lw=ORIGIN_RECT_LINEWIDTH)
    #     ax.add_patch(origin_bands)
    #     ax.text(x_offset + (start + end) / 2 + CHR_HEADER_X_OFFSET, y_offset + ORIGIN_Y_OFFSET + ORIGIN_WIDTH / 2 + ORIGIN_MARK_Y_OFFSET, name,
    #             ha='center', va='center', fontsize=ORIGIN_FONTSIZE, color=text_color, rotation=90, weight=BAND_TEXT_WEIGHT)

    ## Orientations Contig
    for contig in chromosome_data['orientation_contigs']:
        origin_color = chr_color_mapping[contig['origin']]
        if contig['length'] <= MIN_LENGTH_ARROW_WITHOUT_SCALE:
            # triangle
            x1 = x_offset + contig['start'] + CHR_HEADER_X_OFFSET
            x2 = x_offset + contig['end'] + CHR_HEADER_X_OFFSET
            y1 = y_offset + ORIGIN_Y_OFFSET + ORIGIN_WIDTH / 2
            y2a = y_offset + ORIGIN_Y_OFFSET
            y2b = y_offset + ORIGIN_Y_OFFSET + ORIGIN_WIDTH
            if not contig['orientation']:
                # up (left) arrow
                x1 += ORIENTATION_ARROW_X_OFFSET
                vertices = [(x1, y1), (x2, y2a), (x2, y2b)]
            else:
                # down (right) arrow
                x2 -= ORIENTATION_ARROW_X_OFFSET
                vertices = [(x1, y2a), (x1, y2b), (x2, y1)]
        else:
            # pentagon
            vertex_y = y_offset + ORIGIN_Y_OFFSET + ORIGIN_WIDTH / 2
            bottom_y = y_offset + ORIGIN_Y_OFFSET
            top_y = y_offset + ORIGIN_Y_OFFSET + ORIGIN_WIDTH
            if not contig['orientation']:
                vertex_x = x_offset + contig['start'] + CHR_HEADER_X_OFFSET
                x2 = vertex_x + MIN_LENGTH_ARROW_WITHOUT_SCALE
                x3 = x_offset + contig['end'] + CHR_HEADER_X_OFFSET
                vertices = [(vertex_x, vertex_y), (x2, bottom_y), (x3, bottom_y), (x3, top_y), (x2, top_y)]
            else:
                vertex_x = x_offset + contig['end'] + CHR_HEADER_X_OFFSET
                x1 = x_offset + contig['start'] + CHR_HEADER_X_OFFSET
                x2 = vertex_x - MIN_LENGTH_ARROW_WITHOUT_SCALE
                vertices = [(x1, top_y), (x1, bottom_y), (x2, bottom_y), (vertex_x, vertex_y), (x2, top_y)]

        # IDK why it is saying the vertices is of the wrong type
        patch = patches.Polygon(vertices, closed=True, edgecolor='black', facecolor=origin_color, alpha=ORIGIN_ALPHA, lw=ORIGIN_RECT_LINEWIDTH)
        ax.add_patch(patch)

        if contig['length'] > MIN_LENGTH_SHOW_ORIGIN_NAME:
            x = x_offset + (contig['start'] + contig['end']) / 2 + CHR_HEADER_X_OFFSET
            y = y_offset + ORIGIN_Y_OFFSET + ORIGIN_WIDTH / 2 + ORIGIN_MARK_Y_OFFSET
            text_color = get_text_color(origin_color)
            ax.text(x, y, contig['origin'],
                    ha='center', va='center', fontsize=ORIGIN_FONTSIZE, color=text_color, rotation=90, weight=BAND_TEXT_WEIGHT)

    ## Modified Chrom's header
    if chromosome_data['highlight']:
        ax.text(x_offset + CHR_HEADER_HIGHLIGHT_X_OFFSET, y_offset + CHR_HEADER_HIGHLIGHT_Y_OFFSET, "*",
                va='bottom', ha='left',
                fontsize=CHR_HEADER_HIGHLIGHT_FONTSIZE, rotation=90, weight='bold', color=CHR_HEADER_HIGHLIGHT_COLOR)

    ## Add sub-scale ticks
    for i in range(0, math.floor(chromosome_data['length'] / len_scaling) + 1, 2):
        subtick_x_loc = x_offset + (i * len_scaling) + CHR_HEADER_X_OFFSET
        subtick_y_start_loc = y_offset + TICK_Y_OFFSET
        subtick_y_end_loc = y_offset + SUBTICK_END_Y_OFFSET
        ax.plot([subtick_x_loc, subtick_x_loc], [subtick_y_start_loc, subtick_y_end_loc],
                color='grey', linewidth=SUBTICK_THICKNESS, alpha=SUBTICK_ALPHA)
    ## Add scale ticks
    for i in range(0, math.floor(chromosome_data['length'] / len_scaling) + 1, 10):
        tick_x_loc = x_offset + (i * len_scaling) + CHR_HEADER_X_OFFSET
        tick_y_start_loc = y_offset + TICK_Y_OFFSET - SUBTICK_LEN  # to not overlap with the subticks
        tick_y_end_loc = y_offset + TICK_END_Y_OFFSET
        tickmark_x_loc = x_offset + (i * len_scaling) + CHR_HEADER_X_OFFSET + TICK_MARKING_X_OFFSET
        tickmark_y_loc = y_offset + TICK_MARKING_Y_OFFSET
        ax.plot([tick_x_loc, tick_x_loc], [tick_y_start_loc, tick_y_end_loc],
                color='red', linewidth=TICK_THICKNESS, alpha=TICK_ALPHA)
        ax.text(tickmark_x_loc, tickmark_y_loc, str(i),
                ha='center', va='top', fontsize=SCALE_MARKING_FONTSIZE, rotation=90, alpha=TICK_MARKING_ALPHA)

    ## sv_labels
    for sv_label in chromosome_data['sv_labels']:
        pos = sv_label['pos']

        ## format label-text
        labels = []
        for label_idx, label_itr in enumerate(sv_label['label']):
            label_header = label_itr.split(']')[0] + ']'
            if label_idx != len(sv_label['label']) - 1 and (label_idx + 1) % MAX_LABEL_EACH_LINE == 0:
                label_header += '\n'
            labels.append(label_header)
        label = ''.join(labels)

        label_bar_x_loc = x_offset + pos + CHR_HEADER_X_OFFSET
        label_bar_y_start_loc = y_offset + LABEL_BAR_Y_OFFSET
        if len(labels) > 1:
            label_bar_y_end_loc = y_offset + LABEL_BAR_END_Y_OFFSET2 - SUBTICK_LEN
        else:
            label_bar_y_end_loc = y_offset + LABEL_BAR_END_Y_OFFSET1 - SUBTICK_LEN
        label_mark_y_loc = label_bar_y_end_loc + LABEL_MARK_Y_OFFSET
        ax.plot([label_bar_x_loc, label_bar_x_loc], [label_bar_y_start_loc, label_bar_y_end_loc],
                color=LABEL_MARK_COLOR, linewidth=LABEL_BAR_THICKNESS, alpha=LABEL_BAR_ALPHA)
        ax.text(label_bar_x_loc, label_mark_y_loc, label,
                ha='center', va='top', fontsize=LABEL_MARK_FONTSIZE, rotation=90,
                alpha=LABEL_MARK_ALPHA, color=LABEL_MARK_COLOR, weight='normal')


def rotate_image(input_image_path, output_image_path):
    # Open an image file
    with Image.open(input_image_path) as img:
        # Rotate the image by 90 degrees
        rotated_img = img.rotate(270, expand=True)
        # Save the rotated image
        rotated_img.save(output_image_path)


def generate_visualizer_input(events, aligned_haplotypes, segment_to_index_dict):
    index_to_segment_dict = reverse_dict(segment_to_index_dict)
    cyto_path = create_cytoband_path()
    vis_input = []
    for hap_idx, hap in enumerate(aligned_haplotypes):
        segment_list = indexed_segments_to_typed_segments(hap.mt_hap, index_to_segment_dict)
        c_entry = {'hap_id': hap.id,
                   'segment_list': hap.mt_hap,  # this remains un-altered
                   'chr': hap.chrom,
                   'name': hap.chrom,  # remove/reassign
                   'length': get_chr_length(segment_list),
                   'bands': label_cytoband(segment_list, cyto_path),
                   'orientation_contigs': get_orientation_contigs(segment_list),
                   'highlight': chr_is_highlighted(events, hap.id),
                   'sv_labels': []}
        vis_input.append(c_entry)
    assign_sv_labels(events, vis_input, index_to_segment_dict)

    # TODO: assign names to each entry (i.e. Chr1 -> Chr1a)

    return vis_input


def get_orientation_contigs(segment_list):
    """
    generate orientation contigs, continuous in chr_origin and orientation
    :param segment_list:
    :return:
    """
    ## find orientation contigs
    orientation_contigs = []
    current_orientation = True  # True for forward, False for backward
    current_idx = 0
    current_chr_origin = segment_list[0].chr_name
    previous_idx = 0
    for seg in segment_list:
        if len(seg) <= 2:
            continue
        seg_orientation = seg.direction()
        if seg_orientation == current_orientation and seg.chr_name == current_chr_origin:
            current_idx += len(seg)
        else:
            orientation_contigs.append({'start': previous_idx / 1e6, 'end': current_idx / 1e6,
                                        'length': (current_idx - previous_idx + 1) / 1e6, 'orientation': current_orientation,
                                        'origin': current_chr_origin.replace('Chr', '')})
            previous_idx = current_idx
            current_idx += len(seg)
            current_orientation = seg_orientation
            current_chr_origin = seg.chr_name
    # add the last contig
    orientation_contigs.append({'start': previous_idx / 1e6, 'end': current_idx / 1e6,
                                'length': (current_idx - previous_idx + 1) / 1e6, 'orientation': current_orientation,
                                'origin': current_chr_origin.replace('Chr', '')})

    return orientation_contigs



def chr_is_highlighted(input_events, hap_id):
    for event_info in input_events:
        for block in event_info[2]:
            event_path_id = int(block.split('.')[0])
            if hap_id == event_path_id:
                return True
    return False


def assign_sv_labels(input_events, all_vis_input, i_index_to_segment_dict):
    name_abbreviation = {'insertion': 'INS',
                         'deletion': 'DEL',
                         'inversion': 'INV',
                         'tandem_duplication': 'DUP',
                         'left_duplication_inversion': "DUPINV",
                         'right_duplication_inversion': 'DUPINV',
                         'balanced_translocation': 'TRANS'}

    def find_and_assign_single_label(path_id, indexed_seg, label_str):
        path_found = False
        for entry in all_vis_input:
            if path_id == entry['hap_id']:
                if indexed_seg == 'p-ter':
                    pos = 0
                elif indexed_seg == 'q-ter':
                    raise RuntimeError('q-ter found on left bp?')
                else:
                    idx_in_segment_list = entry['segment_list'].index(indexed_seg)
                    pos = 0
                    for seg in range(0, idx_in_segment_list + 1):
                        c_indexed_seg = entry['segment_list'][seg]
                        pos += len(i_index_to_segment_dict[int(c_indexed_seg[:-1])])
                entry['sv_labels'].append({'pos': pos / 1e6, 'label': label_str})
                path_found = True
                break
        if not path_found:
            raise RuntimeError('path not found')

    associated_event = []
    event_id = 1
    for event_idx, event_info in enumerate(input_events):
        event_type = event_info[1]
        if event_type in ['insertion', 'deletion', 'inversion', 'tandem_duplication',
                          'left_duplication_inversion', 'right_duplication_inversion']:
            left_segment = event_info[2][0].split('.')[3]
            event_name = name_abbreviation[event_type]
            c_path_idx = int(event_info[2][0].split('.')[0])
            find_and_assign_single_label(c_path_idx, left_segment, '[{}]{}'.format(event_id, event_name))
        elif event_type.startswith('balanced_translocation_unassociated'):
            left_segment1 = event_info[2][0].split('.')[3]
            left_segment2 = event_info[2][1].split('.')[3]
            event_name = name_abbreviation['balanced_translocation']
            c_path_idx1 = int(event_info[2][0].split('.')[0])
            c_path_idx2 = int(event_info[2][1].split('.')[0])
            find_and_assign_single_label(c_path_idx1, left_segment1, '[{}]{}'.format(event_id, event_name))
            find_and_assign_single_label(c_path_idx2, left_segment2, '[{}]{}'.format(event_id, event_name))
        elif event_type.startswith('balanced_translocation_associated'):
            if event_idx in associated_event:
                continue
            next_event_info = input_events[event_idx + 1]
            associated_event.append(event_idx + 1)
            event_name = name_abbreviation['balanced_translocation']
            if 'mt' in event_info[2][0]:
                left_segment1 = event_info[2][0].split('.')[3]
                c_path_idx1 = int(event_info[2][0].split('.')[0])
            elif 'mt' in event_info[2][1]:
                left_segment1 = event_info[2][1].split('.')[3]
                c_path_idx1 = int(event_info[2][1].split('.')[0])
            else:
                raise RuntimeError('mt string not found')
            if 'mt' in next_event_info[2][0]:
                left_segment2 = next_event_info[2][0].split('.')[3]
                c_path_idx2 = int(next_event_info[2][0].split('.')[0])
            elif 'mt' in next_event_info[2][1]:
                left_segment2 = next_event_info[2][1].split('.')[3]
                c_path_idx2 = int(next_event_info[2][1].split('.')[0])
            else:
                raise RuntimeError('mt string not found')
            find_and_assign_single_label(c_path_idx1, left_segment1, '[{}]{}'.format(event_id, event_name))
            find_and_assign_single_label(c_path_idx2, left_segment2, '[{}]{}'.format(event_id, event_name))
        event_id += 1


def indexed_segments_to_typed_segments(indexed_segment_list, index_to_segment_dict):
    typed_segment_list = []
    for seg in indexed_segment_list:
        direction = seg[-1]
        temp_seg = index_to_segment_dict[int(seg[:-1])].duplicate()
        if direction == '-':
            temp_seg.invert()
        typed_segment_list.append(temp_seg)
    return typed_segment_list


def create_cytoband_path(cyto_file='KT_Reports/Metadata/hg38_400_level_cytoband_updated.tsv'):
    segment_list = []
    with open(cyto_file) as fp_read:
        fp_read.readline()
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            if int(line[0]) == 23:
                chrom = 'ChrX'
            elif int(line[0]) == 24:
                chrom = 'ChrY'
            else:
                chrom = 'Chr' + line[0]
            start = int(line[1])
            end = int(line[2])
            band = line[7]
            stain = line[8]
            if '_' in chrom:
                # non-canonical labeling
                continue
            new_seg = Segment(chrom, start, end, band=band, stain=stain)
            segment_list.append(new_seg)
    return Path(Arm(segment_list, 'cytoband_arm'))


def label_cytoband(input_segment_list, input_cyto_path):
    c_path = Path(Arm(input_segment_list, 'to_label_cytoband'))
    c_path.generate_mutual_breakpoints(other_path=input_cyto_path, mutual=False)
    # now all segs in input_segment_list is a substring of a segment in the cyto_path
    cytoband_assigned = False
    for this_seg in c_path.linear_path.segments:
        for other_seg in input_cyto_path.linear_path.segments:
            if this_seg.chr_name == other_seg.chr_name and this_seg.start >= other_seg.start and this_seg.end <= other_seg.end:
                this_seg.band = other_seg.band
                this_seg.stain = other_seg.stain
                cytoband_assigned = True
                break
        if not cytoband_assigned:
            print(this_seg)
            raise RuntimeError('cytoband not found')

    ## gather band length and merge adjacent segments of the same band
    output_list = []
    p_end = 0
    c_band_len = 0
    p_band = c_path.linear_path.segments[0].band
    p_stain = c_path.linear_path.segments[0].stain
    for seg in c_path.linear_path.segments:
        c_band = seg.band
        if c_band == p_band:
            c_band_len += len(seg)
        else:
            p_band_dict = {'start': p_end,
                           'end': p_end + c_band_len,
                           'band': p_band,
                           'stain': p_stain}
            output_list.append(p_band_dict)
            p_end += c_band_len
            p_band = c_band
            p_stain = seg.stain
            c_band_len = len(seg)
    # last band
    last_band_dict = {'start': p_end,
                      'end': p_end + c_band_len,
                      'band': p_band,
                      'stain': p_stain}
    output_list.append(last_band_dict)

    ## scale band size by Mbp
    scaled_output_list = []
    for band in output_list:
        new_band = {'start': band['start'] / 1e6,
                    'end': band['end'] / 1e6,
                    'band': band['band'],
                    'stain': band['stain']}
        scaled_output_list.append(new_band)
    return scaled_output_list


def get_chr_length(input_segment_list):
    total_length = 0
    for seg in input_segment_list:
        total_length += len(seg)
    return total_length / 1e6


def max_chr_length(vis_input):
    max_length = -1
    for chrom in vis_input:
        c_length = chrom['length']
        if c_length > max_length:
            max_length = c_length
    return math.ceil(max_length)


#########################TESTS#########################

def test_artificial_chr_image():
    # Example data for multiple chromosomes
    chromosomes_data = [
        {
            'name': 'Chr1',
            'length': 125,
            'bands': [
                {'start': 0, 'end': 15, 'name': 'p15', 'stain': 'gneg'},
                {'start': 15, 'end': 25, 'name': 'p14', 'stain': 'gpos25'},
                {'start': 25, 'end': 35, 'name': 'p13', 'stain': 'gpos50'},
                {'start': 35, 'end': 45, 'name': 'p12', 'stain': 'gpos75'},
                {'start': 45, 'end': 55, 'name': 'p11.2', 'stain': 'acen'},
                {'start': 55, 'end': 65, 'name': 'q11.2', 'stain': 'acen'},
                {'start': 65, 'end': 75, 'name': 'q21', 'stain': 'gpos100'},
                {'start': 75, 'end': 85, 'name': 'q22', 'stain': 'gneg'},
                {'start': 85, 'end': 95, 'name': 'q23', 'stain': 'gvar'},
                {'start': 95, 'end': 105, 'name': 'q24', 'stain': 'gneg'},
                {'start': 105, 'end': 115, 'name': 'q25', 'stain': 'stalk'},
                {'start': 115, 'end': 125, 'name': 'q26', 'stain': 'gneg'}
            ],
            'origins': [
                {'start': 0, 'end': 95, 'name': '1'},
                {'start': 95, 'end': 125, 'name': '2'}
            ],
            'highlight': False,
            'sv_labels': [
                {'pos': 24, 'label': '[2]INS'}
            ]
        },
        {
            'name': 'Chr2',
            'length': 120,
            'bands': [
                {'start': 0, 'end': 10, 'name': 'p15', 'stain': 'gneg'},
                {'start': 10, 'end': 20, 'name': 'p14', 'stain': 'gpos25'},
                {'start': 20, 'end': 30, 'name': 'p13', 'stain': 'gpos50'},
                {'start': 30, 'end': 40, 'name': 'p12', 'stain': 'gpos75'},
                {'start': 40, 'end': 50, 'name': 'p11.2', 'stain': 'acen'},
                {'start': 50, 'end': 60, 'name': 'q11.2', 'stain': 'acen'},
                {'start': 60, 'end': 70, 'name': 'q21', 'stain': 'gpos100'},
                {'start': 70, 'end': 80, 'name': 'q22', 'stain': 'gneg'},
                {'start': 80, 'end': 90, 'name': 'q23', 'stain': 'gvar'},
                {'start': 90, 'end': 100, 'name': 'q24', 'stain': 'gneg'},
                {'start': 100, 'end': 110, 'name': 'q25', 'stain': 'stalk'},
                {'start': 110, 'end': 120, 'name': 'q26', 'stain': 'gneg'}
            ],
            'origins': [
                {'start': 0, 'end': 60, 'name': '3'},
                {'start': 60, 'end': 100, 'name': '5'},
                {'start': 100, 'end': 120, 'name': '1'}
            ],
            'highlight': True,
            'sv_labels': [
                {'pos': 58, 'label': '[2]DEL'},
                {'pos': 79, 'label': '[2]DUPINV'}
            ]
        },
        {
            'name': 'Chr2',
            'length': 120,
            'bands': [
                {'start': 0, 'end': 10, 'name': 'p15', 'stain': 'gneg'},
                {'start': 10, 'end': 20, 'name': 'p14', 'stain': 'gpos25'},
                {'start': 20, 'end': 30, 'name': 'p13', 'stain': 'gpos50'},
                {'start': 30, 'end': 40, 'name': 'p12', 'stain': 'gpos75'},
                {'start': 40, 'end': 50, 'name': 'p11.2', 'stain': 'acen'},
                {'start': 50, 'end': 60, 'name': 'q11.2', 'stain': 'acen'},
                {'start': 60, 'end': 70, 'name': 'q21', 'stain': 'gpos100'},
                {'start': 70, 'end': 80, 'name': 'q22', 'stain': 'gneg'},
                {'start': 80, 'end': 90, 'name': 'q23', 'stain': 'gvar'},
                {'start': 90, 'end': 100, 'name': 'q24', 'stain': 'gneg'},
                {'start': 100, 'end': 110, 'name': 'q25', 'stain': 'stalk'},
                {'start': 110, 'end': 120, 'name': 'q26', 'stain': 'gneg'}
            ],
            'origins': [
                {'start': 0, 'end': 60, 'name': '3'},
                {'start': 60, 'end': 100, 'name': '5'},
                {'start': 100, 'end': 120, 'name': '1'}
            ],
            'highlight': True,
            'sv_labels': [
                {'pos': 58, 'label': '[2]DEL'},
                {'pos': 79, 'label': '[2]T'}
            ]
        },
        {
            'name': 'Chr2',
            'length': 120,
            'bands': [
                {'start': 0, 'end': 10, 'name': 'p15', 'stain': 'gneg'},
                {'start': 10, 'end': 20, 'name': 'p14', 'stain': 'gpos25'},
                {'start': 20, 'end': 30, 'name': 'p13', 'stain': 'gpos50'},
                {'start': 30, 'end': 40, 'name': 'p12', 'stain': 'gpos75'},
                {'start': 40, 'end': 50, 'name': 'p11.2', 'stain': 'acen'},
                {'start': 50, 'end': 60, 'name': 'q11.2', 'stain': 'acen'},
                {'start': 60, 'end': 70, 'name': 'q21', 'stain': 'gpos100'},
                {'start': 70, 'end': 80, 'name': 'q22', 'stain': 'gneg'},
                {'start': 80, 'end': 90, 'name': 'q23', 'stain': 'gvar'},
                {'start': 90, 'end': 100, 'name': 'q24', 'stain': 'gneg'},
                {'start': 100, 'end': 110, 'name': 'q25', 'stain': 'stalk'},
                {'start': 110, 'end': 120, 'name': 'q26', 'stain': 'gneg'}
            ],
            'origins': [
                {'start': 0, 'end': 60, 'name': '3'},
                {'start': 60, 'end': 100, 'name': '5'},
                {'start': 100, 'end': 120, 'name': '1'}
            ],
            'highlight': True,
            'sv_labels': [
                {'pos': 58, 'label': '[2]DEL'},
                {'pos': 79, 'label': '[2]DUP'}
            ]
        }
    ]

    plt.rcParams['figure.dpi'] = 250
    n_chrom = len(chromosomes_data)
    n_row = n_chrom // 4 + 1
    fig, i_ax = plt.subplots(figsize=(5 * n_row, 2 * min(4, n_chrom)))

    for chrom_idx, i_chromosome_data in enumerate(chromosomes_data):
        row = chrom_idx // 4
        col = chrom_idx % 4
        plot_chromosome(i_ax, i_chromosome_data, y_offset=col * 3, x_offset=row * 28)

    # plt.show(bbox_inches='tight')
    plt.savefig('test_fig.png', bbox_inches='tight')
    rotate_image('test_fig.png', 'test_fig_rotated.png')


def make_image(vis_input, i_max_length, output_prefix, param_image_len_scale):
    plt.rcParams['figure.dpi'] = IMAGE_DPI

    if i_max_length <= MAX_CHR_LEN_IF_NO_SCALE:
        scaled_image_length = (i_max_length / 200) * 8 * param_image_len_scale
    else:
        scaled_image_length = (MAX_CHR_LEN_IF_NO_SCALE / 200) * 8 * param_image_len_scale

    n_chrom = len(vis_input)
    if n_chrom <= 4:
        image_width = 1.0 * 4
    else:
        image_width = 1.0 * 8
    fig, i_ax = plt.subplots(figsize=(scaled_image_length, image_width))

    ## Scale all Chr in the cluster if at least one Chr is too long to fit
    if i_max_length <= MAX_CHR_LEN_IF_NO_SCALE:
        chr_len_scaling = 1
    else:
        chr_len_scaling = MAX_CHR_LEN_IF_NO_SCALE / i_max_length
        for vis in vis_input:
            apply_scaling_to_vis(vis, chr_len_scaling)

    ## Merge SV-labels if they are too close
    for vis in vis_input:
        merge_sv_labels(vis, SV_LABEL_MIN_DISTANCE / chr_len_scaling)

    ## Limit chrom plot size
    i_ax.set_xlim(0, min(i_max_length, MAX_CHR_LEN_IF_NO_SCALE) + CHR_HEADER_X_OFFSET + 1.5)
    if n_chrom <= 4:
        ylim = 16
    elif n_chrom <= 8:
        ylim = 32
    else:
        ylim = 4 * n_chrom
    i_ax.set_ylim(0, ylim)
    i_ax.axis('off')  # turn on for debugging plotting locations

    ## generate CHR plotting location, depending on the number of chromosomes, fixed distance between two CHR
    Y_INIT_mapping = {1: 6,
                      2: 4,
                      3: 2,
                      4: 0,
                      5: 6,
                      6: 4,
                      7: 2,
                      8: 0}
    if len(vis_input) in Y_INIT_mapping:
        Y_INIT = Y_INIT_mapping[len(vis_input)]
    else:
        Y_INIT = 0
    Y_CONST = 4

    for chrom_idx, i_chromosome_data in enumerate(vis_input):
        row = chrom_idx // 4
        col = chrom_idx % 4
        plot_chromosome(i_ax, i_chromosome_data, Y_INIT + col * Y_CONST, row * 28, chr_len_scaling)

    plt.savefig(output_prefix + '.png', bbox_inches='tight')
    rotate_image(output_prefix + '.png', output_prefix + '_rotated.png')


def merge_sv_labels(vis_entry, min_distance):
    if len(vis_entry['sv_labels']) == 0:
        return
    new_sv_labels = [{'pos': vis_entry['sv_labels'][0]['pos'],
                      'label': [vis_entry['sv_labels'][0]['label']]}]
    last_pos = vis_entry['sv_labels'][0]['pos']
    last_pos_idx = 0
    for sv_label_idx, sv_label in enumerate(vis_entry['sv_labels'][1:]):
        c_pos = sv_label['pos']
        if c_pos - last_pos <= min_distance:
            new_sv_labels[last_pos_idx]['label'].append(sv_label['label'])
        else:
            new_sv_labels.append({'pos': sv_label['pos'],
                                  'label': [sv_label['label']]})
            last_pos = sv_label['pos']
            last_pos_idx = len(new_sv_labels) - 1
    vis_entry['sv_labels'] = new_sv_labels


def apply_scaling_to_vis(vis_entry, scaling_factor):
    """
    When adding any new entry to the vis parameters, update this
    :param vis_entry:
    :param scaling_factor:
    :return:
    """
    vis_entry['length'] = vis_entry['length'] * scaling_factor
    for band in vis_entry['bands']:
        band['start'] = band['start'] * scaling_factor
        band['end'] = band['end'] * scaling_factor
    # for origin in vis_entry['origins']:
    #     origin['start'] = origin['start'] * scaling_factor
    #     origin['end'] = origin['end'] * scaling_factor
    for sv_label in vis_entry['sv_labels']:
        sv_label['pos'] = sv_label['pos'] * scaling_factor
    for orientation in vis_entry['orientation_contigs']:
        orientation['start'] = orientation['start'] * scaling_factor
        orientation['end'] = orientation['end'] * scaling_factor
        orientation['length'] = orientation['length'] * scaling_factor


if __name__ == '__main__':
    omkar_file_path = '/Users/zhaoyangjia/PyCharm_Repos/KarComparator/real_case_data/dremsek_OMKar_output_paths/39.txt'
    mt_indexed_lists, mt_path_chrs, segment_to_index_dict, segment_size_dict = read_OMKar_to_indexed_list(omkar_file_path)
    mt_path_chrs = [info.split(': ')[-1] for info in mt_path_chrs]
    wt_path_dict = generate_wt_from_OMKar_output(segment_to_index_dict)
    wt_indexed_lists = populate_wt_indexed_lists(mt_path_chrs, wt_path_dict)
    events, aligned_haplotypes = interpret_haplotypes(mt_indexed_lists, wt_indexed_lists, mt_path_chrs, segment_size_dict)

    c_vis_input = generate_visualizer_input(events, aligned_haplotypes, segment_to_index_dict)
    vis_input_used = [c_vis_input[1], c_vis_input[3], c_vis_input[5], c_vis_input[5]]
    print(max_chr_length(vis_input_used))
    make_image([c_vis_input[1], c_vis_input[3], c_vis_input[5], c_vis_input[5]], max_chr_length(vis_input_used), 'test_new', IMG_LENGTH_SCALE_VERTICAL_SPLIT)
    # create_cytoband_path()

    # event<0>,type<balanced_translocation_unassociated>,blocks<['44.1.mt(44+).46+.47+', '45.0.wt(44+).p-ter.45+']>

    