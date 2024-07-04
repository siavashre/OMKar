from Structures import Arm, Segment, Path


def read_forbidden_regions(forbidden_region_file) -> Arm:
    segment_list = []
    with open(forbidden_region_file) as fp_read:
        fp_read.readline()  # skip index line
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            new_segment = Segment(str(line[0]), int(line[1]), int(line[2]), str(line[3]))
            segment_list.append(new_segment)
    return Arm(segment_list, 'forbidden_regions')


def output_forbidden_regions_from_arm(input_arm: Arm, output_file_path):
    with open(output_file_path, 'w') as fp_write:
        fp_write.write('Chr\tStartPos\tEndPos\tType\n')
        for segment_itr in input_arm.segments:
            output_str = "{}\t{}\t{}\t{}\n".format(segment_itr.chr_name,
                                                   str(segment_itr.start),
                                                   str(segment_itr.end),
                                                   segment_itr.segment_type)
            fp_write.write(output_str)


def label_path_with_forbidden_regions(input_path_list: [Path], forbidden_region_file):
    forbidden_regions_path = Path(read_forbidden_regions(forbidden_region_file), 'forbidden_regions', 'forbidden_regions')
    for path_itr in input_path_list:
        # breakup into disjoint segments
        path_itr.generate_mutual_breakpoints(forbidden_regions_path, mutual=False)

        # label forbidden regions
        for path_segment_itr in path_itr.linear_path.segments:
            labeled = False
            for forbidden_region_segment_itr in forbidden_regions_path.linear_path.segments:
                if path_segment_itr.same_segment_ignore_dir(forbidden_region_segment_itr):
                    path_segment_itr.segment_type = forbidden_region_segment_itr.segment_type
                    labeled = True
            if not labeled:
                path_segment_itr.segment_type = 'arm_region'


def get_chr_length_from_forbidden_file(input_chr_name, forbidden_region_file='Metadata/acrocentric_telo_cen.bed'):
    with open(forbidden_region_file) as fp_read:
        fp_read.readline()  # skip index line
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            chrom = line[0]
            start_pos = int(line[1])
            end_pos = int(line[2])
            seg_type = line[3]
            if chrom.lower() == input_chr_name.lower():
                if seg_type == 'telomere2':
                    return end_pos
    return None



def test():
    print(read_forbidden_regions('../Metadata/merged_forbidden_regions_unique.bed'))


if __name__ == "__main__":
    test()
