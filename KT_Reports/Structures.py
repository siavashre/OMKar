import numpy as np


class Segment:
    chr_name: str
    start: int
    end: int
    segment_type: str
    kt_index: str
    ordinal: int
    band: str
    stain: str

    def __init__(self, chr_name: str, start: int, end: int, segment_type=None, kt_index=None, band='', stain=''):
        self.chr_name = chr_name
        self.start = start
        self.end = end
        self.segment_type = segment_type
        self.kt_index = kt_index
        self.ordinal = -1
        self.band = band
        self.stain = stain

    def __len__(self):
        return abs(self.end - self.start) + 1

    def __lt__(self, other):
        type_order = ["telomere1",
                      "centromere",
                      "telomere2",
                      "acrocentric",
                      'acrocentric-telomere1',
                      'acrocentric-centromere',
                      "arm_region",
                      "hardmask",
                      "superdup"]

        def get_chr_order(chromosome_name):
            chr_extracted = chromosome_name.replace('Chr', '')
            if chr_extracted == 'X':
                return 23
            elif chr_extracted == 'Y':
                return 24
            else:
                return int(chr_extracted)

        if get_chr_order(self.chr_name) < get_chr_order(other.chr_name):
            return True
        elif get_chr_order(self.chr_name) > get_chr_order(other.chr_name):
            return False
        elif get_chr_order(self.chr_name) == get_chr_order(other.chr_name):
            if max(self.start, self.end) != max(other.start, other.end):
                return max(self.start, self.end) < max(other.start, other.end)
            else:
                self_type_index = type_order.index(self.segment_type)
                other_type_index = type_order.index(other.segment_type)
                return self_type_index < other_type_index

    def __eq__(self, other):
        if isinstance(other, Segment):
            return (self.chr_name, self.start, self.end) == (other.chr_name, other.start, other.end)
        return False
        # return (self.chr_name, self.start, self.end) == (other.chr_name, other.start, other.end)

    def __hash__(self):
        return hash((self.chr_name, self.start, self.end))

    def __str__(self):
        additional_info = ""
        if self.segment_type is not None:
            additional_info += ", " + self.segment_type
        if self.kt_index is not None:
            additional_info += ", " + self.kt_index
        if self.ordinal != -1:
            additional_info += ", " + str(self.ordinal)
        return_str = "({}, {}, {}{})".format(self.chr_name, self.start, self.end, additional_info)
        return return_str

    def thousand_delimited(self):
        return "({}-{}-{})".format(self.chr_name, "{:,}".format(self.start), "{:,}".format(self.end))

    def same_segment_ignore_dir(self, other):
        if self.chr_name != other.chr_name:
            return False
        if self.start != other.start and self.start != other.end:
            return False
        if self.end != other.start and self.end != other.end:
            return False
        return True

    def direction(self):
        """
        :return: 1 for +, 0 for -
        """
        return self.start <= self.end

    def duplicate(self):
        return Segment(self.chr_name, self.start, self.end, self.segment_type, self.kt_index)

    def left_delete(self, bp_to_delete):
        """
        :param bp_to_delete: number of bp deleting
        :return: None
        """
        if self.direction():
            self.start = self.start + bp_to_delete
        else:
            self.start = self.start - bp_to_delete

    def right_delete(self, bp_to_delete):
        """
        :param bp_to_delete: number of bp deleting
        :return: None
        """
        if self.direction():
            self.end = self.end - bp_to_delete
        else:
            self.end = self.end + bp_to_delete

    def invert(self, inplace=True):
        if inplace:
            temp_start = self.start
            self.start = self.end
            self.end = temp_start
            if self.kt_index is not None:
                kt_direction = self.kt_index[-1]
                if kt_direction == '+':
                    self.kt_index = self.kt_index[:-1] + '-'
                else:
                    self.kt_index = self.kt_index[:-1] + '+'
        else:
            new_segment = self.duplicate()
            new_segment.start = self.end
            new_segment.end = self.start
            if new_segment.kt_index is not None:
                kt_direction = new_segment.kt_index[-1]
                if kt_direction == '+':
                    new_segment.kt_index = new_segment.kt_index[:-1] + '-'
                else:
                    new_segment.kt_index = new_segment.kt_index[:-1] + '+'
            return new_segment

    def segment_intersection(self, other_segment):
        """
        Report intersection regardless of directions
        :param other_segment:
        :return:
        """
        duplicate_self = self.duplicate()
        duplicate_other = other_segment.duplicate()

        # put both into + directions
        if not duplicate_self.direction():
            duplicate_self.invert()
        if not duplicate_other.direction():
            duplicate_other.invert()

        if duplicate_self.chr_name != duplicate_other.chr_name:
            return False
        if duplicate_self.start <= duplicate_other.end and duplicate_self.end >= duplicate_other.start:
            return True
        else:
            return False

    def bp_in_interior(self, bp_chromosome, bp_index, bp_type):
        """
        For KarComparator
        Used to see if an index is within the segment, potentially requiring a breakpoint
        :param bp_chromosome:
        :param bp_index:
        :param bp_type: start or end
        :return:
        """
        if self.chr_name != bp_chromosome:
            return False
        if self.direction():
            if bp_type == "start" and bp_index == self.start:
                return False
            elif bp_type == "end" and bp_index == self.end:
                return False

            if bp_index in range(self.start, self.end + 1):
                return True
            else:
                return False
        else:
            if bp_type == "end" and bp_index == self.start:
                return False
            elif bp_type == "start" and bp_index == self.end:
                return False

            if bp_index in range(self.end, self.start + 1):
                return True
            else:
                return False

    def is_continuous(self, other):
        if self.chr_name != other.chr_name:
            return False

        if self.direction():
            if other.start == self.end + 1:
                return True
        else:
            if other.start == self.end - 1:
                return True

        return False

    def assign_cn_bin(self, cn_bins):
        cn = []
        for bin_idx, cn_bin in enumerate(cn_bins):
            cn_value = 0.0
            if self.chr_name.lower() == cn_bin['chrom'].lower():
                # compute inclusion length
                seg_start, seg_end = sorted((self.start, self.end))
                overlap_start = max(seg_start, cn_bin['start'])
                overlap_end = min(seg_end, cn_bin['end'])
                if overlap_start < overlap_end:
                    overlap_length = overlap_end - overlap_start + 1
                    bin_length = cn_bin['end'] - cn_bin['start'] + 1
                    cn_value = overlap_length / bin_length
            cn.append(cn_value)
        return np.array(cn)


class Arm:
    segments: [Segment]
    deleted: bool
    arm_type: str

    def __init__(self, segments: [Segment], arm_type: str):
        self.segments = segments
        self.deleted = False
        self.arm_type = arm_type

    def __len__(self):
        current_sum = 0
        for segment in self.segments:
            current_sum += len(segment)
        return current_sum

    def __contains__(self, item):
        return any(item == e for e in self.segments)

    def __str__(self):
        return_str = ''
        for segment in self.segments:
            return_str += str(segment)
        return return_str

    def duplicate(self):
        new_segments = []
        for segment in self.segments:
            new_segments.append(segment.duplicate())
        return Arm(new_segments, self.arm_type)

    def delete_segments_by_index(self, segment_indices):
        self.segments = [segment for index, segment in enumerate(self.segments) if index not in segment_indices]

    def duplicate_segments_by_index(self, segment_indices):
        # segments come in order, so insert before the first segment
        index_of_insertion = segment_indices[0]
        new_segments = []
        for index in segment_indices:
            new_segments.append(self.segments[index].duplicate())

        self.segments[index_of_insertion:index_of_insertion] = new_segments

    def invert_segments_by_index(self, segment_indices):
        # segments come in order
        index_of_insertion = segment_indices[0]
        new_segments = []
        for index in reversed(segment_indices):
            new_segment = self.segments[index].duplicate()
            new_segment.invert()
            new_segments.append(new_segment)
        self.delete_segments_by_index(segment_indices)
        self.segments[index_of_insertion:index_of_insertion] = new_segments

    def arm_intersection(self, other_arm):
        """
        If two arms have any segments that intersect, regardless of directions
        :param other_arm:
        :return:
        """
        for segment1 in self.segments:
            for segment2 in other_arm.segments:
                if segment1.segment_intersection(segment2):
                    return True
        return False

    def report_arm_intersection(self, other_arm):
        intersecting_segments = []
        for segment1 in self.segments:
            for segment2 in other_arm.segments:
                if segment1.segment_intersection(segment2):
                    intersecting_segments.append(segment2)
        return_str = ''
        for segment in intersecting_segments:
            return_str += segment.annotated_number()
        return return_str

    def gather_boundary_points(self):
        """
        For KarComparator
        :return: a list of all the boundary points for each Segment (two for each)
        """
        return_list = []
        for segment in self.segments:
            if segment.direction():
                return_list.append(tuple([segment.chr_name, segment.start, 'start']))
                return_list.append(tuple([segment.chr_name, segment.end, 'end']))
            else:
                return_list.append(tuple([segment.chr_name, segment.start, 'end']))
                return_list.append(tuple([segment.chr_name, segment.end, 'start']))
        return return_list

    def introduce_breakpoint(self, bp_chromosome, bp_index, bp_type):
        """
        For KarComparator
        Search through the arm and generate the breakpoint, if within an interior of a Segment
        :param bp_chromosome:
        :param bp_index:
        :param bp_type: start or end bp
        :return:
        """
        current_segment_index = 0
        while current_segment_index < len(self.segments):
            current_segment = self.segments[current_segment_index]
            if current_segment.bp_in_interior(bp_chromosome, bp_index, bp_type):
                insertion_index = self.get_segment_index(current_segment)
                if current_segment.direction():
                    if bp_type == "start":
                        left_segment = Segment(current_segment.chr_name, current_segment.start, bp_index - 1,
                                               current_segment.segment_type, current_segment.kt_index)
                        right_segment = Segment(current_segment.chr_name, bp_index, current_segment.end,
                                                current_segment.segment_type, current_segment.kt_index)
                    elif bp_type == "end":
                        left_segment = Segment(current_segment.chr_name, current_segment.start, bp_index,
                                               current_segment.segment_type, current_segment.kt_index)
                        right_segment = Segment(current_segment.chr_name, bp_index + 1, current_segment.end,
                                                current_segment.segment_type, current_segment.kt_index)
                    else:
                        raise ValueError('bp_type must be start OR end')
                else:
                    if bp_type == "start":
                        left_segment = Segment(current_segment.chr_name, current_segment.start, bp_index,
                                               current_segment.segment_type, current_segment.kt_index)
                        right_segment = Segment(current_segment.chr_name, bp_index - 1, current_segment.end,
                                                current_segment.segment_type, current_segment.kt_index)
                    elif bp_type == "end":
                        left_segment = Segment(current_segment.chr_name, current_segment.start, bp_index + 1,
                                               current_segment.segment_type, current_segment.kt_index)
                        right_segment = Segment(current_segment.chr_name, bp_index, current_segment.end,
                                                current_segment.segment_type, current_segment.kt_index)
                    else:
                        raise ValueError('bp_type must be start OR end')

                self.segments.pop(insertion_index)
                self.segments.insert(insertion_index, right_segment)
                self.segments.insert(insertion_index, left_segment)
                current_segment_index += 1  # since we added one more segment in-place

            current_segment_index += 1

    def get_segment_index(self, input_segment):
        """
        For KarComparator
        find the index of segment in list, matching the object using __is__ (not __eq__)
        :param input_segment: Only use segment that is in the Arm
        :return:
        """
        for segment_index in range(0, len(self.segments)):
            current_segment = self.segments[segment_index]
            if current_segment is input_segment:
                return segment_index
        raise RuntimeError('segment not found in Arm')

    def merge_breakpoints(self):
        """
        merge all continuous breakpoints
        :return:
        """
        current_segment_index = 0
        while current_segment_index < len(self.segments) - 1:
            current_segment = self.segments[current_segment_index]
            next_segment = self.segments[current_segment_index + 1]
            if current_segment.chr_name == next_segment.chr_name:
                if current_segment.segment_type == next_segment.segment_type:
                    if current_segment.direction() and current_segment.end + 1 == next_segment.start:
                        new_segment = Segment(current_segment.chr_name, current_segment.start,
                                              next_segment.end, current_segment.segment_type)
                        self.segments.pop(current_segment_index)
                        self.segments.pop(current_segment_index)
                        self.segments.insert(current_segment_index, new_segment)
                        continue
                    elif (not current_segment.direction()) and current_segment.end - 1 == next_segment.start:
                        new_segment = Segment(current_segment.chr_name, current_segment.start,
                                              next_segment.end, current_segment.segment_type)
                        self.segments.pop(current_segment_index)
                        self.segments.pop(current_segment_index)
                        self.segments.insert(current_segment_index, new_segment)
                        continue

            current_segment_index += 1


class Chromosome:
    name: str
    p_arm: Arm
    q_arm: Arm
    centromere: Arm
    t1_len: int
    t2_len: int
    deleted: bool

    def __init__(self, name: str, p_arm: Arm, q_arm: Arm, t1_len: int, t2_len: int, centromere: Arm, deleted=False):
        self.name = name
        self.p_arm = p_arm
        self.q_arm = q_arm
        self.t1_len = t1_len
        self.t2_len = t2_len
        self.centromere = centromere
        self.deleted = deleted

    def __len__(self):
        if self.deleted:
            return 0
        else:
            return self.p_arm_len() + self.q_arm_len()

    def __str__(self):
        if self.deleted:
            return '{}: deleted'.format(self.name)
        return_str = '{}: t1-{} t2-{}\n\tp-arm: {}\n\tq-arm: {}\n\tCEN: {}' \
            .format(self.name, self.t1_len, self.t2_len, str(self.p_arm), str(self.q_arm), str(self.centromere))
        return return_str

    def __iter__(self):
        class ChromosomeIterator:
            def __init__(self, chromosome: Chromosome):
                self.chromosome = chromosome
                self.arms = [chromosome.p_arm, chromosome.centromere, chromosome.q_arm]
                self.current_arm_index = 0
                self.current_segment_index = 0

            def __next__(self):
                if self.current_arm_index < len(self.arms):
                    current_arm = self.arms[self.current_arm_index]
                    if current_arm.deleted:
                        self.current_arm_index += 1
                        return next(self)
                    elif self.current_segment_index < len(current_arm.segments):
                        segment = current_arm.segments[self.current_segment_index]
                        self.current_segment_index += 1
                        return segment
                    else:
                        self.current_arm_index += 1
                        self.current_segment_index = 0
                        return next(self)
                else:
                    raise StopIteration

        return ChromosomeIterator(self)

    def p_arm_len(self):
        if self.p_arm.deleted:
            return 0
        else:
            return len(self.p_arm)

    def q_arm_len(self):
        if self.q_arm.deleted:
            return 0
        else:
            return len(self.q_arm)

    def duplicate(self):
        return Chromosome(self.name, self.p_arm.duplicate(), self.q_arm.duplicate(),
                          self.t1_len, self.t2_len, self.centromere.duplicate())

    def contains_junction(self, segment1, segment2):
        """
        See if the edge/junction between the two segments are in this arm
        Note 1+ -> 2+ is different from 1- -> 2+
        :param segment1:
        :param segment2:
        :return: True if junction is present, false OW
        """
        # TODO: implement
        pass


class Genome:
    full_KT: {str: [Chromosome]}  # has as many slots as there are chromosome type, i.e. 24 for a male, 23 for a female
    motherboard: Arm  # using the Arm object to use generate breakpoint method
    centromere_segments = [Segment]
    histories = []

    def __init__(self, full_KT, motherboard_segments, centromere_segments, histories):
        self.full_KT = full_KT
        self.motherboard = Arm(motherboard_segments, 'motherboard')
        self.centromere_segments = centromere_segments
        self.histories = histories

    def sort_histories(self):
        # sort histories by Chr_from's names, but maintain temporal order
        def sort_key(history_entry):
            event_type = history_entry[0]
            chr_from = history_entry[1]
            info = chr_from.replace('Chr', '')
            numeric = int(info[:-1])
            letter = info[-1]
            if event_type == 'balanced reciprocal translocation':
                # keep balanced translocation together
                return 0, 'a'
            else:
                return numeric, letter

        self.histories = sorted(self.histories, key=sort_key)

    def translate_histories_from_indexing(self):
        def translate_segments_from_indices(indexed_segments):
            actual_segments = []
            for indexed_seg in indexed_segments:
                if indexed_seg.startswith('CEN'):
                    idx = int(indexed_seg[-2])
                    direction = "+"
                    segment = self.centromere_segments[idx - 1].duplicate()
                else:
                    idx = int(indexed_seg[:-1])
                    direction = indexed_seg[-1]
                    segment = self.motherboard.segments[idx - 1].duplicate()
                if direction == '-':
                    segment.invert()
                actual_segments.append(segment)
            return actual_segments

        new_histories = []
        for hist in self.histories:
            event_segments = translate_segments_from_indices(hist[3])
            new_histories.append((hist[0], hist[1], hist[2], event_segments))
        self.histories = new_histories

    def duplicate(self):
        new_full_KT = {}
        for key in self.full_KT:
            new_chr_list = []
            for chr_itr in self.full_KT[key]:
                new_chr_list.append(chr_itr.duplicate())
            new_full_KT[key] = new_chr_list

        new_centromere_segments = []
        for chr_itr in self.centromere_segments:
            new_centromere_segments.append(chr_itr.duplicate())

        return Genome(new_full_KT,
                      self.motherboard.duplicate(),
                      new_centromere_segments,
                      self.histories)

    def __str__(self):
        return_str = ''
        for chromosome in self:
            return_str += str(chromosome) + '\n'
        return return_str

    def __iter__(self):
        def custom_sort_chr(key):
            chr_part = key[3:]  # Extract the part after "Chr"
            if chr_part.isdigit():
                return int(chr_part)
            elif chr_part == "X":
                return int(23)  # Put ChrX at the end
            elif chr_part == "Y":
                return int(24)  # Put ChrY after ChrX
            return key

        class GenomeIterator:
            def __init__(self, genome: Genome):
                self.genome = genome
                self.KT_slots = genome.full_KT
                self.KT_slot_keys = sorted(genome.full_KT.keys(), key=custom_sort_chr)
                self.current_slot_index = 0
                self.current_chromosome_index = 0

            def __next__(self):
                if self.current_slot_index < len(self.KT_slots):
                    current_slot = self.KT_slots[self.KT_slot_keys[self.current_slot_index]]
                    if self.current_chromosome_index < len(current_slot):
                        chromosome = current_slot[self.current_chromosome_index]
                        self.current_chromosome_index += 1
                        return chromosome
                    else:
                        self.current_slot_index += 1
                        self.current_chromosome_index = 0
                        return next(self)
                else:
                    raise StopIteration

        return GenomeIterator(self)

    def get_chromosome_list(self):
        chr_list = []
        for chromosome in self:
            chr_list.append(chromosome)
        return chr_list

    def segment_indexing(self):
        segment_dict = {}
        current_index = 1
        for segment_itr in self.motherboard.segments:
            segment_dict[segment_itr] = str(current_index)
            current_index += 1
        for centromere_itr in self.centromere_segments:
            centromere_name = centromere_itr.chr_name.replace('Chr', 'CEN')
            segment_dict[centromere_itr] = centromere_name
        return segment_dict

    def motherboard_tostring(self):
        segment_dict = self.segment_indexing()
        sorted_segments = sorted(segment_dict)
        return_str = 'index\torigin\tstart\tend\n'
        for segment_itr in sorted_segments:
            return_str += '{}\t{}\t{}\t{}\n'.format(segment_dict[segment_itr], segment_itr.chr_name,
                                                    segment_itr.start, segment_itr.end)
        return return_str

    def KT_tostring(self):
        segment_dict = self.segment_indexing()
        return_str = 'chromosome\tKT\ttelo1_len\ttelo2_len\n'

        for chr_itr in self:
            if chr_itr.deleted:
                return_str += '{}\tdeleted\t0\t0\n'.format(chr_itr.name)
                continue

            tostring_segment_list = []
            for segment_itr in chr_itr:
                if segment_itr.direction():
                    tostring_segment_list.append(segment_dict[segment_itr] + '+')
                else:
                    new_segment_itr = segment_itr.duplicate()
                    new_segment_itr.invert()
                    tostring_segment_list.append(segment_dict[new_segment_itr] + '-')

            return_str += '{}\t{}\t{}\t{}\n'.format(chr_itr.name, ','.join(tostring_segment_list),
                                                    str(chr_itr.t1_len), str(chr_itr.t2_len))
        return return_str

    # def need_breakpoint(self, event_arm: Arm, breakpoint_index: int):
    #     """
    #     split segment such that the breakpoint_index is guaranteed to be the end index of a Segment
    #     :param event_arm: Arm which the event happens on, and the breakpoint_index point at
    #     :param breakpoint_index: the position of break on the current Arm
    #         (left_event_index - 1) OR (right_event_index)
    #     :return: None
    #     """
    #     # TODO: understand what is going on with all the unused
    #     if breakpoint_index == -1:
    #         # this happens when the break point is at the very beginning of the event_arm, no breaking required
    #         return 0
    #
    #     segment_to_break = Segment('temp', -1, -1)
    #     left_delete_len = -1
    #     right_delete_len = -1
    #
    #     current_bp_index = -1  # corrects 0-index off-shift
    #
    #     # locate the Segment to create breakpoint
    #     for segment in event_arm.segments:
    #         current_bp_index += len(segment)
    #         if current_bp_index == breakpoint_index:
    #             # breakpoint exists
    #             return 0
    #         elif current_bp_index > breakpoint_index:
    #             # document the breakpoint location on the current segment
    #             segment_to_break = segment.duplicate()
    #             previous_bp_index = current_bp_index - len(segment)
    #             left_delete_len = breakpoint_index - previous_bp_index
    #             right_delete_len = current_bp_index - breakpoint_index
    #             break
    #         else:
    #             # breakpoint location not yet met
    #             continue
    #     return 1

    def generate_breakpoint(self, event_arm: Arm, breakpoint_index: int):
        """
        split segment such that the breakpoint_index is garenteed to be the end index of a Segment
        :param event_arm: Arm which the event happens on, and the breakpoint_index point at
        :param breakpoint_index: the position of break on the current Arm
            (left_event_index - 1) OR (right_event_index)
        :return: None
        """
        if breakpoint_index == -1:
            # this happens when the break point is at the very beginning of the event_arm, no breaking required
            return

        segment_to_break = Segment('temp', -1, -1)
        left_delete_len = -1
        right_delete_len = -1

        current_bp_index = -1  # corrects 0-index off-shift

        # locate the Segment to create breakpoint
        for segment in event_arm.segments:
            current_bp_index += len(segment)
            if current_bp_index == breakpoint_index:
                # breakpoint exists
                return
            elif current_bp_index > breakpoint_index:
                # document the breakpoint location on the current segment
                segment_to_break = segment.duplicate()
                previous_bp_index = current_bp_index - len(segment)
                left_delete_len = breakpoint_index - previous_bp_index
                right_delete_len = current_bp_index - breakpoint_index
                break
            else:
                # breakpoint location not yet met
                continue

        # create breakpoint on all identical Segments in the genome
        def break_segment(current_arm: Arm):
            left_segment = segment_to_break.duplicate()
            right_segment = segment_to_break.duplicate()
            left_segment.right_delete(right_delete_len)
            right_segment.left_delete(left_delete_len)

            same_direction_match = \
                [index for index, value in enumerate(current_arm.segments) if value == segment_to_break]
            for segment_index_itr in reversed(same_direction_match):
                current_arm.segments.pop(segment_index_itr)
                current_arm.segments.insert(segment_index_itr, left_segment.duplicate())
                current_arm.segments.insert(segment_index_itr + 1, right_segment.duplicate())

            all_match = \
                [index for index, value in enumerate(current_arm.segments)
                 if segment_to_break.same_segment_ignore_dir(value)]
            same_direction_match = \
                [index for index, value in enumerate(current_arm.segments) if value == segment_to_break]
            reversed_direction_match = [element for element in all_match if element not in same_direction_match]
            for segment_index_itr in reversed(reversed_direction_match):
                current_arm.segments.pop(segment_index_itr)
                new_right_segment = right_segment.duplicate()
                new_left_segment = left_segment.duplicate()
                new_right_segment.invert()
                new_left_segment.invert()
                current_arm.segments.insert(segment_index_itr, new_right_segment)
                current_arm.segments.insert(segment_index_itr + 1, new_left_segment)

        break_segment(self.motherboard)
        for slot in self.full_KT:
            for chromosome in self.full_KT[slot]:
                break_segment(chromosome.p_arm)
                break_segment(chromosome.q_arm)

    def output_KT(self, output_file):
        with open(output_file, 'w') as fp_write:
            fp_write.write(self.motherboard_tostring())
            fp_write.write('---\n')
            fp_write.write(self.KT_tostring())


class Path:
    linear_path: Arm
    path_chr: str
    path_name: str

    def __init__(self, linear_path: Arm, path_name=None, path_chr=None):
        self.linear_path = linear_path
        self.path_chr = path_chr
        self.path_name = path_name

    def __str__(self):
        return str("chr_bin: {}, path_name: {}, segments: {}".format(self.path_chr, self.path_name, self.linear_path))

    def concise_str(self):
        segment_str = ""
        for segment in self.linear_path.segments:
            segment_str += str(segment)
        return str("path_name: {}, segments: {}".format(self.path_name, segment_str))

    def reverse(self):
        new_segments = []
        for segment in reversed(self.linear_path.segments):
            new_segment = segment.duplicate()
            new_segment.invert()
            new_segments.append(new_segment)
        self.linear_path.segments = new_segments

    def get_path_notes(self):
        segment_origin_str = ""
        for segment in self.linear_path.segments:
            segment_origin_str += segment.chr_name + " "
        return str("chr_bin: {}, path_name: {}, segment_chr: {}".format(self.path_chr, self.path_name,
                                                                        segment_origin_str))

    def generate_mutual_breakpoints(self, other_path=None, mutual=True):
        """
        make sure all segments within the one/two path/s have mutual breakpoints
        :param other_path: if None, then breaking within itself
        :param mutual: whether to generate breakpoints on the other_path
        :return:
        """
        path1_breakpoints = self.linear_path.gather_boundary_points()

        if other_path is not None:
            path2_breakpoints = other_path.linear_path.gather_boundary_points()
            for breakpoint_itr in path2_breakpoints:
                self.linear_path.introduce_breakpoint(*breakpoint_itr)

            if mutual:
                for breakpoint_itr in path1_breakpoints:
                    other_path.linear_path.introduce_breakpoint(*breakpoint_itr)
        else:
            for breakpoint_itr in path1_breakpoints:
                self.linear_path.introduce_breakpoint(*breakpoint_itr)

    def is_disjoint(self):
        tmp_self = self.duplicate()
        before_breaking_len = len(tmp_self.linear_path.segments)
        tmp_self.generate_mutual_breakpoints()
        after_breaking_len = len(tmp_self.linear_path.segments)
        if before_breaking_len == after_breaking_len:
            return True
        else:
            return False

    def duplicate(self):
        new_arm = self.linear_path.duplicate()
        return Path(new_arm, self.path_chr, self.path_name)

    def get_origins(self):
        origins = set()
        for segment in self.linear_path.segments:
            segment_chr = segment.chr_name.split('-')[0]
            if segment_chr not in origins:
                origins.add(segment_chr)
        return origins

    def tostring_path_by_index(self, segment_to_index):
        output = "{}\t{}\t".format(self.path_name, self.path_chr)
        segment_indices = []
        for segment in self.linear_path.segments:
            current_segment = segment.duplicate()
            if current_segment in segment_to_index:
                segment_indices.append(segment_to_index[current_segment] + '+')
            else:
                current_segment.invert()
                # if this is no-found error, mistake was made during the dict disjoint process
                # print(segment_to_index)
                # print(current_segment)
                segment_indices.append(segment_to_index[current_segment] + '-')

        return output + ','.join(segment_indices)

    def nonforbidden_len(self):
        length = 0
        forbidden_comparison_region_types = ['acrocentric', 'telomere1', 'telomere2', 'acrocentric-telomere1', 'acrocentric-centromere']
        for segment in self.linear_path.segments:
            if segment.segment_type not in forbidden_comparison_region_types:
                length += len(segment)
        return length


def segment_indices_to_segments(segment_index_list, segment_dict):
    """
    :param segment_index_list: a list of segment indices with direction (eg. [1+, 2-, 12+])
    :param segment_dict: key is int, value is Segment
    :return: a list of Segments in the same order
    """
    segment_list = []
    for segment_index_element in segment_index_list:
        segment_index = int(segment_index_element[:-1])
        segment_direction = segment_index_element[-1]

        new_segment = segment_dict[segment_index].duplicate()
        if segment_direction == "-":
            new_segment.invert()
            new_segment.kt_index = str(segment_index) + '-'
        else:
            new_segment.kt_index = str(segment_index) + '+'
        segment_list.append(new_segment)
    return segment_list


def flip_dict(input_dict):
    """
    input dict requires one-to-one correspondence
    :param input_dict:
    :return:
    """
    output_dict = {}
    for key, value in input_dict.items():
        output_dict[value] = key
    return output_dict


if __name__ == "__main__":
    pass
