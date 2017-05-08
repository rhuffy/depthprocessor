#!/usr/bin/env python

# Raymond Huffman 2016
#
# Inputs: (Depth file produced by Samtools Depth), (GTF file produced by CLASS)
# Depth file MUST be sorted

import argparse

import numpy as np

parser = argparse.ArgumentParser(description='This script compiles depth \
                                              statistics from a Depth and \
                                              a GTF file.')
parser.add_argument('depth', help='Depth file path.\
                                   Expected format: [chr]\t[position]\t[depth]')
parser.add_argument('gtf', help='GTF file path')
parser.add_argument('output', help='Output file path')

args = parser.parse_args()

# DEPTH = 'depth.txt'
# GTF = 'class.gtf'
# OUT = 'class.dgtf'


def exon_stack(gtf_file):
    master_gene_list = list()
    current_gene = str()
    i = -1
    for line in gtf_file:
        if line[:1] != '#':
            read_line = line.split()  # [2] - exon/transcript, [3] - start, [4] - end, [9] - ID
            if read_line[9][1:-2] != current_gene:  # start of new gene, must be a transcript, checking is unnecessary
                i += 1  # iterate index in list of lists
                master_gene_list.append(list())  # each item is a bracket list, containing pos, identity, gene_id, transcript_id
                master_gene_list[i].append((int(read_line[3]), 0, read_line[9][1:-2], read_line[11][1:-2]))  # add soft open
                master_gene_list[i].append((int(read_line[4]), 3, read_line[9][1:-2], read_line[11][1:-2]))  # add soft close
                current_gene = read_line[9][1:-2]
            else:
                if read_line[2] == 'exon':
                    master_gene_list[i].append((int(read_line[3]), 2, read_line[9][1:-2], read_line[11][1:-2]))  # add hard open
                    master_gene_list[i].append((int(read_line[4]), 1, read_line[9][1:-2], read_line[11][1:-2]))  # add hard close
                else:
                    master_gene_list[i].append((int(read_line[3]), 0, read_line[9][1:-2], read_line[11][1:-2]))  # add soft open
                    master_gene_list[i].append((int(read_line[4]), 3, read_line[9][1:-2], read_line[11][1:-2]))  # add soft close
    for gene_list in master_gene_list:
        gene_list.sort()
     #   outstr = str()
     #   for item in gene_list:  # visualization
     #       if item[1] == 0:
     #           outstr += '('
     #       if item[1] == 3:
     #           outstr += ')'
     #       if item[1] == 2:
     #           outstr += '['
     #       if item[1] == 1:
     #           outstr += ']'
    return master_gene_list


def old_analyze_gene(gene_list, depth, output):  # gene_list is a sorted list of all exon boundaries of all transcripts
    bracket_level = 0  # represents the depth of nested brackets at current pos
    possible_event = False
    first_pos = True
    possible_event_start = 0
    possible_event_end = 0
    relative_intron_pos = 1
    for bracket in gene_list:  # even brackets are open, odd are close
        if bracket[1] == 2:  # open
            first_pos = False
            bracket_level += 1
            if possible_event:
                possible_event_end = int(bracket[0])
                # run a function that computes the depth of reads at end of possible retention event
                run_get_depth(possible_event_start, possible_event_end, depth, bracket[2], bracket[3], relative_intron_pos, output)
                possible_event = False
                relative_intron_pos += 1
        elif bracket[1] == 1:  # close
            bracket_level -= 1
        if not possible_event and bracket_level == 0 and not first_pos:  # a zero here indicates start of possible intron retention event
            possible_event = True
            possible_event_start = int(bracket[0])  # mark the start of this possible event


def analyze_gene(gene_list, depth, output):
    exon = False
    valid = False
    possible_event_start = 0
    possible_event_end = 0
    relative_intron_pos = 1
    for bracket in gene_list:
        if exon:
            if bracket[1]%2 == 1:  # if odd --> if close bracket
                possible_event_start = bracket[0]  # mark potential end of exon, start of intron event
                exon = False
            # else:  # if even --> if open bracket, do nothing

        else:  # in intron
            if bracket[1]%2 == 0:  # if even --> if open bracket
                if valid:
                    possible_event_end = bracket[0]
                    run_get_depth(possible_event_start, possible_event_end, depth, bracket[2], bracket[3], relative_intron_pos, output)
                    relative_intron_pos += 1
                exon = True
                valid = True
            else:  # if odd --> if close bracket
                possible_event_start = bracket[0]


def gene_id_to_chrom(gene_id):
    chrom = gene_id[:5]
    if chrom[4] == '.':
        chrom = chrom[:-1]
    return chrom


def chrom_to_num(chrom):
    num = chrom[3:]
    if num == 'X':
        num = 23
    elif num == 'Y':
        num = 24
    elif num == 'M':
        num = 25
    else:
        try:
            num = int(num)
        except TypeError:
            num = 0
    return num


def read_depth(file):
    chr_depth_list = list()
    for i in range(26):  # creates lists for all chromosomes, ignore 0, X is 24, Y is 25, M is 26
        chr_depth_list.append(np.zeros(250000000, dtype=np.int))
    for line in file:
        read = line.split()
        chrom_int = (chrom_to_num(read[0]))  # converts chromosome 'chr1' into int '1'
        pos_int = int(read[1])
        chr_depth_list[chrom_int][pos_int] = int(read[2])
    #    if pos_int % 1000000 == 0:
    #        print(line)
    return chr_depth_list


def get_depth(chrom, start, end, depth):
    sum_to_return = 0
    length = end - start
    relevant_length = length
    for i in range(start, end+1):
        sum_to_return += depth[chrom][i]
        if depth[chrom][i] == 0:
            relevant_length -= 1
    return sum_to_return, length, relevant_length


def run_get_depth(start, end, depth, gene_id, trans_id, relative_intron_pos, output):
    chrom = gene_id_to_chrom(gene_id)
    depth_info = get_depth(chrom_to_num(chrom), start, end, depth)
    output_line = chrom + '\t' + str(start) + '\t' + str(end) + '\t' + \
                str(depth_info[0]) + '\t' + str(depth_info[1]) + '\t' + \
                str(depth_info[2]) + '\t'
    if float(depth_info[1]) != 0:
        output_line += str(float(depth_info[0]) / float(depth_info[1])) + '\t'
    else:
        output_line += '0\t'
    if float(depth_info[2]) != 0:
        output_line += str(float(depth_info[0]) / float(depth_info[2])) + '\t'
    else:
        output_line += '0\t'
    output_line += str(relative_intron_pos) + '\t' + gene_id + '\t' + trans_id + '\n'
    output.write(output_line)


def run(depth, gtf, output):
    print('start python')
    depth_file = open(depth, 'r')
    gtf_file = open(gtf, 'r')
    out_file = open(output, 'w')
    print('start exon stack')
    master_gene_list = exon_stack(gtf_file)
    print('start read depth')
    depth_arr = read_depth(depth_file)
    print('finish read depth')
    out_file.write('chrom\tsrt_pos\tend_pos\tdepth\tlength\talign_len\tavg\
                    \talign_avg\tintron_num\tgene_id\ttrans_id\n')
    for gene in master_gene_list:
        analyze_gene(gene, depth_arr, out_file)
    #    print('gene done')
    print('python done')
    depth_file.close()
    gtf_file.close()
    out_file.close()


run(args.depth, args.gtf, args.output)
