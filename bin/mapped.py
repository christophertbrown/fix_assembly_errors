#!/usr/bin/python2.7

"""
script for sorting and filtering a sam file
"""

import sys
import os
from itertools import cycle
from subprocess import Popen, PIPE
import argparse

def sam2fastq(line):
    """
    print fastq from sam
    """
    fastq = []
    fastq.append('@%s' % line[0])
    fastq.append(line[9])
    fastq.append('+%s' % line[0])
    fastq.append(line[10])
    return fastq

def count_mismatches(read):
    """
    look for NM:i:<N> flag to determine number of mismatches
    """
    if read is False:
        return False
    mm = [int(i.split(':')[2]) for i in read[11:] if i.startswith('NM:i:')]
    if len(mm) > 0:
        return sum(mm)
    else:
        return False

def check_mismatches(read, pair, mismatches, mm_option, req_map):
    """
    - check to see if the read maps with <= threshold number of mismatches
    - mm_option = 'one' or 'both' depending on whether or not one or both reads
       in a pair need to pass the mismatch threshold
    - pair can be False if read does not have a pair
    - make sure alignment score is not 0, which would indicate that the read was not aligned to the reference
    """
    # if read is not paired, make sure it is mapped and that mm <= thresh
    if pair is False:
        mm = count_mismatches(read)
        if mm is False:
            return False
        # if no threshold is supplied, return True
        if mismatches is False:
            return True
        # passes threshold?
        if mm <= mismatches:
            return True
    # paired reads
    r_mm = count_mismatches(read)
    p_mm = count_mismatches(pair)
    # if neither read is mapped, return False
    if r_mm is False and p_mm is False:
        return False
    # if no threshold, return True
    if mismatches is False:
        return True
    # if req_map is True, both reads have to map
    if req_map is True:
        if r_mm is False or p_mm is False:
            return False
    ## if option is 'one,' only one read has to pass threshold
    if mm_option == 'one':
        if (r_mm is not False and r_mm <= mismatches) or (p_mm is not False and p_mm <= mismatches):
            return True
    ## if option is 'both,' both reads have to pass threshold 
    if mm_option == 'both':
        ## if one read in pair does not map to the scaffold, 
        ## make sure the other read passes threshold
        if r_mm is False:
            if p_mm <= mismatches:
                return True
        elif p_mm is False:
            if r_mm <= mismatches:
                return True
        elif (r_mm is not False and r_mm <= mismatches) and (p_mm is not False and p_mm <= mismatches):
            return True 
    return False

def reads_from_mapping(mapping, contigs, mismatches, mm_option, req_map):
    c = cycle([1, 2])
    for line in mapping:
        line = line.strip().split()
        if line[0].startswith('@'): # get the sam header
            if line[0].startswith('@SQ'):
                contig = line[1].split('SN:', 1)[1]
                if contigs is not False:
                    if contig in contigs:
                        yield [0, line]
                else:
                    yield [0, line]
            else:
                yield [0, line]
            continue
        if int(line[1]) <= 20: # is this from a single read?
            if contigs is False or line[2] in contigs:
                if check_mismatches(line, False, mismatches, mm_option, req_map) is True:
                    yield [1, sam2fastq(line)]
                    yield [10, line]
        else:
            n = c.next()
            if n == 2:
                if contigs is False:
                    if prev[2] != '*' or line[2] != '*':
                        if check_mismatches(line, prev, mismatches, mm_option, req_map) is True:
                            yield [2, sam2fastq(prev)]
                            yield [20, prev]
                            yield [2, sam2fastq(line)]
                            yield [20, line]
                else:
                    if prev[2] in contigs or line[2] in contigs:
                        if check_mismatches(line, prev, mismatches, mm_option, req_map) is True:
                            yield [2, sam2fastq(prev)]
                            yield [20, prev]
                            yield [2, sam2fastq(line)]
                            yield [20, line]
            prev = line

def get_reads(sam, \
        contigs = False, mismatches = False, mm_option = False, sort_sam = True, req_map = False):
    """
    get mapped reads (and their pairs) from an unsorted sam file
    """
    if sort_sam is True:
        mapping = '%s.sorted.sam' % (sam.rsplit('.', 1)[0])
        if sam != '-':
            if os.path.exists(mapping) is False:
                os.system("cat %s | \
                    sort -k1 --buffer-size=5000M -o %s\
                    " % (sam, mapping)) 
        else:
            mapping = 'stdin-sam.sorted.sam'
            p = Popen("sort -k1 --buffer-size=5000M -o %s" \
                    % (mapping), stdin = sys.stdin, shell = True) 
            p.communicate()
        mapping = open(mapping)
    else:
        if sam == '-':
            mapping = sys.stdin
        else:
            mapping = open(sam)
    for read in reads_from_mapping(mapping, contigs, mismatches, mm_option, req_map):
        yield read

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# filter sam file based on mismatches')
    parser.add_argument(\
            '-s', required = True, help = 'path to sam file')
    parser.add_argument(\
            '-m', required = True, help = 'maximum number of mismatches (or False to include all mapped)')
    parser.add_argument(\
            '-p', required = True, help = 'require that "one" or "both" reads in pair have <= m mismatches')
    parser.add_argument(\
            '--require-mapping', action = 'store_true', help = 'require both reads are mapped')
    parser.add_argument(\
            '-o', default = False, help = 'name for new sam file')
    parser.add_argument(\
            '-f', default = False, help = 'filter based on scaffold name (single name or list from stdin (-))')
    parser.add_argument(\
            '-r', action = 'store_true', help = 'print mapped paired reads to stdout and single reads to stderr')
    parser.add_argument(\
            '-n', action = 'store_false', help = 'do not sort sam file')
    args = vars(parser.parse_args())
    # is -o or -r specified? If not, script won't output anything
    if args['o'] is False and args['r'] is False:
        print '# specify -o and/or -r'
        exit()
    sam, contigs, mismatches, mm_option, new_sam = args['s'], args['f'], args['m'], args['p'], args['o']
    print_reads, sort_sam, req_map = args['r'], args['n'], args['require_mapping']
    if mismatches == 'False' or mismatches == 'FALSE' or mismatches == 'false':
        mismatches = False
    if mismatches is False: # just to make sure mm_option is not specified if mismatches is not specified
        mm_option = False
    else:
        mismatches = int(mismatches)
        if mm_option != 'one' and mm_option != 'both':
            print >> sys.stderr, '# specify one or both for mismatch option'
            print >> sys.stderr, '# i.e. should the mismatches threshold apply to one or both reads in a pair?'
            exit()
    if contigs == '-':
        contigs = [i.strip() for i in sys.stdin]
    if new_sam is not False:
        new_sam = open(new_sam, 'w')
    for type, read in get_reads(sam, contigs, mismatches, mm_option, sort_sam, req_map):
        if type == 1:
            if print_reads is True:
                print >> sys.stderr, '\n'.join(read)
        elif type == 2:
            if print_reads is True:
                print '\n'.join(read)
        elif new_sam is not False and (type == 0 or type == 10 or type == 20):
            print >> new_sam, '\t'.join(read)
    if new_sam is not False:
        new_sam.close()
