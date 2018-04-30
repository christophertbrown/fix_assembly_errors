#!/usr/bin/env python3

"""
Script for fixing scaffolding errors.

 1. Errors are identified based on regions of zero coverage by stringently
    mapped reads (stringent mapping for both reads in a pair).
 2. Errors are re-assembled using reads stringently mapped to
    error-containing region (stringent mapping for one read in a pair).
 3. Errors are replaced with the re-assembled sequence, if successful.
 4. Errors that could not be fixed are either left alone or replaced with
    Ns (option).
 5. Scaffolds are split if insert coverage is zero.

Chris Brown
ctb@berkeley.edu

to-do
* add length of errors to report
* check bowtie insert length parameters
* save log of bowtie commands
"""

# python modules
import re
import os
import sys
import json
import shutil
import argparse
import itertools
import subprocess
from glob import glob as glob

# ctb modules
import ctbBio.mapped as map_tool
from ctbBio.fastq_split import split as fastq_split
from ctbBio.nr_fasta import de_rep as fix_fasta
from ctbBio.fasta import iterate_fasta as parse_fasta
from ctbBio.rc import reverse_complement as rc
from ctbRA.assemble import velvet as velvet

def run_bowtie(assembly, sam, pr, pr_split, sr, threads, multiple, bt_dir):
    """
    map reads to assembly with bowite2
    """
    if os.path.exists(sam) is False: # if sam file does not exist, run bowtie2
        db_check = '%s/%s.1.bt2' % (bt_dir, assembly)
        sr_command = pr_command = matches_command = ''
        if os.path.exists(db_check) is False: # check for db file
            p = subprocess.Popen('bowtie2-build -q %s %s/%s' \
                    % (assembly, bt_dir, assembly.rsplit('/', 1)[-1]), shell = True)
            p.communicate()
        if sr is not False:
            sr_command = '-U %s' % (sr)
        if pr_split is not False:
            r1, r2 = pr_split
            pr_command = '-1 %s -2 %s' % (r1, r2)
        if pr is not False:
            base = '%s/%s' % (bt_dir, pr.rsplit('/', 1)[-1].rsplit('.', 1)[0])
            r1 = '%s.R1.fastq' % (base)
            r2 = '%s.R2.fastq' % (base)
            pr_split = [r1, r2]
            if os.path.exists(r1) is False:
                fastq_split(open(pr), base)
            pr_command = '-1 %s -2 %s' % (r1, r2)
        if multiple is True: # should reads be allowed to map to more than one position or contig?
            matches_command = '-a '
        p = subprocess.Popen('\
                bowtie2 --very-fast --reorder --quiet %s -p %s -x %s/%s %s %s | shrinksam -k %s' \
                % (matches_command, threads, bt_dir, assembly.rsplit('/', 1)[-1], pr_command, sr_command, sam)\
                , shell = True)
        p.communicate()
    return sam, pr_split

def map_reads(assembly, scaffolds, pr, threads, multiple, bt_dir = 'bt2', pr_split = False, sr = False):
    """
    map reads to assembly with bowite2
    # calls run_bowtie for running the mapping
    """
    if '/' in assembly:
        bt_dir = '%s/%s' % (assembly.rsplit('/', 1)[0], bt_dir)
        assembly_dir, assembly_name = assembly.rsplit('/', 1)
    else:
        assembly_dir, assembly_name = './', assembly
    os.system('mkdir -p %s' % (bt_dir))
    sam = '%s/%s.bt.sam' % (bt_dir, assembly_name.rsplit('.', 1)[0])
    # check to make sure assembly has sequences
    if os.path.getsize(assembly) == 0:
        return False, pr_split
    return run_bowtie(assembly, sam, pr, pr_split, sr, threads, multiple, bt_dir) # run bowtie, return sam file

def check_cov(c, cov_thresh):
    """
    check:
    1) read coverage
    2) passing of overlap test
    3) insert coverage
    c = [coverage, overlap = True or False, insert coverage = True or False]
    """
    if c[0] >= cov_thresh and c[1] is True:
        return True
    return False

def errors_from_cov(s2c, cov_thresh):
    """
    identify errors based on coverage threshold
    """
    errors = {} # errors[scaffold] = [position of errors]
    for scaf, cov in list(s2c.items()):
        errors[scaf] = []
        prev = False
        for pos, c in enumerate(cov):
            if check_cov(c, cov_thresh) is False:
                if prev is True:
                    continue
                prev = True
                errors[scaf].append(pos)
            else:
                prev = False
            cov[pos] = c
    return errors

def mm_positions_from_md(sam, read_length):
    """
    get read positions from MD flag in sam file
    """
    positions = False
    md = [i.rsplit(':', 1)[1] for i in sam if i.startswith('MD:Z:')][0]
    if '^' in md or '+' in md: # do not count reads with indels 
        return True
    if md == str(read_length - 1):
        return positions
    last = len([i for i in re.split('[0-9]', md) if len(i) > 0])
    for i, pos in enumerate([i for i in re.split('[A-Z]', md) if i.isdigit()]):
        if i >= last:
            continue
        pos = int(pos) + 1
        if positions is False:
            positions = [pos]
        else:
            positions.append(positions[-1] + pos)
    if positions is False:
        return True
    return set(positions)

def check_mm(sam, window, read_length):
    """
    make sure mismatches are not in window at beginning or end of read
    if mismatches are in the beginning or end of the read, return False
    """
    mm = map_tool.count_mismatches(sam)
    if mm is False:
        return True
    if mm == 0:
        return False
    mm_positions = mm_positions_from_md(sam, read_length)
    if mm_positions is False:
        return False
    elif mm_positions is True:
        return True
    for pos in mm_positions:
        if pos <= window or pos >= (read_length - window):
            return True
    return False

def add_coverage(scaffold, scaffolds, overlap, s2c, sam, window):
    """
    add coverage to scaffolds for read region defined by overlap
    cov_info = [coverage, overlap_coverage, insert coverage]
    """
    read_length = overlap[1] - overlap[0] + 1
    mm = check_mm(sam, window, read_length)
    for i in range(overlap[0] - 1, overlap[1]):
        if scaffold not in s2c:
            return s2c
        try:
            s2c[scaffold][i][0] += 1
            if i < overlap[1] - window and mm is False:
                s2c[scaffold][i][1] = True
        except IndexError:
            break
    if sam[6] != '=': # paired read mapped to different scaffold (don't calc. insert cov.)
        return s2c
    insert = [min(overlap), min(overlap) + int(sam[8])]
    for i in range(min(insert) - 1, max(insert)):
        try:
            s2c[scaffold][i][2] = True
        except IndexError:
            break
    return s2c

def id_errors(pairs, header, assembly, scaffolds, \
        cov_thresh, mismatches, allow_orphan = False, \
        allow_orphan_ends = False, window = 3, \
        orphan_window = 1000, save_mapping = False):
    """
    identify regions with zero coverage by stringently mapped paired reads
    * pairs[read] = [bit, mate, mappping[scaffold] = [map, map2, ...], fastq]
    *  map = [overlap, mismatches, sam_info]
    *  sam_info = all sam lines except read + quality
    """
    if save_mapping is True:
        # open sam file for writing
        out = open('%s.both.sam' % (assembly.rsplit('.', 1)[0]), 'w')
        for line in header:
            print(line, file=out)
    # s2c[scaffold] = [[coverage per position, # reads connecting base to the next base], [p2, pn2]]
    s2c = {id: [[0, False, False] for i in range(0, info[1])] for id, info in list(scaffolds.items())}
    # filter reads
    for read in list(pairs.values()):
        bit, mate, maps, fastq = read
        if mate not in pairs:
            continue
        mate = pairs[mate]
        for scaffold, mappings in list(maps.items()):
            for mapping in mappings:
                overlap, mm, sam_info = mapping
                if scaffold not in mate[2]:
                    continue
                for mate_map in mate[2][scaffold]:
                    mate_overlap, mate_mm, mate_sam = mate_map
                    if mm is False and mate_mm is False: # make sure at least one pair mapped to scaffold
                        continue
                    if (allow_orphan is False and allow_orphan_ends is False) \
                            and (mm is False or mate_mm is False): # make sure both pairs mapped
                        continue
                    elif allow_orphan_ends is True and (mm is False or mate_mm is False):
                        pos = [sam_info[0][3], mate_sam[0][3]]
                        if min(pos) > orphan_window and max(pos) < scaffolds[scaffold][1] - orphan_window:
                            continue
                    if mm > mismatches or mate_mm > mismatches: # both pairs must pass mismatches
                        continue
                    sam = mapping[2][0] + [read[3][1]] + [read[3][3]] + mapping[2][1]
                    if save_mapping is True:
                        print('\t'.join(sam), file=out)
                    s2c = add_coverage(scaffold, scaffolds, overlap, s2c, sam, window)
    errors = errors_from_cov(s2c, cov_thresh)
    return s2c, errors

def define_windows(scaffolds, s2errors, window, combine_windows):
    """
    determine what window to use for collecting reads
    if windows overlap, combine them
    """
    # define windows
    s2windows = {} # s2windows[scaffold][error] = [window]
    for scaffold, errors in list(s2errors.items()):
        if len(errors) == 0:
            continue
        s2windows[scaffold] = {}
        for error in errors:
            start = (error - int(window/2))
            if start < 0:
                start = 0
                stop = window
                if stop > scaffolds[scaffold][1]:
                    stop = scaffolds[scaffold][1]
            else:
                stop = (error + int(window/2))
                if stop > scaffolds[scaffold][1]:
                    stop = scaffolds[scaffold][1]
                    start = stop - window
                    if start < 0:
                        start = 0
            s2windows[scaffold][error] = [start, stop]
    if combine_windows is False:
        return s2windows
    # combine overlapping windows
    updated = {}
    for scaffold in s2windows:
        updated[scaffold] = {}
        errors = sorted(list(s2windows[scaffold].items()), key = lambda x: x[0])
        if len(errors) > 1:
            i = 1
            for error in errors[1:]:
                prev = errors[i - 1]
                if prev[1][1] >= error[1][0]:
                    pos = '_'.join([str(prev[0]), str(error[0])])
                    start, stop = prev[1][0], error[1][1]
                    merged = (pos, [start, stop])
                    del errors[i - 1]
                    errors[i - 1] = merged
                else:
                    i += 1
        for error in errors:
            updated[scaffold][error[0]] = error[1]
    return updated

def get_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def check_overlap(overlap, window):
    """
    check to see if overlap is within window
    """
    if overlap is False:
        return False
    if get_overlap(overlap, window) > 0:
        return True
    return False

def map2window(scaffold, s2windows, s2errors, overlap, m_overlap): 
    """
    determine if reads maps within window
    return errors that reads map to, or False
    """
    matches = []
    if scaffold not in s2windows:
        return False
    errors = s2windows[scaffold]
    for error, window in list(errors.items()):
        if check_overlap(overlap, window) is True or check_overlap(m_overlap, window) is True:
            matches.append(error)
    if len(matches) == 0:
        return False
    return matches

def collect_reads(pairs, assembly, scaffolds, mismatches, \
        prefix, s2errors, s2windows, max_pairs = 1000):
    """
    collect reads mapped to scaffold within window
    - only one read in pair has to pass mismatch critera
    * pairs[read] = [bit, mate, mappping[scaffold] = map, fastq]
       map = [overlap, mismatches, sam_info]
       sam_info = all sam lines except read + quality
    * reads[scaffold][error] = [[sequences-fastq], filename]
    """
    # make files for reads
    reads = {} # reads[scaffold][error] = pe
    for scaffold, errors in list(s2windows.items()):
        if len(errors) == 0:
            continue
        reads[scaffold] = {}
        for error in errors:
            dir = '%s/s_%s_e_%s' % (prefix, scaffold, error)
            dir = dir[0:100] # make sure file name is not too long
            os.system('mkdir -p %s' % (dir))
            reads[scaffold][error] = [{}, '%s/reads.pe.fastq' % (dir)]
    # get reads
    for id, read in list(pairs.items()):
        id, num = id.rsplit('_', 1)
        bit, mate_id, maps, fastq = read
        if mate_id not in pairs:
            continue
        mate = pairs[mate_id]
        for scaffold, mappings in list(maps.items()):
            for mapping in mappings:
                overlap, mm, sam_info = mapping
                m_fastq, m_num = mate[3], mate_id.rsplit('_', 1)[1]
                if scaffold not in mate[2]:
                    m_info = [[False, False, False]]
                else:
                    m_info = [map for map in mate[2][scaffold]]
                for map in m_info:
                    m_overlap, m_mm, m_sam_info = map
                    if mismatches is not False and mm > mismatches and m_mm > mismatches:
                        # skip if both reads have too many mismatches
                        continue
                    if (mm is False and m_mm > mismatches) or (m_mm is False and mm > mismatches):
                        continue
                    if mm is False and m_mm is False:
                        continue
                    else:
                        errors = map2window(scaffold, s2windows, s2errors, overlap, m_overlap)
                    if errors is False:
                        continue
                    for error in errors:
                        fq = [None, None]
                        fq[int(num) - 1] = '\n'.join(fastq)
                        fq[int(m_num) - 1] = '\n'.join(m_fastq)
                        reads[scaffold][error][0][id] = '\n'.join(fq)
    # save reads to files
    for scaffold in reads:
        for error in reads[scaffold]:
            seqs, name = reads[scaffold][error]
            out = open(name, 'w')
            for seq in list(seqs.values())[0:max_pairs]:
                print(seq, file=out)
            out.close()
            reads[scaffold][error] = out.name
    return reads

def break_by_coverage(assembled, s2c, cov_thresh, ignore_ends = False):
    """
    break scaffolds where coverage is below cov_thresh
    """
    fragments = [] # fragments = [[[len, header], [sequence]], ...]
    if ignore_ends is False:
        for id, seq in list(assembled.items()):
            sequence = []
            for i, base in enumerate(s2c[id]):
                if check_cov(base, cov_thresh) is True:
                    if sequence == []:
                        start = i
                    sequence.append(seq[2][i].upper())
                else:
                    if sequence != []:
                        fragments.append([[len(sequence), '>%s_e:%s' % (id, start)], ''.join(sequence)])
                    sequence = []
            if sequence != []:
                fragments.append([[len(sequence), '>%s_e:%s' % (id, start)], ''.join(sequence)])
    else:
        for id, seq in list(assembled.items()):
            sequence = []
            for i, base in enumerate(s2c[id]):
                if check_cov(base, cov_thresh) is True or (i < 100 or i > (seq[1] - 100)):
                    if sequence == []:
                        start = i
                    sequence.append(seq[2][i].upper())
                else:
                    if sequence != []:
                        fragments.append([[len(sequence), '>%s_e:%s' % (id, start)], ''.join(sequence)])
                    sequence = []
            if sequence != []:
                fragments.append([[len(sequence), '>%s_e:%s' % (id, start)], ''.join(sequence)])
    fragments.sort(key = lambda x: x[0][0], reverse = True)
    return fragments

def re_assemble_velvet(pr, prefix, scaffold, error, re_assembled_fasta, scaffolding, min_contig):
    """
    re-assemble reads using velvet
    """
    out = '%s/s_%s_e_%s' % (prefix, scaffold, error)
    out = out[0:100] # make sure file name is not too long
    velvet(paired = [pr], out = out, scaffolding = scaffolding, \
        silent = True, min_contig = min_contig, kmer_min = 21, kmer_max = 71, kmer_increase = 10)
    re_assembled_seqs = {}
    for fasta in glob('%s/*.fasta' % (out)):
        for seq in parse_fasta(fasta):
            if seq[0] == []:
                continue
            seq[0] = '>%s_e:%s_%s' % (scaffold, error, seq[0].split('>')[1])
            re_assembled_seqs[seq[0].split('>')[1]] = [seq[0], len(seq[1]), seq[1]]
            print('\n'.join(seq), file=re_assembled_fasta)
    return re_assembled_seqs

def re_assemble_errors(reads, prefix, pr, pr_split, threads, cov_thresh, \
        mismatches, multiple, save_mapping = False, scaffolding = True, \
        assembler = 'velvet', min_contig = '200'):
    """
    - re-assemble mis-assembled regions
    - return re-assembled fragments (sorted longest to shortest)
    * assembler == 'velvet'
    """
    re_assembled_errors = {}
    re_assembled_seqs = {}
    re_assembled_fasta = open('%s/re_assembled_errors.fa' % (prefix), 'w')
    id2error = {}
    for scaffold, errors in list(reads.items()):
        re_assembled_errors[scaffold] = {}
        for error, e_pr in list(errors.items()):
            if assembler == 'velvet':
                ra = re_assemble_velvet(e_pr, prefix, scaffold, error, \
                                        re_assembled_fasta, scaffolding, \
                                        min_contig)
            else:
                print('# specify valid assembler', file=sys.stderr)
                exit()
            if len(ra) == 0:
                re_assembled_errors[scaffold][error] = False
                continue
            re_assembled_errors[scaffold][error] = []
            for i in list(ra.keys()):
                id2error[i] = [scaffold, error]
            for seq, info in list(ra.items()):
                re_assembled_seqs[seq] = info
    re_assembled_fasta.close
    mapping, pr_split = map_reads(re_assembled_fasta.name, \
            re_assembled_seqs, pr, threads, multiple, \
            pr_split = pr_split)
    stringent_pairs, header = \
            parse_mapping_stringent(mapping, re_assembled_seqs, \
            mismatches)
    s2c, s2errors = id_errors(stringent_pairs, header, \
            re_assembled_fasta.name, re_assembled_seqs, \
            cov_thresh, mismatches, allow_orphan = True, save_mapping = save_mapping)
    fragments = break_by_coverage(re_assembled_seqs, s2c, cov_thresh, ignore_ends = False)
    f_out = open('%s/re_assembled_fragments.fa' % (prefix), 'w')
    for f in fragments:
        print('\n'.join([f[0][1], f[1]]), file=f_out)
        scaffold, error = id2error[f[0][1].split('>')[1].rsplit('_', 1)[0]]
        re_assembled_errors[scaffold][error].append(f)
    f_out.close()
    return re_assembled_errors

def patch_start(orig, cov, error, n_error, patches, k, cov_thresh, buffer):
    """
    use re-assembled fragments to extend scaffold from start
    """
    length = 0
    trimmed = False
    start = []
    for i, c in enumerate(cov):
        if length < k + buffer:
            if check_cov(c, cov_thresh) is False:
                start = []
            else:
                start.append(orig[i])
                length = len(start)
        else:
            start = ''.join(start).upper()
            if n_error is False:
                trimmed = orig[(i - k - buffer):]
            else:
                trimmed = orig[(i - k - buffer):n_error]
            break
    if trimmed is False: # no region of original sequence passed cov. threshold
        return False
    attempts = []
    for patch in patches:
        patch = patch[1].upper()
        i = patch.find(start)
        if i == -1:
            patch = rc(['', patch])[1]
            i = patch.find(start)
        if i != -1:
            patch = patch[0:i]
            if len(patch) > 0:
                attempts.append([len(patch), patch])
    if len(attempts) > 0:
        best = sorted(attempts, key = lambda x: x[0], reverse = True)[0][1]
        return [[error, 'e', best], [error, 'o', trimmed]]
    else: # what if you could not extend the scaffold?
        return [[error, 'o', trimmed]]

def patch_end(orig, cov, error, n_error, patches, k, origA, origB, start, stop, buffer):
    """
    use re-assembled fragments to extend scaffold from end
    """
    attempts = []
    for patch in patches:
        patch = patch[1].upper()
        i = patch.find(start)
        if i == -1:
            patch = rc(['', patch])[1]
            i = patch.find(start)
        if i != -1:
            patch = patch[(i + k + buffer):]
            attempts.append([len(patch), patch])
    if len(attempts) > 0:
        best = sorted(attempts, key = lambda x: x[0], reverse = True)[0][1]
        return [error, 'e', best]
    else: # what if you could not extend the scaffold?
        return False

def patch_middle(orig, cov, error, n_error, patches, k, origA, origB, origM, start, stop, buffer):
    """
    use re-assembled fragments to patch mis-assembly in middle of scaffold
    """
    attemptsA = []
    for patch in patches:
        patch = patch[1].upper()
        p_start = patch.find(start)
        if p_start == -1:
            patch = rc(['', patch])[1]
            p_start = patch.find(start) 
        if p_start == -1:
            continue
        p = patch[(p_start + k + buffer):]
        attemptsA.append([len(p), p])
        p_stop = patch.find(stop) 
        if p_stop == -1:
            patch = rc(['', patch])[1]
            p_stop = patch.find(stop) 
        if p_stop == -1:
            continue
        p = patch[(p_start + k + buffer):(p_stop)]
        if len(p) == 0:
            continue
        return [0, p]
    attemptsB = []
    for patch in patches: # extend origB if possible
        patch = patch[1].upper()
        p_start = patch.find(start) 
        if p_start == -1:
            patch = rc(['', patch])[1]
            p_start = patch.find(start) 
        if p_start != -1:
            continue
        p_stop = patch.find(stop) 
        if p_stop == -1:
            patch = rc(['', patch])[1]
            p_stop = patch.find(stop) 
        if p_stop == -1:
            continue
        p = patch[:p_stop]
        attemptsB.append([len(p), p])
    # what if only the start kmer could be found?
    if len(attemptsA) > 0:
        bestA = sorted(attemptsA, key = lambda x: x[0], reverse = True)[0][1]
        if len(bestA) > 20:
            if len(attemptsB) > 0:
                bestB = sorted(attemptsB, key = lambda x: x[0], reverse = True)[0][1]
                if len(bestB) > 20:
                    return [12, [bestA, bestB]]
            return [1, bestA]
    # what if only the stop kmer could be found?
    if len(attemptsB) > 0:
        bestB =  sorted(attemptsB, key = lambda x: x[0], reverse = True)[0][1]
        if len(bestB) > 20:
            return [2, bestB]
    return False # could not find start and stop in any fragment
 
def find_start_stop_error(orig, cov, error, n_error, patches, k, cov_thresh, buffer):
    """
    find regions of length k flanking region with error
    """
    start = orig[(error - k - buffer):(error - buffer)].upper()
    origA = orig[:error].upper()
    length = 0
    for i, c in enumerate(cov[error:], error): # find stop sequence
        if length < k + buffer:
            if check_cov(c, cov_thresh) is False:
                stop = []
            else:
                stop.append(orig[i])
                length = len(stop)
        else:
            break
    stop = ''.join(stop[buffer:]).upper()
    if n_error is False:
        origB = orig[(i - k):]
    else:
        origB = orig[(i - k):n_error]
    origM = orig[error:(i - k)] # mismatched sequence
    return origA, origB, origM, start, stop

def check_overlap_seqs(mA, mB, k):
    """
    see if extensions overlap with one another
    """
    kmer = mA[(len(mA) - k):]
    start = mB.find(kmer)
    if start == -1:
        return False
    return mA + mB[(start + k):]

def patch_contig(seq, cov, errors, scaffold_patches, cov_thresh, k = 10, buffer = 5):
    """
    patch original assembly with re-assembled fragments
    buffer = extra space before/after low coverage region
    k = length of overlap requried for combining re-assembly with original assembly
    merged = [[code, p1], [code, p2], [code, p2], ... [code, pn]]
    codes:
    o = original
    e = extension
    f = fixed
    n = not fixed
    """
    merged = []
    for i, error in enumerate(errors):
        patches = scaffold_patches[error] \
                # error-check re-assembled fragments to replace errors
        # determine the next error
        if (i + 1) >= len(errors):
            n_error = False
        else:
            n_error = errors[i + 1]
        # is the error at the beginning of the scaffold (this is an extension)?
        if error < 100 and patches is not False:
            patch = patch_start(seq[2], cov, error, n_error, patches, k, cov_thresh, buffer)
            if patch is False:
                return False
            merged.extend(patch)
            continue
        origA, origB, origM, start, stop = \
                find_start_stop_error(seq[2], cov, error, n_error, patches, k, cov_thresh, buffer)
        if merged == []: \
                # in case beginning was not added (no beginning error)
            merged.append([error, 'o', origA])
        # check to make sure re-assembly was successful
        if patches is False:
            merged.append([error, 'n', origM])
            merged.append([error, 'o', origB])
            continue
        # is the error at the end of the scaffold (this is an extension)?
        if len(stop) < k:
            patch = \
                    patch_end(seq[2], cov, error, n_error, patches, \
                    k, origA, origB, start, stop, buffer)
            if patch is False:
                continue
            merged.append(patch)
            continue
        # error is in the middle of the scaffold
        middle = \
                patch_middle(seq[2], cov, error, n_error, patches, \
                k, origA, origB, origM, start, stop, buffer)
        if middle is False or middle[0] == 1 or middle[0] == 2: \
                # Failed
            if len(origB) == 0:
                continue
            merged.append([error, 'n', origM])
            merged.append([error, 'o', origB])
            continue
        if middle[0] == 0: \
                # was able to find start and stop sequences in the re-assembly
            merged.append([error, 'f', middle[1]])
            merged.append([error, 'o', origB])
            continue
        if middle[0] == 12: \
                # was able to find the start and the stop, but not on the same re-assembled scaffold
            mA, mB = middle[1][0], middle[1][1]
            combined = check_overlap_seqs(mA, mB, k)
            if combined is False:
                merged.append([error, 'n', origM])
                merged.append([error, 'o', origB])
                continue
            else:
                merged.append([error, 'f', combined])
                merged.append([error, 'o', origB])
                continue
    return merged

def replace_N(merged, min_Ns):
    """
    replace errors that could not be fixed with Ns
    """
    for i, s in enumerate(merged):
        if s[1] == 'n':
            length = len(s[2])
            if length < min_Ns:
                length = min_Ns
            merged[i] = [s[0], s[1], "N"*length]
    return merged

def mask_errors(merged):
    """
    mask errors that could not be fixed
    """
    for i, s in enumerate(merged):
        if s[1] == 'n':
            merged[i] = [s[0], s[1], s[2].lower()]
    return merged

def break_scaf(scaffold, merged, break_scaffolds, ignore_insert, coverage):
    """
    1) if break_scaffolds is True - break all scaffolds where errors could not be fixed
    2) if break_scaffolds is False - break all scaffolds where errors could not be fixed
        only if insert coverage is zero (unless ignore_insert is True)
    """
    errors = [i[0] for i in merged if i[1] == 'n']
    if break_scaffolds is True:
        breaks = errors
    elif ignore_insert is True:
        breaks = []
    else:
        breaks = [i for i in errors if coverage[i][2] is False]
    for i, s in enumerate(merged):
        if s[0] in breaks:
            if s[2] == '':
                del merged[i]
                continue
            if s[1] == 'o': # don't print the error
                merged[i] = [s[0], 'b', '\n>%s_e:%s\n%s' % (scaffold, s[0], s[2])]
            else:
                merged[i] = [s[0], 'b', '']
    return merged

def merge_assemblies(assembly, scaffolds, s2c, s2errors,
        re_assembled, combine_windows, prefix, cov_thresh,
        break_scaffolds, extend_scaffolds, add_Ns, min_Ns, mask, ignore_insert):
    """
    replace errors with re-assembly
    report = [[scaffold, error, code], ... [scaffold, error, code]]
    """
    report = []
    merged_assembly = open('%s/re_assembled.fa' % (prefix), 'w')
    merged_report = open('%s/re_assembled.report.txt' % (prefix), 'w')
    if combine_windows is True:
        for s in re_assembled:
            for w, i in list(re_assembled[s].items()):
                if type(w) is str:
                    for r in w.split('_'):
                        re_assembled[s][int(r)] = i
                    del re_assembled[s][w]
    for id, seq in list(scaffolds.items()):
        if id not in s2errors: # no errors, print scaffold
            report.append([id, 'none', 'n/a'])
            print('\n'.join([seq[0], seq[2]]), file=merged_assembly)
            continue
        errors = sorted(s2errors[id])
        merged = patch_contig(seq, s2c[id], errors, re_assembled[id], cov_thresh)
        if merged is False: # contig has been discarded
            report.append([id, 'n/a', 'removed'])
            continue
        if extend_scaffolds is False:
            merged = [i for i in merged if i[1] != 'e']
        if add_Ns is True:
            merged = replace_N(merged, min_Ns)
        if mask is True:
            merged = mask_errors(merged)
        merged = break_scaf(id, merged, break_scaffolds, ignore_insert, s2c[id])
        # remove empty entries in merged
        merged = [i for i in merged if i[2] != '']
#        # remove scaffolds that are >50% Ns
#        merged = [i for i in merged if (len([j for j in i[2] if j.lower() != "n"])/len(i[2])) > 0.5]
        for s in merged:
            report.append([id, s[0], s[1]])
        print('\n'.join(['>%s' % (id), ''.join([i[2] for i in merged])]), file=merged_assembly)
    for i in report:
        if i[2] == 'o':
            continue
        print('\t'.join([str(j) for j in i]), file=merged_report)
    merged_assembly.close()
    merged_report.close()
    return merged_assembly.name

def add_read(pairs, read, info, r, bit, mate, fastq, mismatches, scaffold):
    """
    add read to paired-read dictionary
    """
    if read not in pairs:
        pairs[read] = [bit, mate, {}, fastq]
    if scaffold not in pairs[read][2]:
        pairs[read][2][scaffold] = []
    pairs[read][2][scaffold].append([r, mismatches, info])
    return pairs

def parse_mapping_stringent(mapping, assembly, mm, \
        ends = False, scaffolds = False, max_cov = 100):
    """
    - create a paired-read dictionary from sam files
    - only include stringently mapped reads
    - use max_cov to limit the number of reads stored
    * pairs[read] = [bit, mate, mappping[scaffold] = [map, map2, ...], fastq]
    *  map = [overlap, mismatches, sam_info]
    *  sam_info = all sam lines except read + quality
    """
    pairs = {}
    header = []
    # make sure that mapping was successful
    if mapping is False:
        return pairs, header
    # s2c[scaffold] = [[coverage per position, # reads connecting base to the next base], [p2, pn2]]
    s2c = {id: [[0, False, False] for i in range(0, info[1])] for id, info in list(assembly.items())}
    for line in open(mapping):
        if line.startswith('@'):
            header.append(line.strip())
            continue
        line = line.strip().split()
        # only include stringently mapped reads
        mismatches = map_tool.count_mismatches(line)
        if mismatches > mm:
            continue
        read, bit, scaffold, start = line[0:4]
        bit, start = int(bit), int(start)
        r = [start, start + len(line[9]) - 1]
        # make sure read is > 10 bp long
        if r[1] - r[0] < 10:
            continue
        fastq = map_tool.sam2fastq(line)
        info = [line[0:9], line[11:]]
        if '/' in read:
            read = read.rsplit('/', 1)[0]
        if bin(bit)[-7] == '1': # first sequence in pair
            read = '%s_1' % (read)
            mate = '%s_2' % (read.rsplit('_', 1)[0])
        else:
            read = '%s_2' % (read)
            mate = '%s_1' % (read.rsplit('_', 1)[0])
        if ends is not False and (r[0] > ends and r[1] < scaffolds[scaffold][1] - ends):
            continue
        if min([i[0] for i in s2c[scaffold][r[0]:r[1]]]) >= max_cov:
            continue
        s2c = add_coverage(scaffold, assembly, r, s2c, line, window = 0)
        pairs = add_read(pairs, read, info, r, bit, mate, fastq, mismatches, scaffold)
    return pairs, header

def parse_mapping_errors(mapping, s2errors, s2windows):
    """
    - create a paired-read dictionary from sam files
    - only include reads mapping to error window
    * pairs[read] = [bit, mate, mappping[scaffold] = [map, map2, ...], fastq]
    *  map = [overlap, mismatches, sam_info]
    *  sam_info = all sam lines except read + quality
    """
    pairs = {}
    for line in open(mapping):
        if line.startswith('@'):
            continue
        line = line.strip().split()
        read, bit, scaffold, start = line[0:4]
        bit, start = int(bit), int(start)
        r = [start, start + len(line[9]) - 1]
        m_scaffold = line[6]
        if scaffold != m_scaffold:
            mate_r = [False, False]
        else:
            mstart = int(line[7])
            mate_r = [mstart, mstart + len(line[9]) - 1]
        # make sure read or mate overlaps with an error window
        if map2window(scaffold, s2windows, s2errors, r, mate_r) is False:
            continue
        mismatches = map_tool.count_mismatches(line)
        fastq = map_tool.sam2fastq(line)
        info = [line[0:9], line[11:]]
        if '/' in read:
            read = read.rsplit('/', 1)[0]
        if bin(bit)[-7] == '1': # first sequence in pair
            read = '%s_1' % (read)
            mate = '%s_2' % (read.rsplit('_', 1)[0])
        else:
            read = '%s_2' % (read)
            mate = '%s_1' % (read.rsplit('_', 1)[0])
        pairs = add_read(pairs, read, info, r, bit, mate, fastq, mismatches, scaffold)
    return pairs

def check_assembly(assembly, pr, threads, cov_thresh, \
        mismatches, collection_mismatches, multiple, \
        prefix, window, combine_windows, \
        pr_split, allow_orphan, allow_orphan_ends, save_mapping):
    """
    identify assembly errors
    """
    # read assembly into memory
    scaffolds = {i[0].split('>')[1]: [i[0], len(i[1]), i[1]] for i in parse_fasta(assembly) if i != [[], []]}
    # map reads to assembly
    mapping, pr_split = map_reads(assembly, scaffolds, pr, threads, multiple, pr_split = pr_split)
    # identfy errors as zero coverage regions based on stringent mapping
    stringent_pairs, header = parse_mapping_stringent(mapping, scaffolds, mismatches)
    # * s2c[scaffold] = [per scaffold coverage at each position]
    # * s2errors = [positions with errors]
    s2c, s2errors = id_errors(stringent_pairs, header, assembly, \
            scaffolds, cov_thresh, mismatches, allow_orphan = allow_orphan, \
            allow_orphan_ends = allow_orphan_ends, save_mapping = save_mapping)
    # define windows
    s2windows = define_windows(scaffolds, s2errors, window, combine_windows)
    error_pairs = parse_mapping_errors(mapping, s2errors, s2windows)
    # collect reads that map to window
    # * reads[scaffold][error] = pe
    reads = collect_reads(error_pairs, assembly, scaffolds, collection_mismatches, prefix, s2errors, s2windows)
    return scaffolds, mapping, s2c, s2errors, reads, pr_split

def format_assembly(fasta, prefix):
    """
    simplify fasta headers and create lookup
    """
    fixed = open('%s/%s' % (prefix, fasta.split('/')[-1]), 'w')
    lookup = open('%s/%s.lookup' % (prefix, fasta.split('/')[-1]), 'w')
    for seq in fix_fasta([fasta], append_index = True, return_original = True):
        original_header, seq = seq
        header = seq[0].split()[0]
        print('\n'.join([header, seq[1].upper()]), file=fixed)
        print('\t'.join(['>%s' % (' '.join(original_header)), header]), file=lookup)
    fixed.close()
    lookup.close()
    return fixed.name

def cleanup(prefix):
    """
    remove intermediate files
    """
    save = ['re_assembled.fa', 're_assembled.report.txt']
    for root, dirs, files in os.walk(prefix, topdown = False):
        save.extend([i for i in files if i.endswith('.lookup')])
        for name in files:
            if name not in save:
                os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))

def curate_assembly(assembly, pr, pr_split, prefix,
        threads = 6, mismatches = 1, collection_mismatches = 2,
        window = 10000, combine_windows = True, cov_thresh = 1,
        multiple = True, assembler = 'velvet',
        break_scaffolds = False, extend_scaffolds = False,
        save_mapping = False, save_int = False, add_Ns = False, min_Ns = 0,
        mask = False, ignore_insert = False):
    """
    identify and correct assembly errors
    """
    os.system('mkdir -p %s' % (prefix))
    assembly = format_assembly(assembly, prefix)
    scaffolds, mapping, s2c, s2errors, mapped_reads, pr_split = \
            check_assembly(assembly, pr, threads, cov_thresh, \
            mismatches, collection_mismatches, multiple, prefix, \
            window, combine_windows, pr_split = pr_split, \
            allow_orphan = False, allow_orphan_ends = False, save_mapping = save_mapping)
    re_assembled = re_assemble_errors(mapped_reads, prefix, pr, pr_split, \
            threads, cov_thresh, mismatches, multiple, assembler = assembler, \
            save_mapping = save_mapping)
    merged = merge_assemblies(assembly, scaffolds, s2c, s2errors, re_assembled, \
            combine_windows, prefix, cov_thresh, break_scaffolds, extend_scaffolds, \
            add_Ns, min_Ns, mask, ignore_insert)
    if save_int is False:
        cleanup(prefix)

def find_reads(one, two, inter, js, reads, read_list):
    """
    define reads based on input file names or supplied json
    """
    if js is not False:
        s2reads = {}
        with open(js) as handle:
            json_data = json.load(handle)
        for i in json_data:
            s2reads.update(json_data[i]['reads'])
    if read_list is True:
        print('samples:')
        samples = ['   %s' % (i) for i in list(s2reads.keys())]
        print('\n'.join(samples))
        exit()
    if js is not False and reads is False:
        print('specify reads with -reads', file=sys.stderr)
        print('use --read_list for list of samples in %s' % (js), file=sys.stderr)
        exit()
    if js is not False:
        for sample in reads.split(','):
            if sample not in s2reads:
                print('%s not in %s' % (sample, js))
                print('use --read_list for list of samples in %s' % (js))
                exit()
        one = ','.join([s2reads[i][0] for i in reads.split(',')])
        two = ','.join([s2reads[i][1] for i in reads.split(',')])
    elif (one is False and two is not False) or \
            (one is False and inter is False):
        print('specify reads with -1 and -2, -12, or -json', file=sys.stderr)
        exit()
    return one, two, inter

def check_previous(directory, overwrite, keep):
    """
    check to see if output directory exists
    - if it does
        - remove it if -f is provided (overwrite)
        - leave it if -ff is provied (keep)
    """
    if os.path.exists(directory) is False:
        return
    if overwrite is False and keep is False:
        print('output directory found: %s' % (directory), file=sys.stderr)
        print('use -f to overwrite (recommended) or -ff to use existing files (not recommended)', file=sys.stderr)
        exit()
    if overwrite is True and keep is True:
        print('output directory found: %s' % (directory), file=sys.stderr)
        print('use -f to overwrite (recommended) or -ff to use existing files (not recommended)', file=sys.stderr)
        exit()
    if overwrite is True:
        shutil.rmtree(directory, ignore_errors = True)
        return
    if keep is True:
        return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '# curate assembly errors')
    parser.add_argument(\
            '-i', required = True, \
            help = 'fasta file to curate')
    parser.add_argument(\
            '-1', required = False, default = False, \
            help = 'fastq for the first read-pair')
    parser.add_argument(\
            '-2', required = False, default = False, \
            help = 'fastq for the second read-pair')
    parser.add_argument(\
            '-12', required = False, default = False, \
            help = 'fastq for interleaved paired-reads (not recommended; generates split fastq files)')
    parser.add_argument(\
            '-json', required = False, default = False, \
            help = 'json file with paths to read data')
    parser.add_argument(\
            '-reads', required = False, default = False, \
            help = 'comma separated list of sample names')
    parser.add_argument(\
            '--read_list', action = 'store_true', \
            help  = 'list sample names from json file')
    parser.add_argument(\
            '-m', required = False, default = 1, type = int, \
            help = 'number of mismatches tolerated during mapping (default = 1)')
    parser.add_argument(\
            '-c', required = False, default = 2, type = int, \
            help = 'number of mismatches tolerated during re-assembly read recruitment (applies to one read in pair; default = 2)')
    parser.add_argument(\
            '-w', required = False, default = 1000, type = int, \
            help = 'size of re-assembly window around errors (default = 1000)')
    parser.add_argument(\
            '--mask', action = 'store_true', \
            help = 'mask errors that could not be corrected')
    parser.add_argument(\
            '--add-Ns', action = 'store_true', \
            help = 'replace errors that could not be corrected with Ns')
    parser.add_argument(\
            '-N', required = False, default = 50, type = int, \
            help = 'min. number of Ns for filling gaps (default = 50)')
    parser.add_argument(\
            '--break', action = 'store_true', \
            help = 'break scaffold if error could not be corrected')
    parser.add_argument(\
            '--extend', action = 'store_true', \
            help = 'extend ends of scaffolds')
    parser.add_argument(\
            '--ignore-insert-cov', action = 'store_true', \
            help = 'do not break scaffolds when insert coverage is zero')
    parser.add_argument(\
            '--mapping', action = 'store_true', \
            help = 'save filtered intermediate sam files (not recommended)')
    parser.add_argument(\
            '--save_int', action = 'store_true', \
            help = 'save intermediate files (not recommended)')
    parser.add_argument(\
            '-f', action = 'store_true', \
            help = 'force overwrite of output directory')
    parser.add_argument(\
            '-ff', action = 'store_true', \
            help = 'force continue from previous run (not recommended)')
    parser.add_argument(\
            '-t', required = False, default = 6, type = int, \
            help = 'number of threads (default = 6)')
    args = vars(parser.parse_args())
    one, two, inter = find_reads(args['1'], args['2'], args['12'], \
            args['json'], args['reads'], args['read_list'])
    fasta = args['i']
    window = args['w']
    prefix = '%s.curated' % (fasta.split('/')[-1].rsplit('.', 1)[0])
    check_previous(prefix, args['f'], args['ff'])
    curate_assembly(fasta, inter, [one, two], prefix,
            threads = args['t'], mismatches = args['m'],
            collection_mismatches = args['c'], break_scaffolds = args['break'],
            extend_scaffolds = args['extend'], save_mapping = args['mapping'],
            save_int = args['save_int'], add_Ns = args['add_Ns'], min_Ns = args['N'],
            mask = args['mask'], ignore_insert = args['ignore_insert_cov'],
            window = window)
