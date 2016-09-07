#!/usr/bin/env python3

"""
script for fixing scaffolding errors by:
 1. identifying errors based on stringent read mapping (stringent mapping for both reads in a pair)
 2. collecting mapped reads (stringent mapping for one read in pair)
 3. re-assembling collected reads
 4. merging the new assembly with the old assembly
 5. doing a final check of the re-assembled scaffolds
(this could be made iterative, or combined with a re-assembly step using the final assembly and 
 reads mapping to the final assembly - stringent mapping for one read in the pair, ideally)
"""

import sys
import os
import itertools
from glob import glob as glob
import subprocess
import re

# ctb
sys.path.append('/home/cbrown/bioscripts/bin')
import mapped as map_tool
from fastq_split import split as fastq_split
from fix_fasta import fix_fasta as fix_fasta
from fasta import iterate_fasta as parse_fasta
from assemble import velvet as velvet
from rc import reverse_complement as rc

def fastq2fasta(fastq, paired = True):
    """
    convert fastq file to fasta
    """
    c = itertools.cycle([1, 2, 3, 4])
    p = itertools.cycle([1, 2])
    for line in open(fastq):
        n = next(c)
        if n == 1:
            if paired is True:
                s = ['%s/%s' % (line.strip().replace('@', '>'), next(p))]
            else:
                s = [line.strip().replace('@', '>')]
        elif n == 2:
            s.append(line.strip())
            yield s

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
        elif pr is not False:
            base = '%s/%s' % (bt_dir, pr.rsplit('/', 1)[-1].rsplit('.', 1)[0])
            r1 = '%s.R1.fastq' % (base)
            r2 = '%s.R2.fastq' % (base)
            pr_split = [r1, r2]
            if os.path.exists(r1) is False:
                fastq_split(open(pr), base)
            pr_command = '-1 %s -2 %s' % (r1, r2)
        if multiple is True: # should reads be allowed to map to more than one position or contig?
            matches_command = '-a '
#        print 'bowtie2 --very-fast --reorder --quiet %s -p %s -x %s/%s %s %s | shrinksam -k %s' \
#                % (matches_command, threads, bt_dir, assembly.rsplit('/', 1)[-1], pr_command, sr_command, sam)
        p.communicate()

        p = subprocess.Popen('\
                bowtie2 --very-fast --reorder --quiet %s -p %s -x %s/%s %s %s | shrinksam -k %s' \
                % (matches_command, threads, bt_dir, assembly.rsplit('/', 1)[-1], pr_command, sr_command, sam)\
                , shell = True)
        p.communicate()
    return sam, pr_split

def map_reads(assembly, scaffolds, pr, threads, multiple, bt_dir = 'bt2', pr_split = False, sr = False):
    """
    map reads to assembly with bowite2
    """
    if '/' in assembly:
        bt_dir = '%s/%s' % (assembly.rsplit('/', 1)[0], bt_dir)
        assembly_dir, assembly_name = assembly.rsplit('/', 1)
    else:
        assembly_dir, assembly_name = './', assembly
    os.system('mkdir -p %s' % (bt_dir))
    sam = '%s/%s.bt.sam' % (bt_dir, assembly_name.rsplit('.', 1)[0])
    return run_bowtie(assembly, sam, pr, pr_split, sr, threads, multiple, bt_dir) # run bowtie, return sam file

def check_cov(c, cov_thresh):
    """
    check for both read coverage and passing overlap test
    c = [coverage, overlap = True or False]
    """
    if c[0] >= cov_thresh and c[1] is True:
        return True
    return False

def errors_from_cov(s2c, cov_thresh):
    """
    identify errors from low coverage region
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
    if mismatches are not in the beginning or end of the read, return False
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
    cov_info = [coverage, overlap_coverage]
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
    return s2c

def id_errors(pairs, header, assembly, scaffolds, cov_thresh, mismatches, allow_orphan = False, allow_orphan_ends = False, window = 3, orphan_window = 1000):
    """
    identify regions with zero coverage by stringently mapped paired reads
    * pairs[read] = [bit, mate, mappping[scaffold] = [map, map2, ...], fastq]
    *  map = [overlap, mismatches, sam_info]
    *  sam_info = all sam lines except read + quality
    """
    # open sam file for writing
    out = open('%s.both.sam' % (assembly.rsplit('.', 1)[0]), 'w')
    for line in header:
        print(line, file=out)
    # s2c[scaffold] = [[coverage per position, # reads connecting base to the next base], [p2, pn2]]
    s2c = {id: [[0, False] for i in range(0, info[1])] for id, info in list(scaffolds.items())}
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
                    if (allow_orphan is False and allow_orphan_ends is False) and (mm is False or mate_mm is False): # make sure both pairs mapped
                        continue
                    elif allow_orphan_ends is True and (mm is False or mate_mm is False):
                        pos = [sam_info[0][3], mate_sam[0][3]]
                        if min(pos) > orphan_window and max(pos) < scaffolds[scaffold][1] - orphan_window:
                            continue
                    if mm > mismatches or mate_mm > mismatches: # both pairs must pass mismatches
#                    if mm + mate_mm > mismatches: # total number of mismatches must be less than mismatches
                        continue
                    sam = mapping[2][0] + [read[3][1]] + [read[3][3]] + mapping[2][1]
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
                if stop > scaffolds[scaffold]:
                    stop = scaffolds[scaffold]
            else:
                stop = (error + int(window/2))
                if stop > scaffolds[scaffold]:
                    stop = scaffolds[scaffold]
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

def collect_reads(pairs, assembly, scaffolds, mismatches, prefix, s2errors = False, window = False, combine_windows = False):
    """
    collect reads mapped to scaffold within window, or for the entire contig
    - only one read in pair has to pass mismatch critera
    * pairs[read] = [bit, mate, mappping[scaffold] = map, fastq]
       map = [overlap, mismatches, sam_info]
       sam_info = all sam lines except read + quality
    * reads[scaffold][error] = [[sequences-fastq], filename]
    """
    # define windows
    if window is not False:
        s2windows = define_windows(scaffolds, s2errors, window, combine_windows)
    else:
        s2windows = False
    # make files for reads
    reads = {} # reads[scaffold][error] = pe
    if window is not False:
        for scaffold, errors in list(s2windows.items()):
            if len(errors) == 0:
                continue
            reads[scaffold] = {}
            for error in errors:
                dir = '%s/s_%s_e_%s' % (prefix, scaffold, error)
                os.system('mkdir -p %s' % (dir))
                reads[scaffold][error] = [{}, '%s/reads.pe.fastq' % (dir)]
    else:
        reads[False] = {}
        reads[False][False] = [{}, '%s.one.pe.fastq' % (assembly.rsplit('.', 1)[0])]
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
                    if window is False:
                        errors = [False]
                        scaffold = False
                    else:
                        errors = map2window(scaffold, s2windows, s2errors, overlap, m_overlap)
                    if window is not False and errors is False:
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
            for seq in list(seqs.values()):
                print(seq, file=out)
            out.close()
            reads[scaffold][error] = out.name
    return reads

def break_by_coverage(assembled, s2c, cov_thresh, ignore_ends = False):
    """
    look through re-assembled contigs and break them where stringent coverage is low
    """
    fragments = [] # fragments = [[[len, header], [sequence]], ...]
    if ignore_ends is False:
        for id, seq in list(assembled.items()):
            sequence = []
            for i, base in enumerate(s2c[id]):
                if check_cov(base, cov_thresh) is True:
                    if sequence == []:
                        start = i
                    sequence.append(seq[2][i])
                else:
                    if sequence != []:
                        fragments.append([[len(sequence), '>%s_f:%s' % (id, start)], ''.join(sequence)])
                    sequence = []
            if sequence != []:
                fragments.append([[len(sequence), '>%s_f:%s' % (id, start)], ''.join(sequence)])
    else:
        for id, seq in list(assembled.items()):
            sequence = []
            for i, base in enumerate(s2c[id]):
                if check_cov(base, cov_thresh) is True or (i < 100 or i > (seq[1] - 100)):
                    if sequence == []:
                        start = i
                    sequence.append(seq[2][i])
                else:
                    if sequence != []:
                        fragments.append([[len(sequence), '>%s_f:%s' % (id, start)], ''.join(sequence)])
                    sequence = []
            if sequence != []:
                fragments.append([[len(sequence), '>%s_f:%s' % (id, start)], ''.join(sequence)])
    fragments.sort(key = lambda x: x[0][0], reverse = True)
    return fragments

def re_assemble_velvet(pr, prefix, scaffold, error, scaffolding, min_contig):
    """
    re-assemble reads using velvet
    """
    out = '%s/s_%s_e_%s' % (prefix, scaffold, error)
    velvet(paired = [pr], out = out, scaffolding = scaffolding, \
        silent = True, min_contig = min_contig, kmer_min = 21, kmer_max = 71, kmer_increase = 10)
    re_assembled_seqs = {}
    assembled_fasta = open('%s/velvet_s-%s.e_%s.fa' % (out, scaffold, error), 'w')
    for fasta in glob('%s/*.fasta' % (out)):
        for seq in parse_fasta(fasta):
            if seq[0] == []:
                continue
            re_assembled_seqs[seq[0].split('>')[1]] = [seq[0], len(seq[1]), seq[1]]
            print('\n'.join(seq), file=assembled_fasta)
    assembled_fasta.close()
    return assembled_fasta.name, re_assembled_seqs

def re_assemble_shorty(pr, prefix, scaffold, error, min_contig):
    """
    re-assemble reads using shorty
    """
    out = '%s/s_%s_e_%s' % (prefix, scaffold, error) # shorty
    os.system('mkdir -p %s' % out) # shorty
    assembled_fasta = '%s/shorty_s-%s.e-%s.fa' % (out, scaffold, error) # shorty
    p = subprocess.Popen('shorty -s %s -q %s | fix_fasta.py - > %s' \
        % (min_contig, pr, assembled_fasta), shell = True)
    p.communicate()
    re_assembled_seqs = {seq[0].split('>')[1]: [seq[0], len(seq[1]), seq[1]] for seq in parse_fasta(assembled_fasta)}
    return assembled_fasta, re_assembled_seqs
 
def re_assemble_minimo(pr, prefix, scaffold, error, min_contig):
    """
    re-assemble reads using shorty
    """
    current = os.getcwd()
    out = '%s/s_%s_e_%s' % (prefix, scaffold, error)
    os.system('mkdir -p %s' % out)
    assembled_fasta = '%s/minimo_s_%s.e_%s.fa' % (out, scaffold, error)
    if os.path.exists(assembled_fasta):
        re_assembled_seqs = {seq[0].split('>')[1]: [seq[0], len(seq[1]), seq[1]] for seq in parse_fasta(assembled_fasta) if seq[0] is not list}
        return assembled_fasta, re_assembled_seqs
    fasta = open('%s.fa' % (pr.rsplit('.', 2)[0]), 'w')
    for seq in fastq2fasta(pr, True):
        print('\n'.join(seq), file=fasta)
    fasta.close()
    p = subprocess.Popen('cd %s; Minimo reads.fa \
            -D MIN_LEN=50 -D MIN_IDENT=98 -D FASTA_EXP=1 -D OUT_PREFIX=minimo -D ACE_EXP=0 > minimo.log; \
            fix_fasta.py minimo-contigs.fa | fasta_length.py - %s > %s; cd %s' \
        % (out, min_contig, assembled_fasta.rsplit('/', 1)[1], current), shell = True)
    p.communicate()
    re_assembled_seqs = {seq[0].split('>')[1]: [seq[0], len(seq[1]), seq[1]] for seq in parse_fasta(assembled_fasta) if type(seq[0]) is str}
    return assembled_fasta, re_assembled_seqs

def idba_ud_seqs(fasta_name, idba_dir):
    """
    get fasta file from idba_ud assembly directory, save to dictionary
    """
    out = '%s/scaffold.fa' % (idba_dir)
    if os.path.exists(out) is False:
        out = '%s/contig.fa' % (idba_dir)
    re_assembled_seqs = \
            {seq[0].split('>')[1].replace(' ', '_'): [seq[0].replace(' ', '_'), len(seq[1]), seq[1]] \
            for seq in parse_fasta(out) if type(seq[0]) is str}
    os.system('cat %s | tr " " "_" > %s' % (out, fasta_name))
    return re_assembled_seqs

def re_assemble_idba_ud(pr, prefix, scaffold, error, min_contig, threads):
    """
    re-assemble reads using idba_ud
    """
    current = os.getcwd()
    out = '%s/s_%s_e_%s' % (prefix, scaffold, error)
    idba_out = '%s/%s' % (out, 'idba_ud')
    os.system('mkdir -p %s' % out)
    assembled_fasta = '%s/idba_ud_s_%s.e_%s.fa' % (out, scaffold, error)
    if os.path.exists(idba_out):
        re_assembled_seqs = idba_ud_seqs(assembled_fasta, idba_out)
        return assembled_fasta, re_assembled_seqs
    fasta = open('%s.fa' % (pr.rsplit('.', 2)[0]), 'w')
    for seq in fastq2fasta(pr, True):
        print('\n'.join(seq), file=fasta)
    fasta.close()
    p = subprocess.Popen('idba_ud --step 10 --pre_correction --similar 0.98  -r %s -o %s --num_threads %s --min_contig %s \
            >> idba_ud.log 2>> idba_ud.log' \
        % (fasta.name, idba_out, threads, min_contig), shell = True)
    out, error = p.communicate()
    return assembled_fasta, idba_ud_seqs(assembled_fasta, idba_out)
  
def re_assemble_fragments(reads, prefix, threads, cov_thresh, mismatches, multiple, scaffolding = True, assembler = 'velvet', min_contig = '200'):
    """
    re-assemble misassembled regions, return fragments (sorted longest to shortest) representative of re-assembled
    regions with coverage by stringently mapped reads
    * assembler == 'velvet' or 'shorty'
    """
    re_assembled = {}
    for scaffold, errors in list(reads.items()):
        re_assembled[scaffold] = {}
        for error, pr in list(errors.items()):
            if assembler == 'velvet':
                assembled_fasta, re_assembled_seqs = \
                        re_assemble_velvet(pr, prefix, scaffold, error, scaffolding, min_contig)
            elif assembler == 'shorty':
                assembled_fasta, re_assembled_seqs = \
                        re_assemble_shorty(pr, prefix, scaffold, error, min_contig)
            elif assembler == 'minimo':
                assembled_fasta, re_assembled_seqs = \
                        re_assemble_minimo(pr, prefix, scaffold, error, min_contig)
            elif assembler == 'idba_ud':
                assembled_fasta, re_assembled_seqs = \
                        re_assemble_idba_ud(pr, prefix, scaffold, error, min_contig, threads)
            else:
                print('# please specify valid assembler', file=sys.stderr)
                exit()
            if len(re_assembled_seqs) == 0:
                re_assembled[scaffold][error] = False
                continue
            # map reads to re-assembled scaffolds and identify zero coverage regions
            mapping, pr_split = map_reads(assembled_fasta, re_assembled_seqs, pr, threads, multiple)
            pairs, header = parse_mapping(mapping)
            s2c, s2errors = id_errors(pairs, header, assembled_fasta, re_assembled_seqs, cov_thresh, mismatches, allow_orphan = True)
            fragments = break_by_coverage(re_assembled_seqs, s2c, cov_thresh, ignore_ends = False)
            f_out = open('%s/s_%s_e_%s/fragments.fa' % (prefix, scaffold, error), 'w')
            for f in fragments:
                print('\n'.join([f[0][1], f[1]]), file=f_out)
            f_out.close()
            re_assembled[scaffold][error] = fragments
    return re_assembled

def patch_start(orig, cov, error, n_error, patches, merked, k, cov_thresh, buffer):
    """
    use re-assembled fragments to extend scaffold from start
    """
    length = 0
    trimmed = False
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
        return [best, trimmed]
    else: # what if you could not extend the scaffold?
        return [trimmed]

def patch_end(orig, cov, error, n_error, patches, merged, k, origA, origB, start, stop, buffer):
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
        return best
    else: # what if you could not extend the scaffold?
        return False

def patch_middle(orig, cov, error, n_error, patches, merged, k, origA, origB, origM, start, stop, buffer):
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
 
def find_start_stop_error(orig, cov, error, n_error, patches, merged, k, cov_thresh, buffer):
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

def patch_contig(orig, cov, error, n_error, patches, merged, scaffold, cov_thresh, k = 10, buffer = 5):
    """
    patch original assembly with re-assembled fragment
    buffer = extra space before/after low coverage region
    """
    # check to see if error is at the beginning of scaffold, if so - extend it
    if error == 0 and patches is not False:
        merged = patch_start(orig, cov, error, n_error, patches, merged, k, cov_thresh, buffer)
        return merged
    origA, origB, origM, start, stop = find_start_stop_error(orig, cov, error, n_error, patches, merged, k, cov_thresh, buffer)
    if merged == []:
        merged = [origA] # if not, origA is already there
    # check to make sure re-assembly was successful
    if patches is False:
        merged.append('\n>%s_e:%s\n' % (scaffold, error))
        merged.append(origB)
        return merged
    if len(stop) < k: # error is at the end of the scaffold
        extended = patch_end(orig, cov, error, n_error, patches, merged, k, origA, origB, start, stop, buffer)
        if extended is not False and merged is not False:
            merged.append(extended)
            return merged
        else:
            return merged
    # error is in the middle of the scaffold
    middle = patch_middle(orig, cov, error, n_error, patches, merged, k, origA, origB, origM, start, stop, buffer)
    if middle is False: # was not able to find start or stop sequence in re-assembly
        if len(origB) > 0:
            merged.append('\n>%s_e:%s\n' % (scaffold, error))
            merged.append(origB)
        return merged
    if middle[0] == 0: # was able to find start and stop sequences in the re-assembly
        middle = middle[1]
        merged.append(middle)
        merged.append(origB)
        print('\t\t fixed')
        return merged
    if middle[0] == 1: # was only able to find the start (not the stop) sequence in the re_assembly
        middle = middle[1]
        merged.append(middle)
        if len(origB) > 0:
            merged.append('\n>%s_e:%s\n' % (scaffold, error))
            merged.append(origB)
        return merged
    if middle[0] == 2: # was only able to find the stop (not the start) sequence in the re_assembly
        middle = middle[1]
        merged.append('\n>%s_e:%s\n' % (scaffold, error))
        merged.append(middle)
        if len(origB) > 0:
            merged.append(origB)
        return merged
    if middle[0] == 12: # was able to find the start and the stop, but not on the same re-assembled scaffold
        mA, mB = middle[1][0], middle[1][1]
        combined = check_overlap_seqs(mA, mB, k)
        if combined is False:
            merged.append(mA)
            merged.append('\n>%s_e:%s\n' % (scaffold, error))
            merged.append(mB)
            if len(origB) > 0:
                merged.append(origB)
            return merged
        else:
            merged.append(combined)
            print('\t\t fixed')
            return merged
        
def merge_assemblies(assembly, scaffolds, s2c, s2errors, re_assembled, combine_windows, prefix, cov_thresh):
    """
    merge asseembly and re-assembly
    """
    merged_assembly = open('%s/re_assembled-draft.fa' % (prefix), 'w')
    if combine_windows is True:
        for s in re_assembled:
            for w, i in list(re_assembled[s].items()):
                if type(w) is str:
                    for r in w.split('_'):
                        re_assembled[s][int(r)] = i
                    del re_assembled[s][w]
    for id, seq in list(scaffolds.items()):
        print(id)
        merged = []
        if id not in s2errors:
            print('\n'.join([seq[0], seq[2]]), file=merged_assembly)
        errors = sorted(s2errors[id])
        print('\terrors: %s' % (errors))
        for i, error in enumerate(errors):
            print('\terror: %s' % error)
            if (i + 1) >= len(errors):
                n_error = False
            else:
                n_error = errors[i + 1]
            merged = patch_contig(seq[2], s2c[id], error, n_error, re_assembled[id][error], merged, id, cov_thresh)
        if merged is not False:
            merged.insert(0, '>%s\n' % id)
            print(''.join(merged), file=merged_assembly)
    merged_assembly.close()
    return merged_assembly.name

def final_check(assembly, pr_split, threads, cov_thresh, mismatches, collection_mismatches, prefix, multiple):
    """
    - check re-assembled contigs for zero coverage regions, split at these points
    - re-map to final version
    """
    scaffolds, mapping, s2c, s2errors, reads, pr_split = \
            check_assembly(assembly, False, threads, cov_thresh, mismatches, collection_mismatches, multiple, prefix, pr_split = pr_split, allow_orphan = False, allow_orphan_ends = True) 
    final_fasta = open('%s/re_assembled-final.fa' % (prefix), 'w')
    for seq in break_by_coverage(scaffolds, s2c, cov_thresh, ignore_ends = False):
        print('\n'.join([seq[0][1], seq[1]]), file=final_fasta)
    final_fasta.close()
    # check final assembly
    return check_assembly(final_fasta.name, False, threads, cov_thresh, mismatches, collection_mismatches, multiple, prefix, pr_split = pr_split, allow_orphan = False, allow_orphan_ends = True)

def parse_mapping(mapping, ends = False, scaffolds = False):
    """
    create a paired-read dictionary from sam files
    * pairs[read] = [bit, mate, mappping[scaffold] = [map, map2, ...], fastq]
    *  map = [overlap, mismatches, sam_info]
    *  sam_info = all sam lines except read + quality
    """
    pairs = {}
    header = []
    for line in open(mapping):
        if line.startswith('@'):
            header.append(line.strip())
            continue
        line = line.strip().split()
        read, bit, scaffold, start = line[0:4]
        bit, start = int(bit), int(start)
        r = [start, start + len(line[9]) - 1]
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
        if ends is not False and (r[0] > ends and r[1] < scaffolds[scaffold][1] - ends):
            continue
        if read not in pairs:
            pairs[read] = [bit, mate, {}, fastq]
        if scaffold not in pairs[read][2]:
            pairs[read][2][scaffold] = []
        pairs[read][2][scaffold].append([r, mismatches, info])
    return pairs, header

def check_assembly(assembly, pr, threads, cov_thresh, mismatches, collection_mismatches, multiple, prefix, window = False, combine_windows = False, pr_split = False, allow_orphan = False, allow_orphan_ends = False):
    """
    check assembly for mismatches
    """
    # read assembly into memory
    scaffolds = {i[0].split('>')[1]: [i[0], len(i[1]), i[1]] for i in parse_fasta(assembly) if i != [[], []]}
    # map reads to assembly
    mapping, pr_split = map_reads(assembly, scaffolds, pr, threads, multiple, pr_split = pr_split)
    # identfy errors (no coverage regions) based on stringent mapping
    # * s2errors = [positions with errors]
    # * s2c[scaffold] = [per scaffold coverage at each position]
    # * filtered_mapping_both = sam file requiring both reads to map within mismatch criteria
    pairs, header = parse_mapping(mapping)
    s2c, s2errors = id_errors(pairs, header, assembly, scaffolds, cov_thresh, mismatches, allow_orphan = allow_orphan, allow_orphan_ends = allow_orphan_ends)
    # collect reads that map to either window or entire contigs
    # * filtered_mapping_one = sam file requiring one of paired reads to map within mismatch criteria
    # * reads[scaffold][error] = pe
    if window is False:
        reads = collect_reads(pairs, assembly, scaffolds, collection_mismatches, prefix)
    else:
        reads = collect_reads(pairs, assembly, scaffolds, collection_mismatches, prefix, s2errors, window, combine_windows)
    return scaffolds, mapping, s2c, s2errors, reads, pr_split

def format_assembly(fasta, prefix):
    """
    get rid of strange characters in assembly, make new file in curated directory
    """
    if '/' in fasta:
        fixed = open('%s/%s' % (prefix, fasta.rsplit('/', 1)[1]), 'w')
    else:
        fixed = open('%s/%s' % (prefix, fasta), 'w')
    for seq in fix_fasta(fasta):
        print('\n'.join(seq), file=fixed)
    fixed.close()
    return fixed.name

def curate_assembly(assembly, pr, pr_split, prefix, \
        threads, mismatches = 1, collection_mismatches = 2, window = 1000, combine_windows = True, cov_thresh = 1, multiple = True, assembler = 'velvet'):
    """
    identify and correct scaffolding errors in assembly
    """
    os.system('mkdir -p %s' % (prefix))
    assembly = format_assembly(assembly, prefix)
    scaffolds, mapping, s2c, s2errors, mapped_reads, pr_split = \
            check_assembly(assembly, pr, threads, cov_thresh, mismatches, collection_mismatches, multiple, prefix, window, combine_windows, pr_split = pr_split, allow_orphan = False, allow_orphan_ends = False) 
    re_assembled = re_assemble_fragments(mapped_reads, prefix, threads, cov_thresh, mismatches, multiple, assembler = assembler)
    merged = merge_assemblies(assembly, scaffolds, s2c, s2errors, re_assembled, combine_windows, prefix, cov_thresh)
    return final_check(merged, pr_split, threads, cov_thresh, mismatches, collection_mismatches, prefix, multiple)

if __name__ == '__main__':
    if len(sys.argv) != 7:
        print('specify threads, fasta for assembly, mismatches tolerated, paired reads (fastq), pair 1 (fastq), and pair 2 (fastq)', file=sys.stderr)
        exit()
    for i, v in enumerate(sys.argv[3:], 3):
        if v == 'False' or v == 'FALSE' or v == 'false':
            sys.argv[i] = False
    threads, assembly, mm, pr, pr1, pr2 = sys.argv[1:]
    if pr1 is False or pr2 is False:
        pr_split = False
    else:
        pr_split = [pr1, pr2]
    if '/' in assembly:
        prefix = '%s.curated' % (assembly.rsplit('.', 1)[0].rsplit('/', 1)[1])
    else:
        prefix = '%s.curated' % (assembly.rsplit('.', 1)[0])
    curate_assembly(assembly, pr, pr_split, prefix, mismatches = int(mm), threads = int(threads))
