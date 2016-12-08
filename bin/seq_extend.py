#!/usr/bin/env python3

"""
extend contig(s) and/or fill gap(s)
"""

import os
import sys
import argparse
import networkx as nx
from tqdm import tqdm
import pickle as Pickle
from itertools import cycle
from multiprocessing import Pool
from numpy import random as Random
from matplotlib import pyplot as plt
from Levenshtein import distance as ldist

# ctb
from fasta import iterate_fasta as parse_fasta

def parse_fastq(fastq):
    """
    parse fastq file and return list for each sequence
    containing each line describing the sequence
    """
    c = cycle([1, 2, 3, 4])
    for line in fastq:
        n = next(c)
        if n == 1:
            s = []
        s.append(line.strip())
        if n == 4:
            yield s

def parse_readsRC(reads_file):
    """
    parse reads file and return both
    the sequence and its reverse complement
    """
    COMP = dict(zip(['A', 'T', 'G', 'C', 'N', '-'], ['T', 'A', 'C', 'G', 'N', '-']))
    for line in reads_file:
        if line.startswith('>'):
            file_type = 'fasta'
            break
        file_type = 'fastq'
        break
    reads_file.seek(0)
    if file_type == 'fasta':
        for seq in parse_fasta(reads_file):
            # return sequence and reverse complement
            seq = seq[1].upper()
            yield seq
            yield ''.join([COMP[base] for base in seq[::-1]])
    else:
        for seq in parse_fastq(reads_file):
            # return sequence and reverse complement
            seq = seq[1]
            yield seq
            yield ''.join([COMP[base] for base in seq[::-1]])

def trim_seqs(seqs, length):
    """
    trim sequence next to Ns
    """
    # remove first l bases from all but the first fragment
    # remove last l bases from all but the last fragment
    trimmed = []
    for i, seq in enumerate(seqs):
        if i == 0:
            seq = seq[0:-length]
        elif i == len(seqs) - 1:
            seq = seq[length:]
        else:
            seq = seq[length:-length]
        trimmed.append(seq)
    return trimmed

def parse_gaps(contig, pars):
    """
    parse sequence into gapped regions:
        return [seqA, gap, seqB], [seqB, gap, seqA], ....
    """
    # return False if there are no gaps
    if contig.find('N') == -1:
        return False
    seqs, gaps = [], []
    next_seq = True
    for base in contig:
        if base != 'N':
            if next_seq is True:
                seqs.append('')
                gaps.append('')
                next_seq = False
            seqs[-1] += base
        else:
            gaps[-1] += base
            next_seq = True
    # remove bases close to Ns because they could be incorrect
    seqs = trim_seqs(seqs, pars['buffer'])
    fragments = []
    for seq, gap in zip(seqs, gaps):
        if len(seq) < pars['min_overlap']:
            continue
        fragments.append(seq)
        fragments.append(gap)
    if fragments[-1] == '':
        fragments.pop(-1)
    return fragments

def add_sequence(readG, seq, position, pars):
    """
    add sequence to read graph
    """
    add = []
    if position == 'left':
        for read in readG:
            add.append((seq, read))
    elif position == 'right':
        for read in readG:
            add.append((read, seq))
    elif position == 'both':
        for read in readG:
            add.append((seq, read))
            add.append((read, seq))
    readG.add_node(seq)
    pool = Pool(pars['threads'])
    num_compare = len(add)
    for pair in pool.imap_unordered(align_pair, \
            [(pair, pars) for pair in add]):
        if pair is not False:
            score, a, b, extension, len_extension, overlap = pair
            readG.add_edge(a, b, score = score, \
                    extension = extension, len_extension = len_extension, \
                    overlap = overlap, inv_overlap = (1/overlap))
    return readG

def draw_net(G, node_list, edge_list, source, target):
    """
    draw network graph
    """
    # draw network
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, nodelist = node_list, node_color = 'b', \
                            node_size = 500, alpha = 0.5)
    nx.draw_networkx_nodes(G, pos, nodelist = [source], node_color = 'g', \
                            node_size = 750, alpha = 0.5)
    nx.draw_networkx_nodes(G, pos, nodelist = [target], node_color = 'r', \
                            node_size = 750, alpha = 0.5)
    nx.draw_networkx_edges(G, pos, edgelist = edge_list)
    plt.show()

def shortest_path(readG, source, target, trim = True):
    """
    get sequence corresponding with path through reads graph
    """
    path   = nx.shortest_path(readG, source = source, target = target, weight = 'inv_overlap')
    edges  = [(path[i], path[i+1]) for i, node in enumerate(path) if i + 1 < len(path)]
    if trim is True:
        seq = ''.join([readG[a][b]['extension'] for a, b in edges]).rsplit(target, 1)[0]
    else:
        seq = ''.join([readG[a][b]['extension'] for a, b in edges])
    length = len(seq)
    return path, edges, seq, length

def subset_graph(G, source, target):
    """
    cut graph at source and target to only
    include nodes in between
    """
    subG = G.copy()
    ## remove connections to source node
    for a, b in subG.edges(nbunch = source):
        if a != source:
            subG.remove_edge(a, b)
    ## remove connection from target
    for a, b in subG.edges(nbunch = target):
        if b != target:
            subG.remove_edge(a, b)
    ## remove other nodes from graph
    keep = set(nx.node_connected_component(subG.to_undirected(), source))
    remove = set(subG.nodes()).difference(keep)
    subG.remove_nodes_from(remove)
    return subG

def reverse_graph(readG):
    """
    reverse graph and re-determine the
    extension
    """
    # reverse the graph
    RreadG = nx.DiGraph()
    for edge in readG.edges():
        a, b = edge[0], edge[1]
        extension  = readG[a][b]['extension']
        # determine extension for reverse direction
        Rlength_extension = len(a) - b.rfind(extension)
        Rextension = a[0:Rlength_extension]
        RreadG.add_edge(b, a, extension = Rextension, length_extension = Rlength_extension)
    return RreadG

def find_consensus(readG, source, target, trim = True):
    """
    find consensus sequence between two
    sequences in directed read graph
    """
    # use shortest path as scaffold
    path, path_edges, path_seq, path_len = shortest_path(readG, source, target, trim = trim)
    bases      = {pos:{'A':0, 'T':0, 'G':0, 'C':0, 'N':0, '-':0} \
                    for pos, base in enumerate(path_seq, 0)}
    # make sub-graph for gap by disconnecting edges outside region
    subG = subset_graph(readG, source, target)
    # get positions of bases in each read based on position relative to source
    for node in subG:
        try:
            # path from source to node
            nodePath, nodePath_edges, nodePath_seq, nodePath_len = \
                    shortest_path(subG, source, node, trim = False)
            # path from node to target
            targetPath, targetPath_edges, targetPath_seq, targetPath_len = \
                    shortest_path(subG, node, target, trim = trim)
        except:
            continue
        # make sure path for node is same length as scaffold
        if nodePath_len + targetPath_len != path_len:
            continue
        for pos, base in enumerate(nodePath_seq, 0):
            try:
                bases[pos][base] += 1
            except:
                break
    # call consensus bases
    consensus = []
    for pos, base_counts in sorted(bases.items()):
        base_counts   = sorted([(count, base) for base, count in base_counts.items()])
        highest_count = base_counts[-1][0]
        # use path scaffold if coverage is 0
        if highest_count == 0:
            consensus.append(path_seq[pos])
        # choose randomly if there is a tie
        else:
            best_bases     = [base[1] for base in base_counts if base[0] == highest_count]
            consensus.append(Random.choice(best_bases))
    return ''.join(consensus)

def extend_sequence(readG, source, direction):
    """
    extend from sequence in specified direction
    """
    if direction == 'right':
        node_index = 0
        next_index = 1
    elif direction == 'left':
        node_index = 1
        next_index = 0
    edges = [(readG[edge[0]][edge[1]]['len_extension'], edge) for \
                edge in readG.edges() if edge[node_index] == source]
    path = []
    while len(edges) > 0:
        best_edge = sorted(edges)[-1][1]
        node = best_edge[next_index]
        path.append(best_edge)
        edges = [(readG[edge[0]][edge[1]]['len_extension'], edge) for \
                    edge in readG.edges() if edge[node_index] == node]
    if len(path) == 0:
        return ''
    target = path[-1][next_index]
    if direction == 'left':
        RreadG = reverse_graph(subset_graph(readG, target, source))
        return find_consensus(RreadG, source, target, trim = False)
    else:
        return find_consensus(readG, source, target, trim = False)

def fill_gap(A, gap, B, readG, pars):
    """
    fill gap in sequence using reads
    """
    readG = add_sequence(readG, A, 'left', pars)
    readG = add_sequence(readG, B, 'right', pars)
    # fill gap if path exists between A and B
    if nx.has_path(readG, A, B) is True:
        return find_consensus(readG, A, B, trim = True)
    # otherwise, extend inner ends of contig
    else:
        Aext = extend_sequence(readG, A, 'right')
        Bext = extend_sequence(readG, B, 'left')
        return ''.join([Aext, gap, Bext])

def trimNs(seq):
    """
    remove Ns from beginning and end of sequence
    """
    bases = ['A', 'T', 'G', 'C']
    if seq[0] in bases:
        start = 0
    else:
        start = min([seq.find(base) for base in bases])
    if seq[-1] in bases:
        end = len(seq)
    else:
        end = -1 * min([seq[::-1].find(base) for base in bases])
    return ''.join(list(seq)[start:end])

def pw_reads(reads, pars):
    """
    generator for comparing all reads with one another
    """
    min_length = pars['min_overlap'] + pars['min_extend']
    for A in reads:
        for B in reads:
            if A != B and len(A) >= min_length and len(B) >= min_length:
                yield (A, B)

def align_pair(pars):
    """
    align directed pair of reads A -> B
    by sliding read B over A to find the best
    extension
    """
    pair, pars = pars
    A, B = pair
    Astart = len(A) - pars['min_overlap']
    Bstart = 0
    Bend = pars['min_overlap'] - 1
    Blen = len(B)
    Lext = Blen - Bend
    alignments = []
    while Astart >= 0 and Lext >= pars['min_extend']:
        a = A[Astart:]
        b = B[Bstart:Bend]
        overlap   = Bend - Bstart
        extension = B[Bend:]
        # score alignment
        mm     = ldist(a, b)
        pident = float((overlap - mm)/overlap) * 100
        score  = overlap - (mm * 2)
        # save alignment if passing thresholds
        if \
                pident >= pars['min_pident']:
                    alignments.append((score, A, B, extension, Lext, overlap))
        elif \
                mm <= pars['max_mm']:
                    alignments.append((score, A, B, extension, Lext, overlap))
        # update variables for next alignment
        Astart -= 1
        Bend += 1
        Lext = Blen - Bend
    # return the best alignment, or FALSE if no alignments passed thresholds
    try:
        return sorted(alignments)[-1]
    except:
        return False

def read_graph(reads, pars):
    """
    creat directed network representating pairwise
    read comparisons
    """
    G = nx.DiGraph()
    pool = Pool(pars['threads'])
    num_compare = len(reads) * len(reads)
    print('# building overlap graph', file = sys.stderr)
    for pair in tqdm(pool.imap_unordered(align_pair, \
            [(pair, pars) for pair in pw_reads(reads, pars)]), total = num_compare):
        if pair is not False:
            score, A, B, extension, len_extension, overlap = pair
            G.add_edge(A, B, score = score, \
                    extension = extension, len_extension = len_extension, \
                    overlap = overlap, inv_overlap = (1/overlap))
    return G

def seq_extend(contigs, reads, no_extend, ignore_gaps, pars):
    """
    extend contig(s) and/or fill gap(s)
    """
    # graph of pair-wise read comparisons
    if pars['pickle_in'] is not False:
        pickle_in = open(pars['pickle_in'], 'rb')
        readG = Pickle.load(pickle_in)
    else:
        readG = read_graph(reads, pars)
    if pars['pickle_out'] is not False:
        pickle_out = open(pars['pickle_out'], 'wb')
        Pickle.dump(readG, pickle_out)
        pickle_out.close()
    # remove Ns from the ends of contigs
    for contig, seq in contigs.items():
        contigs[contig] = trimNs(seq)
    if ignore_gaps is not False:
        # parse contig fragments
        for contig, seq in contigs.items():
            fragments = parse_gaps(seq, pars)
            if fragments is False:
                continue
            num_fragments = len(fragments)
            filled = ''
            for fragment in range(0, num_fragments, 2):
                try:
                    A, gap, B = \
                            fragments[fragment], \
                            fragments[fragment + 1], \
                            fragments[fragment + 2]
                    filled += A
                    filled += fill_gap(A, gap, B, readG, pars)
                except:
                    break
            filled += B
            contigs[contig] = filled
    if no_extend is not False:
        for contig, seq in contigs.items():
            # add contig to readG
            readGC = add_sequence(readG, seq, 'both', pars)
            # extend left
            left  = extend_sequence(readG, seq, 'left')
            # extend right
            right = extend_sequence(readG, seq, 'right')
            seq = left + seq + right
            contigs[contig] = seq
    return contigs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = \
            '# extend contig(s) and/or fill gap(s)')
    parser.add_argument(\
            '-f', required = True, \
            help = 'contig(s) (fasta)')
    parser.add_argument(\
            '-r', required = True, \
            help = 'reads (fastq or fasta)')
    parser.add_argument(\
            '-k', default = False, \
            help = 'save overlap graph to pickle file (optional)')
    parser.add_argument(\
            '-g', default = False, \
            help = 'load overlap graph from pickle file (optional)')
    parser.add_argument(\
            '--no-extend', action = 'store_false', \
            help = 'do not extend')
    parser.add_argument(\
            '--ignore-gaps', action = 'store_false', \
            help = 'do not fill gaps')
    parser.add_argument(\
            '-p', default = float(98), type = float, \
            help = 'min. percent identity for alignments (default: 98)')
    parser.add_argument(\
            '-m', default = int(1), type = int, \
            help = 'max. number of mismatches allowed if identity < p (default: 1)')
    parser.add_argument(\
            '-o', default = int(20), type = int, \
            help = 'min. bp overlap between reads (default: 20)')
    parser.add_argument(\
            '-e', default = int(1), type = int, \
            help = 'min. bp for read extension (default: 1)')
    parser.add_argument(\
            '-b', default = int(10), type = int, \
            help = 'bp for error region buffer (default = 10)')
    parser.add_argument(\
            '-t', default = int(6), type = int, \
            help = 'number of threads (default: 6)')
    args = vars(parser.parse_args())
    seq_file, reads_file, no_extend, ignore_gaps = \
            args['f'], args['r'], args['no_extend'], args['ignore_gaps']
    contigs = {seq[0]:seq[1].upper() for seq in parse_fasta(open(seq_file))}
    reads = [read for read in parse_readsRC(open(reads_file))]
    pars = {'threads':args['t'], 'min_pident':args['p'], 'max_mm':args['m'], \
            'min_overlap':args['o'], 'min_extend':args['e'], 'buffer':args['b'], \
            'pickle_in':args['g'], 'pickle_out':args['k']}
    for contig, seq in seq_extend(contigs, reads, no_extend, ignore_gaps, pars).items():
        print(contig)
        print(seq)
