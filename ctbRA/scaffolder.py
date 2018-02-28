#!/usr/bin/env python3

"""
script for scaffolding contigs based on paired read coverage
"""

import os
import sys
import networkx as nx
from operator import itemgetter

import ctbRA.re_assemble_errors as curate
from ctbBio.search import search as search
from ctbBio.fasta import iterate_fasta as parse_fasta
from ctbBio.rc import reverse_complement as rc

def combine_fastas(fastas, out):
    """
    print multiple fastas to single file
    """
    out = open(out, 'w')
    for f in fastas:
        for seq in list(f.values()):
            for s in seq:
                print('\n'.join([seq[0], seq[2]]), file=out)
    out.close()
    return out.name

def assemble_extensions(contigs, mapping, mismatches, collection_mismatches, threads, prefix):
    """
    find reads that map to ends of scaffolds and re-assemble them
    """
    pairs, header = curate.parse_mapping(mapping[0], ends = 1000, scaffolds = contigs[1])
    s2pos = {id: [0, info[1]] for id, info in list(contigs[1].items())} # s2pos[contig] = [0, length of contig]
    reads = curate.collect_reads(\
            pairs, contigs[0], contigs[1], collection_mismatches, prefix, s2errors = s2pos, window = 1000, combine_windows = True)
    extensions = curate.re_assemble_fragments(\
            reads, prefix, threads, cov_thresh = 1, mismatches = mismatches, \
            multiple = True, scaffolding = True, assembler = 'velvet', min_contig = '200')
    # put extensions into same format as contigs
    ext = {}
    for contig, assemblies in list(extensions.items()):
        for pos, assembly in list(assemblies.items()):
            if assembly is False:
                continue
            for seq in assembly:
                info, seq = seq
                length, header = info
                id = '%s-%s-%s' % (contig, pos, header.split('>', 1)[1])
                ext[id] = ['>%s' % (id), length, seq]
    return ext

def check_rc(seq, coords):
    """
    reverse complement sequence, if necessary
    """
    if coords[0] > coords[1]:
        return rc(['', seq])[1]
    return seq

def trim_for_scaffolding(seqs):
    """
    trim contigs so they can be concatenated
    - check to see if sequence needs to be reverse complemented
    """
    first, second = seqs[0][1], seqs[1][1]
    fid, fseq, flen, faln = first
    fstart, fstop = min(faln), max(faln)
    ftrimmed = check_rc(fseq[0:fstop], faln)
    sid, sseq, slen, saln = second
    sstart, sstop = min(saln), max(saln)
    strimmed = check_rc(sseq[sstop:], saln)
    return ftrimmed, strimmed

def cat(query, target, pident, buffer = 25, min_extension = 200):
    """
    query/target = [id, seq, length, [start, stop]]
    - alignment at beginning or end of sequence?
    - update pident in case alignment does not extend to end of the read
    - combine sequences, reverse compliment added fragment if necessary
    """
    seqs = []
    for contig in query, target:
        id, seq, length, aligned = contig
        start, stop = min(aligned), max(aligned)
        left = start - 0
        right = length - stop
        if min([left, right]) > buffer:
            return False
        if left > right: # alignment is at end of sequence 
            pos = 1
        else:
            pos = 0 # alignment is at beginning of sequence
        seqs.append([pos, contig])
    sorted(seqs, key = itemgetter(0))
    # make sure direction is correct
    if seqs[0][0] == seqs[1][0]:
        return False
    id = '%s-cat-%s' % (seqs[0][1][0], seqs[1][1][0])
    first, second = trim_for_scaffolding(seqs)
    scaffold = [first, second]
    join = len(scaffold[0])
    added = min([len(scaffold[1]), len(scaffold[0])])
    if added < min_extension:
        return False
    return [seqs[0][1][0], seqs[1][1][0], id, join, added, scaffold]

def contig_join(join):
    """
    get last 1kb of first contig and first 1kb of second contig
    """
    a = join[1][0][-1000:]
    b = join[1][1][:1000]
    header = '>%s' % (join[0])
    length = len(a + b)
    return [header, length, '\n'.join([a, b])]

def network_contigs(contigs, extensions, blast, overlap, prefix):
    """
    create network of contig connections based on overlap and paired reads
    """
    g = nx.DiGraph()
#    g = nx.Graph()
    id2scaffold = {}
    nodes = dict(list(contigs.items()) + list(extensions.items()))
    for hit in open(blast):
        hit = hit.strip().split()
        q, t = hit[0], hit[1]
        if q == t:
            continue
        if q not in nodes or t not in nodes:
            continue
        qstart, qstop = int(hit[6]), int(hit[7]) 
        tstart, tstop = int(hit[8]), int(hit[9])
        qlength, tlength = nodes[q][1], nodes[t][1]
        qseq, tseq = nodes[q][2], nodes[t][2]
        pident, bit = float(hit[2]), float(hit[11])
        if pident < overlap:
            continue
        scaffold = cat([q, qseq, qlength, [qstart, qstop]], [t, tseq, tlength, [tstart, tstop]], pident)
        if scaffold is not False:
           a, b, dir_id, join, added, scaffold = scaffold
           id = ''.join(sorted([a, b]))# id = uniq to pair, dir_id reflects directionality
           if id in id2scaffold:
               continue
           id2scaffold[id] = [dir_id, scaffold]
           g.add_edge(a, b, dir_id = dir_id, id = id, pident = pident, \
                   bit = bit, join_pos = join, bases_added = added)
    joins = {i[0]: contig_join(i) for i in list(id2scaffold.values())}
    joins_out = open('%s/contig_joins.fa' % (prefix), 'w')
    for s in list(joins.values()):
        print('\n'.join([s[0], s[2]]), file=joins_out)
    joins_out.close()
    nx.write_edgelist(g, '%s/0.all-connections.network' % (prefix), \
            data = ['dir_id', 'id', 'join_pos', 'bases_added', 'pident', 'bit'])
    return g, [joins_out.name, joins], id2scaffold

def no_errors(scaffold, errors, end = 200):
    """
    check to see if there is an error in the scaffold - ignore ends
    """
    filtered = []
    for e in errors:
        if e == 0:
            continue
        if e >= scaffold[1] - end:
            continue
        filtered.append(e)
    if len(filtered) == 0:
        return True
    return False

def filter_graph_mapping(joins, graph, pr, threads, prefix, pr_split, mismatches):
    """
    check joind contigs with paired read mapping, remove edges in graph for unsupported joins
    """
    mapping, pr_split = curate.map_reads(joins[0], joins[1], pr, threads, pr_split = pr_split, multiple = True)
    pairs, header = curate.parse_mapping(mapping)
    s2c, s2errors = curate.id_errors(\
            pairs, header, joins[0], joins[1], cov_thresh = 1, mismatches = mismatches, allow_orphan = True)
    for node, connections in graph.adjacency():
        for c, info in list(connections.items()):
            if 'dir_id' not in info:
                continue
            dir_id = info['dir_id']
            # remove connection if scaffolding does not pass error check
            if no_errors(joins[1][dir_id], s2errors[dir_id]) is False:
                graph.remove_edge(node, c)
    nx.write_edgelist(graph, '%s/1.supported-connections.network' % (prefix), \
            data = ['dir_id', 'id', 'join_pos', 'bases_added', 'pident', 'bit'])
    return graph

def follow_path(c, graph, added):
    """
    determine the length of the path from node c until next junction
    """
    connections = len(list(graph[c].keys()))
    p = c
    while connections == 1:
        n = graph.successors(p)
        connections = len(n)
        if connections != 1:
            break
        added += graph[p][n[0]]['bases_added']
        p = n[0]
    return [added, c]

def filter_graph_best(graph, prefix):
    """
    find best path through graph (remove additional edges)
    - based on furthest distance (nucleotides) before arriving at another junction
    """
    filtered = nx.DiGraph()
    # filter in forward direction
    for node, connections in graph.adjacency():
        successors = graph.successors(node)
        # how far does each successor extend through the graph?
        length_connections = [follow_path(c, graph, graph[node][c]['bases_added']) for c in successors]
        ranked = sorted(length_connections, key = itemgetter(0), reverse = True)
        if len(ranked) == 0:
            continue
        best = ranked[0][1]
        filtered.add_edge(node, best, **graph[node][best])
    # remove multiple nodes
    for node, connections in filtered.adjacency():
        if filtered.in_degree(node) <= 1:
            continue
        prev = [[filtered[i][node]['bases_added'], i] for i in filtered.predecessors(node)]
        ranked = sorted(prev, key = itemgetter(0), reverse = True)
        for i in ranked[1:]:
            i = i[1]
            filtered.remove_edge(i, node)
    nx.write_edgelist(filtered, '%s/2.best-connections.network' % (prefix), \
            data = ['dir_id', 'id', 'join_pos', 'bases_added', 'pident', 'bit'])
#    nx.draw_graphviz(filtered, node_size = 21, width = 0.1, font_size = 0.9)
#    plt.savefig('%s/2.best-connections.network.pdf' % (prefix))
    return filtered

def order_group(graph, group):
    """
    order elements in group based on directed graph
    """
    ordered = []
    first = [n for n in group if graph.in_degree(n) == 0][0]
    prev = first
    while prev is not False:
        successor = [i for i in graph.successors(prev)]
        if len(successor) == 0:
            prev = False
            break
        successor = successor[0]
        ordered.append([prev, successor, graph[prev][successor]['id']])
        prev = successor
    return ordered

def check_contig(group, contigs):
    """
    make sure at least one contig is in the group
    """
    for g in group:
        if g in contigs:
            return True
    return False

def scaffold_contigs(contigs, id2scaffolds, graph, prefix):
    """
    combine contigs into scaffolds based on error-checked overlap graph
    """
    fasta = open('%s/scaffolds.fa' % (prefix), 'w')
    scaffolds = {}
    scaffolded = []
    for group in nx.connected_components(graph.to_undirected()):
        if len(group) == 1:
            continue
        if check_contig(group, contigs) is False:
            continue
        group = order_group(graph, group)
        for g in group:
            scaffolded.append(g[0])
            scaffolded.append(g[1])
        id = '_'.join([g[2] for g in group])
        scaffolds[id] = []
        for i in id2scaffolds[group[0][2]][1]:
            scaffolds[id].append(i)
        for g in group[1:]:
            seq = id2scaffolds[g[2]][1][1]
            scaffolds[id].append(seq)
    for contig, info in list(contigs.items()):
        if contig not in scaffolded:
            scaffolds[contig] = info[2]
    for id, seq in list(scaffolds.items()):
        print('>%s' % (id), file=fasta)
        print(''.join(seq), file=fasta)
    fasta.close()
    return [fasta.name, scaffolds]

def scaffolder(contigs, pr, pr_split, prefix, threads, mismatches, collection_mismatches = 2, overlap = 0.90):
    """
    scaffold contigs based on overlap, paired read support, and paired read confirmation of scaffold
    1) map reads to contigs (best mapping for each read)
    2) re-assemble reads that map to ends of conitgs in order to find extensions
    3) blast contigs and extensions against one another and find where ends overlap
    4) create network of contig connections (require both paired read connections and overlap)
    5) check all possible contig joins with paired read mapping
    6) break edges in network if joinn is not supported by paired read mapping
    7) traverse network to build scaffolds
    """
    # prepare scaffolder directory
    os.system('mkdir -p %s' % (prefix))
    contigs = curate.format_assembly(contigs, prefix)
    # save contigs to dictionary: contigs = [file_path, contigs[id] = [header, length, sequence]]
    contigs = [contigs, {i[0].split('>')[1]: [i[0], len(i[1]), i[1]] for i in parse_fasta(contigs) if i != [[], []]}]
    # map reads to contigs
    mapping = curate.map_reads(\
            contigs[0], contigs[1], pr, threads, multiple = False, pr_split = pr_split)
    # find possible contig extensions
    extensions = assemble_extensions(contigs, mapping, mismatches, collection_mismatches, threads, prefix)
    combined = combine_fastas([contigs[1], extensions], '%s/combined.fa' % (prefix))
    # find overlap between contigs and extensions
    blast_out = search(combined, combined, method = 'blast', \
            max_hits = 100, threads = threads, prefix = '%s/' % (prefix))
    # create network of contig connections
    graph, joins, id2scaffolds = network_contigs(contigs[1], extensions, blast_out, overlap, prefix)
    # check joined contigs with paired read mapping
    graph = filter_graph_mapping(joins, graph, pr, threads, prefix, pr_split, mismatches)
    # find best path through graph
    graph = filter_graph_best(graph, prefix)
    # scaffold contigs based on graph
    scaffolds = scaffold_contigs(contigs[1], id2scaffolds, graph, prefix)

if __name__ == '__main__':
    if len(sys.argv) != 8:
        print('specify threads, # mismatches tolerated, % id of overlaps, fasta for scaffolding, paired reads (fastq), pair 1 (fastq), and pair 2 (fastq)', file=sys.stderr)
        exit()
    for i, v in enumerate(sys.argv[3:], 3):
        if v == 'False' or v == 'FALSE' or v == 'false':
            sys.argv[i] = False
    threads, mm, overlap, assembly, pr, pr1, pr2 = sys.argv[1:]
    if pr1 is False or pr2 is False:
        pr_split = False
    else:
        pr_split = [pr1, pr2]
    prefix = '%s.scaffolded' % (assembly.rsplit('.', 1)[0])
    scaffolder(assembly, pr, pr_split, prefix, int(threads), mismatches = int(mm), overlap = float(overlap))
