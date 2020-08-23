#!/usr/bin/env python

"""
For each allele cluster, this script extracts the shortest genomic regions harbouring all the alleles from
given genome assemblies (contigs) and cluster them using CD-HIT-EST.

Python versions 2 and 3 compatible.

Copyright (C) 2018 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
First version: 25 Nov 2018; the latest edition: 23 Aug 2020
"""

from __future__ import print_function
from __future__ import division
import os
import sys
sys.dont_write_bytecode = False  # Set True to disable generation of .pyc file
import argparse
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from collections import namedtuple, defaultdict


# Public data types
Path = namedtuple("Path", ["contig", "start", "end", "length"])  # the shortest path embedding all alleles in a contig; length = end - start + 1


def get_arguments():
    parser = argparse.ArgumentParser(description = "Searching allele clusters in contigs")
    
    # inputs
    parser.add_argument("--clusters", type = str, required = True, help = "A two-column tab-delimited file defining allele clusters")
    parser.add_argument("--allele_db", type = str, required = True, help = "A FASTA file for allele sequences")
    parser.add_argument("--assemblies", nargs = "+", type = str, required = True, help = "FASTA files for assembled contigs")
    parser.add_argument("--suffix", type = str, required = False, default = ".fasta", help = "Characters to be chopped off from the end of assembly filename in order to get a sample name")
    
    # outputs
    parser.add_argument("--outdir", type =str, required = False, default = ".", help = "Output directory")
    parser.add_argument("--prefix", type = str, required = False, default = "BLAST", help = "A prefix for output files")
    
    # environmental settings
    parser.add_argument("--makeblastdb", default = "makeblastdb", required = False, help = "Path to call makeblastdb")
    parser.add_argument("--blast", default = "blastn", required = False, help = "Path to call BLAST")
    parser.add_argument("--identity", type = int, required = False, default = 100, help = "The minimum nucleotide identity used for calling a hit.")
    parser.add_argument("--algorithm", default = "megablast", required = False, help = "BLAST algorithm (the -task argument)")
    parser.add_argument("--cdhit", default = "cd-hit-est", required = False, help = "Path to call CD-HIT-EST")
    parser.add_argument("--cdhit_args", default = "", required = False, help = "A string of arguments for CD-HIT-EST")
    
    # output control
    parser.add_argument("--skip", action = "store_true", required = False, help = "Flag it to skip existing output files")
    parser.add_argument("--clean", action = "store_true", required = False, help = "Flag it to delete original BLAST output files")
    
    return parser.parse_args()


def main():
    args = get_arguments()
    subdirs = make_output_dirs(args.outdir)
    cls = read_cluster_contents(args.clusters)
    queries = extract_allele_seqs(cluster_contents = cls, db = args.allele_db, outdir = subdirs["allele"],\
                                  prefix = args.prefix, skip = args.skip)
    assemblies = makeblastdb(prog = args.makeblastdb, assemblies = args.assemblies, db_dir = subdirs["database"],\
                             suffix = args.suffix, skip = args.skip)
    hits = blast(prog = args.blast, task = args.algorithm, queries = queries, dbs = assemblies,\
                 outdir = subdirs["blast"], prefix = args.prefix, skip = args.skip)
    exact_hits = concatenate_blast_output(hits = hits, prefix = args.prefix, outdir = subdirs["blast"], skip = args.skip,\
                                          clean = args.clean, min_id = args.identity)  # exact hits of alleles
    paths = find_shortest_paths(hits = exact_hits, cluster_contents = cls, outdir = subdirs["region"],\
                                prefix = args.prefix, skip = args.skip)  # assuming there is only a single copy per allele in a contig
    region_dirs = extract_cluster_seq(paths = paths, assemblies = assemblies, outdir = subdirs["region"], prefix = args.prefix,\
                                      skip = args.skip)
    cluster_regions(region_dirs = region_dirs, cdhit_path = args.cdhit, cdhit_args = args.cdhit_args,\
                    outdir = subdirs["cluster"], prefix = args.prefix, skip = args.skip)
    print("The cluster localisation process has been run through successfully.")


def make_output_dirs(outdir):
    # Create output directories
    print("Check output directories.")
    subdirs = {"allele" : os.path.join(outdir, "Allele"),
               "database" : os.path.join(outdir, "Database"),
               "blast" : os.path.join(outdir, "BLAST"),
               "region" : os.path.join(outdir, "Region"),
               "cluster" : os.path.join(outdir, "Cluster")}
    check_dir(outdir)
    for d in list(subdirs.values()):
        check_dir(d)
    
    return subdirs


def check_dir(d):
    # This is a subordinate function of make_output_dirs and extract_cluster_seq.
    if not os.path.exists(d):
        try:
            subprocess.check_call(["mkdir", d])  # stay compatible with older Python versions, although subprocess.run is available for Python 3.5+.
        except subprocess.CalledProcessError:
            print("* Error: the target directory cannot be created.")
            raise
    
    return


def read_cluster_contents(cluster_def):
    # Parameter cluster_def is the path to a tab-delimited file, which must not have a header line.
    print("Read cluster contents.")
    cls = {}
    with open(cluster_def) as f:
        lines = f.read().splitlines()
    for c in lines:
        cid, alleles = c.split("\t")
        alleles = alleles.split(",")
        cls[cid] = alleles
    print("* Totally %i clusters are defined." % len(cls))
    
    return cls


def extract_allele_seqs(cluster_contents, db, outdir, prefix, skip):
    # Create a FASTA file for each cluster to store their allele sequences
    print("Extract allele sequences of each cluster.")
    allele_db = list(SeqIO.parse(db, "fasta"))  # import the allele database
    queries = {}
    file_count = 0
    for cid, alleles in cluster_contents.items():
        if prefix != "":
            out_fasta = os.path.join(outdir, "__".join([prefix, cid, "alleles.fna"]))
        else:
            out_fasta = os.path.join(outdir, cid + "__alleles.fna")
        queries[cid] = out_fasta
        if not (os.path.exists(out_fasta) and skip):
            out_handle = open(out_fasta, "w")  # create a new file or override an existing file
            sum_alleles = len(alleles)
            allele_count = 0
            for rec in allele_db:
                if rec.id in alleles:
                    SeqIO.write(rec, out_handle, "fasta")
                    allele_count += 1
            out_handle.close()
            if allele_count < sum_alleles:
                print("* Warning: not all alleles of cluster %s are found in the allele database." % cid)
        file_count += 1
    print("* Totally %i FASTA files are generated (or already exist)." % file_count)
    
    return queries


def makeblastdb(prog, assemblies, db_dir, suffix, skip):
    print("Creating BLAST database from every assembly file.")
    Assembly = namedtuple("Assembly", ["fasta", "blast_db"])
    dbs = {}  # a dictionary of namedtuples
    n = 0
    for a in assemblies:
        a_base = os.path.basename(a)
        sample = a_base.rstrip(suffix)
        new_db = os.path.join(db_dir, sample)
        dbs[sample] = Assembly(fasta = a, blast_db = new_db)
        if not (os.path.exists(new_db + ".nhr") and os.path.exists(new_db + ".nin") and os.path.exists(new_db + ".nsq") and skip):
            cmd = [prog, "-in", a, "-dbtype", "nucl", "-out", new_db]
            try:
                subprocess.check_call(cmd)
            except subprocess.CalledProcessError:
                print("Runtime error: cannot make BLAST database for sample %s." % sample)
                raise
        n += 1
    print("* Totally %i BLAST databases are created (or already exist)." % n)
    
    return dbs


def blast(prog, task, queries, dbs, outdir, prefix, skip):
    print("Search allele clusters in contigs.")
    header = "qseqid sseqid sstart send pident qlen length gaps"
    blast_outputs = {}
    hit_num = 0
    nohit_num = 0
    n = 0
    
    if prefix != "":
        success_mark = os.path.join(outdir, "__".join([prefix, "blast.success"]))
    else:
        success_mark = os.path.join(outdir, "blast.success")
    
    for cid, query_fasta in queries.items():
        for sample, assembly in dbs.items():
            # Set output filename
            if prefix != "":
                output_file = os.path.join(outdir, "%s__%s__%s.txt" % (prefix, cid, sample))
            else:
                output_file = os.path.join(outdir, "%s__%s.txt" % (cid, sample))
                
            # Run BLAST
            if os.path.exists(output_file) and skip:  # only record the output filename
                if sample in list(blast_outputs.keys()):
                    blast_outputs[sample].append(output_file)
                else:
                    blast_outputs[sample] = [output_file]
            elif (not os.path.exists(success_mark)) or (not skip):  # Do nothing else.
                cmd = [prog, "-task", task, "-db", assembly.blast_db, "-query", query_fasta,\
                       "-outfmt", "6 " + header]
                process = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                out, err = process.communicate()
                out = out.decode()
                err = err.decode()
                n += 1  # job count + 1
                if len(err) > 0:
                    print("* Error: sequence search for cluster %s failed for the sample %s:" % (cid, sample),\
                          file = sys.stderr)
                    print(err, file = sys.stderr)
                    quit()
                if len(out) > 0:
                    hit_num += 1
                    with open(output_file, "w") as f:
                        f.write(out)
                    if sample in list(blast_outputs.keys()):
                        blast_outputs[sample].append(output_file)
                    else:
                        blast_outputs[sample] = [output_file]
                else:
                    nohit_num += 1  # Do not record the output filename as the file does not exist.
                    print("No hit of cluster %s in sample %s." % (cid, sample))
    subprocess.call(["touch", success_mark])
                
    print("* There are %i BLAST jobs launched, %i returned hits and %i did not have any hits." % (n, hit_num, nohit_num))
    print("* Header of every output file: " + header)
    
    return blast_outputs


def concatenate_blast_output(hits, prefix, outdir, skip, clean, min_id):
    """
    Concatenate BLAST outputs into a TSV file and delete the original outputs (*.txt).
    This function only keeps exact hits of alleles.
    """
    print("Parse and concatenate BLAST outputs.")
    Hit = namedtuple("Hit", ["allele", "contig", "start", "end", "identity", "allele_len", "hit_len", "gap_num"])
    hits_conc = defaultdict(dict)  # {cluster id : {sample : [Hit1, Hit2, ...]}}
    with_prefix = prefix != ""
    if with_prefix:
        f_out = os.path.join(outdir, prefix + "__BLAST.tsv")
    else:
        f_out = os.path.join(outdir, "BLAST.tsv")
    
    if os.path.exists(f_out) and skip:
        # Scenario 1: directly read outputs from a previous run.
        print("* Read %s for exact hits of alleles." % f_out)
        f = open(f_out, "rU")
        lines = f.read().splitlines()
        n = 0
        for h in lines:
            cid, sample, allele, contig, start, end, identity, allele_len, hit_len, gap_num = h.split("\t")
            hit = Hit(allele = allele, contig = contig, start = int(start), end = int(end), identity = float(identity),\
                      allele_len = int(allele_len), hit_len = int(hit_len), gap_num = int(gap_num))
            if cid not in list(hits_conc.keys()):
                hits_conc[cid] = {sample : [hit]}
            elif sample not in list(hits_conc[cid].keys()):
                hits_conc[cid][sample] = [hit]
            else:
                hits_conc[cid][sample].append(hit)
            n += 1
        print("* %i hits have been imported." % n)
        f.close()
    else:
        # Scenario 2: compile BLAST outputs
        for sample, out_files in hits.items():
            for out_f in out_files:
                if with_prefix:
                    _, cid, _ = out_f.split("__")
                else:
                    cid, _ = out_f.split("__")
            
                # Parse a single output file
                with open(out_f, "rU") as f_blast:
                    lines = f_blast.read().splitlines()
                for h in lines:
                    qseqid, sseqid, sstart, send, pident, qlen, length, gaps = h.split("\t")
                    hit = Hit(allele = qseqid, contig = sseqid, start = int(sstart), end = int(send),\
                              identity = float(pident), allele_len = int(qlen), hit_len = int(length),\
                              gap_num = int(gaps))  # a hit of an allele in the current cluster in the current sample
                    if hit.identity >= min_id and hit.gap_num == 0 and hit.allele_len == hit.hit_len:
                        # Add this hit into the dictionary; else, discard.
                        if cid not in list(hits_conc.keys()):
                            hits_conc[cid] = {sample : [hit]}
                        elif sample not in list(hits_conc[cid].keys()):
                            hits_conc[cid][sample] = [hit]
                        else:
                            hits_conc[cid][sample].append(hit)
        
        # Write concatenated BLAST hits into a file
        print("Print exact hits of alleles.")
        f = open(f_out, "w")
        for cid, hits_dict in hits_conc.items():
            for sample, hit_list in hits_dict.items():
                for h in hit_list:
                    f.write("\t".join([cid, sample, h.allele, h.contig, str(h.start), str(h.end), str(h.identity),\
                                       str(h.allele_len), str(h.hit_len), str(h.gap_num)]) + "\n")
        f.close()
                
    # Finally, delete all *.txt files
    if clean:
        print("Delete original BLAST outputs.")
        subprocess.call(["rm", os.path.join(outdir, "*.txt")])
    
    return hits_conc

    
def find_shortest_paths(hits, cluster_contents, outdir, prefix, skip):
    # Assume all hits are exact to corresponding alleles.
    print("Find the shortest path embedding all alleles per cluster in each contig.")
    Path = namedtuple("Path", ["contig", "start", "end", "length"])  # the final output of this function
    paths = defaultdict(dict)  # {cid : {sample : [Path(contig, start, end), ...]}}
    if prefix != "":
        out_file = os.path.join(outdir, prefix + "__paths.tsv")
    else:
        out_file = os.path.join(outdir, "paths.tsv")
    
    if os.path.exists(out_file) and skip:
        print("* Read previously determined shortest paths.")
        f = open(out_file, "rU")
        lines = f.read().splitlines()
        for p in lines:
            cid, sample, contig, start, end, length = p.split("\t")
            if cid not in list(paths.keys()):
                paths[cid] = {sample : [Path(contig = contig, start = int(start), end = int(end), length = int(length))]}
            elif sample not in list(paths[cid].keys()):
                paths[cid][sample] = [Path(contig = contig, start = int(start), end = int(end), length = int(length))]
            else:
                paths[cid][sample].append(Path(contig = contig, start = int(start), end = int(end), length = int(length)))
        f.close()
    else:  # for new runs
        print("* Determine and save shortest paths to file %s." % out_file)
        f = open(out_file, "w")
        
        # Run path determination for each cluster and each sample.
        for cid, hits_dict in hits.items():
            cluster_size = len(cluster_contents[cid])  # number of alleles per cluster
            for sample, hit_list in hits_dict.items():
                if len(hit_list) >= cluster_size:  # The number of hits must not be smaller than the number of alleles of the current cluster.
                    contigs = get_hits_in_contigs(hit_list)  # a two-dimensional dictionary for the current cluster and sample
                    for contig, allele_dict in contigs.items():  # alleles: {allele : Allele(start = [...], end = [...])}
                        """
                        Below, I show the only condition that starts the procedure of looking for the shortest path per contig.
                        len(allele_dict.keys()) < cluster_size when the current contig does not harbour all the alleles.
                        get_shortest_path_per_contig is the most sophisticated and perhaps the slowest function of this script.
                        """
                        if len(list(allele_dict.keys())) == cluster_size:  # Cases where not all alleles are found in a contig are skipped.
                            new_path = get_shortest_path_per_contig(contig_name = contig, alleles = allele_dict,\
                                                                    cluster_size = cluster_size, cid = cid,\
                                                                    sample = sample)
                            
                            # Store the identified path into a two-dimensional dictionary.
                            if cid not in list(paths.keys()):
                                paths[cid] = {sample : [new_path]}
                            elif sample not in list(paths[cid].keys()):
                                paths[cid][sample] = [new_path]
                            else:
                                paths[cid][sample].append(new_path)
                            f.write("\t".join([cid, sample, new_path.contig, str(new_path.start), str(new_path.end),\
                                               str(new_path.length)]) + "\n")
        f.close()
    
    return paths


def get_hits_in_contigs(hit_list):
    """
    This is the first subordinate function of get_shortest_paths. The data structure of contigs can be displayed as:
    contig name -- level-1 key
        allele  -- level-2 key
            [(start1, end1)
            (start2, end2)
            ...] -- a list of namedtuples
    """
    Allele_coords = namedtuple("Allele_coords", ["start", "end"])  # Attributes "start" and "end" are lists.
    contigs = defaultdict(dict)  # {contig name : {allele : [Allele_coords(start, end), Allele_coords(start, end), ...]}}
    for h in hit_list:
        if h.contig not in list(contigs.keys()):
            contigs[h.contig] = {h.allele : [Allele_coords(start = h.start, end = h.end)]}
        elif h.allele not in list(contigs[h.contig].keys()):
            contigs[h.contig][h.allele] = [Allele_coords(start = h.start, end = h.end)]
        else:  # Another copy of a known allele is found.
            contigs[h.contig][h.allele].append(Allele_coords(start = h.start, end = h.end))  # append a new tuple into the list
    
    return contigs


def get_shortest_path_per_contig(contig_name, alleles, cluster_size, cid, sample):
    """
    This is the second subordinate function of get_shortest_paths. It returns a Path object.
    alleles: {allele : [Allele(start, end), Allele(start, end), ...]}. This function addresses the challenge in determining the
    shortest path when one or more alleles having multiple exact copies in the same contig.
    """
    cn = count_copy_numbers(alleles)  # returns a namedtuple of copy numbers for alleles in the current contig
    only_single_copies = check_singularity(cn, contig_name, cid, sample)  # returns a boolean value showing whether there is only a single copy per allele
    if only_single_copies:
        p = find_shortest_path_single_copy(alleles, contig_name)
    else:
        p = find_shortest_path_multi_copy(cn, alleles, contig_name)
        
    return p


def count_copy_numbers(allele_dict):
    # A subordinate function of get_shortest_path_per_contig, and it counts the number of copies per allele.
    CN = namedtuple("CN", ["alleles", "copies"])
    cn = CN(alleles = [], copies = [])
    for allele_name, coordinate_list in allele_dict.items():
        cn.alleles.append(allele_name)  # [allele1, allele2, allele3, ...]
        cn.copies.append(len(coordinate_list))  # [2, 1, 3, ...]
    
    return cn


def check_singularity(cn, contig_name, cid, sample):
    # Determine whether there is only a single copy per allele
    is_singular = True
    for i in list(range(0, len(cn.alleles))):
        n = cn.copies[i] 
        if n > 1:
            print("* Notice: allele %s of the cluster %s in the contig %s of sample %s has %i copies." %\
                  (cn.alleles[i], cid, contig_name, sample, n))
            is_singular = False  # change the Boolean value forever
    
    return is_singular


def find_shortest_path_single_copy(alleles, contig_name):
    # A subordinate function of get_shortest_path_per_contig    
    low_coords = []
    high_coords = []
    for coordinate_list in list(alleles.values()):
        c = coordinate_list[0]  # There must be only one element in the list.
        s = c.start
        e = c.end
        if s > e:  # in complementary orientaion
            low_coords.append(e)
            high_coords.append(s)
        else:
            low_coords.append(s)
            high_coords.append(e)
    low = min(low_coords)
    high = max(high_coords)
    p = Path(contig = contig_name, start = low, end = high, length = high - low + 1)
        
    return p


def find_shortest_path_multi_copy(cn, alleles, contig_name):
    perms = permutation_generator(cn.copies)
    region_lengths = []  # to store lengths of all possible regions so that we can find out the shortest one
    regions = []  # a list of Path objects for all possible regions
    p = None  # initial value
    for pm in perms:  # pm is a list of indices specifying which copy of each allele will be used for calculating the region length.
        regions.append(calc_region_length(indices = pm, allele_names = cn.alleles, allele_dict = alleles))
    
    # Find out the shortest length
    n = len(regions)
    for i in list(range(0, n)):
        region_lengths.append(regions[i].length)
    len_min = min(region_lengths)
    
    # Determine which region has the shortest length
    for i in list(range(0, n)):
        r = regions[i]
        if r.length == len_min:  # only take the first region when there is a tie
            p = Path(contig = contig_name, start = r.start, end = r.end, length = len_min)
            break

    return p


def permutation_generator(ns):
    """
    An elegant generator of permutations of m steps, for which the i-th step has n_i choices. Therefore,
    this function returns a list of n_1 * n_2 * ... * n_m elements for all permutations. Assuming there
    are alleles [a1, a2, a3, ...] in an array, then perms = [i1, i2, i3, ...] for indices of copies of
    these alleles.
    """
    perms = []
    if len(ns) > 1:
        perms_stack = permutation_generator(ns[1 : ])  # permutations from the stack of steps
        for i in list(range(0, ns[0])):
            for j in perms_stack:  # j is a list in perms_prev
                perms.append([i] + j)  # [[1, 2, 3, ...], [2, 5, 3, ...], ...], where each element list represents a permutation
    else:  # The function reaches the bottom of the stack, when there is only one element in the list ns.
        for i in list(range(0, ns[0])):
            perms.append([i])  # produces [[0], [1], [2], ..., [ns[0]]]
    
    return perms


def calc_region_length(indices, allele_names, allele_dict):
    """
    Calculate a region width given a list of copy indices of alleles
    alleles: a dictionary; copy_numbers: a list of CN objects; indices: a list of indices following the same order
    as the argument "alleles". In theory, find_shortest_path_single_copy can be replaced by this function. However,
    the former is faster than the latter.
    """
    low_coords = []
    high_coords = []
    for i in list(range(0, len(indices))):  # len(indices): the number of alleles in the current contig
        j = indices[i]  # the (j + 1)-th copy of the (i + 1)-th allele in copy_numbers.alleles
        a = allele_names[i]
        pos = allele_dict[a][j]  # take the (j + 1)-th copy of this allele
        s = pos.start
        e = pos.end
        if s > e:  # in complementary orientaion
            low_coords.append(e)
            high_coords.append(s)
        else:
            low_coords.append(s)
            high_coords.append(e)
    low = min(low_coords)
    high = max(high_coords)
    region = Path(contig = None, start = low, end = high, length = high - low + 1)

    return region


def extract_cluster_seq(paths, assemblies, outdir, prefix, skip):
    # Create output directories
    cids = list(paths.keys())
    dirs = {}  # {cid : directory path}
    for cid in cids:
        d = os.path.join(outdir, cid)  # parental dir/Region/1, ...
        check_dir(d)
        dirs[cid] = d
        
    # Extract sequences of each cluster
    for cid, samples in paths.items():
        for sample, paths in samples.items():
            if prefix != "":
                f_out = os.path.join(dirs[cid], "%s__%s.fna" % (prefix, sample))
            else:
                f_out = os.path.join(dirs[cid], sample + ".fna")
            
            if not (os.path.exists(f_out) and skip):
                fasta = assemblies[sample].fasta  # FASTA file of the original genome assembly
                with open(fasta, "rU") as fasta_handle:
                    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))  # key: sequence ID; value: SeqRecord object.
                f_out_handle = open(f_out, "w")
                for p in paths:
                    contig = fasta_dict[p.contig]  # a SeqRecord object
                    contig_seq = str(contig.seq)
                    contig_seq = contig_seq[(p.start - 1) : p.end]
                    """
                    Do not change attributes of contig directly, as it will affect fasta_dict[p.contig]. It is better to create
                    a new record from contig. Further, there is no need to add an index to the sequence ID by far, as currently,
                    only one region is reported per contig by this script. It is necessary to append the sample name to each
                    sequence for the next step -- sequence clustering.
                    """
                    new_record = SeqRecord(Seq(contig_seq, generic_dna),\
                                           id = contig.id + "|" + sample,\
                                           description = ",".join([str(p.start), str(p.end), str(p.length)]))
                    SeqIO.write(new_record, f_out_handle, "fasta")
                f_out_handle.close()
    
    return dirs


def cluster_regions(region_dirs, cdhit_path, cdhit_args, outdir, prefix, skip):
    if cdhit_path != "":
        print("Perform clustering to extracted region sequences.")
        for cid, seq_dir in list(region_dirs.items()):
            print("* Concatenate FASTA files into a single temporary file for cluster %s." % cid)
            
            if prefix != "":
                seq_files = os.path.join(seq_dir, "%s__*.fna" % prefix)  # Region/[cid]/[prefix]__*.fna
                tmp_fasta = os.path.join(seq_dir, "%s__fasta.tmp" % prefix)
                out_file = os.path.join(outdir, "%s__%s__clusters.fna" % (cid, prefix))
                success_flag = os.path.join(outdir, "%s__%s__clustering.success" % (cid, prefix))
            else:
                seq_files = os.path.join(seq_dir, "*.fna")  # Region/[cid]/*.fna
                tmp_fasta = os.path.join(seq_dir, "fasta.tmp")
                out_file = os.path.join(outdir, cid + "__clusters.fna")
                success_flag = os.path.join(outdir, cid + "__clustering.success")
            
            if os.path.exists(success_flag) and skip:
                print("* Skip clustering for the cluster %s." % cid)
            else:
                proc = subprocess.Popen("cat %s > %s" % (seq_files, tmp_fasta), shell = True)  # Sample names have been attached to sequences previously.
                proc.wait()  # wait until the command finishes
                cmd = [cdhit_path, "-i", tmp_fasta, "-o", out_file] + cdhit_args.split(" ")
                try:
                    subprocess.check_call(cmd)
                except subprocess.CalledProcessError:
                    print("Command line '%s' failed." % " ".join(cmd))
                    raise
                subprocess.call(["rm", tmp_fasta])
                subprocess.call(["touch", success_flag])
    else:
        print("Sequence cluster is skipped as CD-HIT-EST is not specified.")
    
    return


if __name__ == "__main__":
    main()
