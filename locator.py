#!/usr/bin/env python

"""
For each allele cluster, this script extracts the shortest genomic regions harbouring all the alleles from
given genome assemblies (contigs) and cluster them using CD-HIT-EST.

Python versions 2.7 and 3 compatible.

Copyright (C) 2018 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
First version: 25 Nov 2018; the latest edition: 26 Nov 2018
"""

from __future__ import print_function
from __future__ import division
import os
import sys
sys.dont_write_bytecode = True
import argparse
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from collections import namedtuple, defaultdict


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
                                          clean = args.clean)  # exact hits of alleles
    paths = find_shortest_paths(hits = exact_hits, cluster_contents = cls, outdir = subdirs["region"],\
                                prefix = args.prefix, skip = args.skip)  # assuming there is only a single copy per allele in a contig
    region_dirs = extract_cluster_seq(paths, assemblies, subdirs["region"], args.prefix)
    

def extract_cluster_seq(paths, assemblies, outdir, prefix):
    # Create output directories
    cids = list(paths.keys())
    dirs = {}
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
                only one region is reported per contig by this script.
                """
                new_record = SeqRecord(Seq(contig_seq, generic_dna),\
                                       id = contig.id + "|" + sample,\
                                       description = ",".join([str(p.start), str(p.end), str(p.length)]))
                SeqIO.write(new_record, f_out_handle, "fasta")
            f_out_handle.close()
    
    return dirs
    
    
def find_shortest_paths(hits, cluster_contents, outdir, prefix, skip):
    # Assume all hits are exact to corresponding alleles.
    print("Find the shortest path embedding all alleles per cluster in each contig.")
    Path = namedtuple("Path", ["contig", "start", "end", "length"])
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
    else:
        print("* Determine and save shortest paths to file %s." % out_file)
        f = open(out_file, "w")
        for cid, hits_dict in hits.items():
            cluster_size = len(cluster_contents[cid])  # number of alleles per cluster
            for sample, hit_list in hits_dict.items():
                if len(hit_list) >= cluster_size:  # The number of hits must not be smaller than the number of alleles of the current cluster.
                    contigs = get_hits_in_contigs(hit_list)  # for the current cluster and sample
                    for contig, alleles in contigs.items():  # alleles: {allele : Allele(start = [...], end = [...])}
                        if len(list(alleles.keys())) == cluster_size:  # len(alleles.keys()) < cluster_size when the current contig does not harbour all the alleles
                            new_path = get_shortest_path_per_contig(contig, alleles, cid, sample)
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
    # This is the first subordinate function of get_shortest_paths.
    Allele = namedtuple("Allele", ["start", "end"])  # Attributes "start" and "end" are lists.
    contigs = defaultdict(dict)  # {contig name : {allele : Allele(start = [...], end = [...])}}
    for h in hit_list:
        if h.contig not in list(contigs.keys()):
            contigs[h.contig] = {h.allele : Allele(start = [h.start], end = [h.end])}
        elif h.allele not in list(contigs[h.contig].keys()):
            contigs[h.contig][h.allele] = Allele(start = [h.start], end = [h.end])
        else:  # Another copy of a known allele is found.
            contigs[h.contig][h.allele].start.append(h.start)
            contigs[h.contig][h.allele].end.append(h.end)
    
    return contigs


def get_shortest_path_per_contig(contig_name, alleles, cid, sample):
    """
    This is the second subordinate function of get_shortest_paths. It returns a Path object.
    alleles: {allele : Allele(start = [...], end = [...])}.
    Currently, this function can only handle the scenario where there is only a single copy per
    allele in a contig.
    """
    Path = namedtuple("Path", ["contig", "start", "end", "length"])  # the shortest path embedding all alleles in a contig; length = end - start + 1
    low_coords = []
    high_coords = []
    for allele_name, coords in alleles.items():
        coords_len = len(coords.start)
        if coords_len > 1:
            print("* Warning: allele %s of cluster %s has more than one copy in contig %s of sample %s. " % (allele_name,\
            cid, contig_name, sample))
        s = coords.start[0]
        e = coords.end[0]
        if s > e:  # in complementary orientaion
            low_coords.append(e)
            high_coords.append(s)
        else:
            low_coords.append(s)
            high_coords.append(e)
        low = min(low_coords)
        high = max(high_coords)
        
    return Path(contig = contig_name, start = low, end = high, length = high - low + 1)


def concatenate_blast_output(hits, prefix, outdir, skip, clean):
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
                    if hit.identity == 100 and hit.gap_num == 0 and hit.allele_len == hit.hit_len:
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
    # cluster_def is the path to a tab-delimited file, which must not have a header line.
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


if __name__ == "__main__":
    main()
