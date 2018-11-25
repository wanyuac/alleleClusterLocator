#!/usr/bin/env python

"""
For each allele cluster, this script extracts the shortest genomic regions harbouring all the alleles from
given genome assemblies (contigs) and cluster them using CD-HIT-EST.

Python versions 2.7 and 3 compatible.

Copyright (C) 2018 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
First version: 25 Nov 2018; the latest edition: 25 Nov 2018
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
from collections import namedtuple


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
    

def blast(prog, task, queries, dbs, outdir, prefix, skip):
    print("Search allele clusters in contigs.")
    header = "qseqid sseqid sstart send pident qlen length bitscore"
    blast_outputs = {}
    n = 0
    for cid, query_fasta in queries.items():
        for sample, assembly in dbs.items():
            if prefix != "":
                output_file = os.path.join(outdir, "%s__%s__%s.tsv" % (prefix, cid, sample))
            else:
                output_file = os.path.join(outdir, "%s__%s.tsv" % (cid, sample))
            if not (os.path.exists(output_file) and skip):
                cmd = [prog, "-task", task, "-db", assembly.blast_db, "-query", query_fasta,\
                       "-outfmt", "6 " + header]
                process = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                out, err = process.communicate()
                out = out.decode()
                err = err.decode()
                if len(err) > 0:
                    print("* Error: sequence search for cluster %s failed for the sample %s:" % (cid, sample),\
                          file = sys.stderr)
                    print(err, file = sys.stderr)
                    quit()
                with open(output_file, "w") as f:
                    f.write(out)
                n += 1
            blast_outputs[sample] = output_file
    print("* %i BLAST jobs are implemented." % n)
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
            alleles = alleles.split(",")
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
    if not os.path.exists(outdir):
        try:
            subprocess.check_call(["mkdir", outdir])  # stay compatible with older Python versions, although subprocess.run is available for Python 3.5+.
        except subprocess.CalledProcessError:
            print("* Error: the parental output directory cannot be created.")
            raise
    for d in list(subdirs.values()):
        if not os.path.exists(d):
            subprocess.call(["mkdir", d])
    
    return subdirs


def read_cluster_contents(cluster_def):
    # cluster_def is the path to a tab-delimited file, which must not have a header line.
    print("Read cluster contents.")
    cls = {}
    with open(cluster_def) as f:
        lines = f.read().splitlines()
    for c in lines:
        cid, alleles = c.split("\t")
        cls[cid] = alleles
    print("* Totally %i clusters are defined." % len(cls))
    
    return cls


if __name__ == "__main__":
    main()
