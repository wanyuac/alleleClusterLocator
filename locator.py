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
    parser.add_argument("--blast", default = "blastn", required = False, help = "Path to call BLAST")
    parser.add_argument("--algorithm", default = "megablast", required = False, help = "BLAST algorithm (the -task argument)")
    parser.add_argument("--cdhit", default = "cd-hit-est", required = False, help = "Path to call CD-HIT-EST")
    parser.add_argument("--cdhit_args", default = "", required = False, help = "A string of arguments for CD-HIT-EST")
    
    return parser.parse_args()


def main():
    args = get_arguments()
    subdirs = make_output_dirs(args.outdir)
    cls = read_cluster_contents(args.clusters)
    queries = extract_allele_seqs(cluster_contents = cls, db = args.allele_db, outdir = subdirs["allele"],\
                                  prefix = args.prefix)
    assemblies = makeblastdb(args.assemblies, subdirs["database"])
    

def makeblastdb(assemblies, db_dir):
    print("Creating BLAST database from every assembly file.")
    return


def extract_allele_seqs(cluster_contents, db, outdir, prefix):
    # Create a FASTA file for each cluster to store their allele sequences
    print("Extract allele sequences of each cluster.")
    allele_db = list(SeqIO.parse(db, "fasta"))  # import the allele database
    queries = {}
    file_count = 0
    for cid, alleles in cluster_contents.items():
        if prefix != "":
            out_fasta = os.path.join(outdir, "__".join([prefix, cid, "alleles.fasta"]))
        else:
            out_fasta = os.path.join(outdir, cid + "__alleles.fasta")
        queries[cid] = out_fasta
        out_handle = open(out_fasta, "w")  # create a new file or override an existing file
        alleles = alleles.split(",")
        sum_alleles = len(alleles)
        allele_count = 0
        for rec in allele_db:
            if rec.id in alleles:
                SeqIO.write(rec, out_handle, "fasta")
                allele_count += 1
        out_handle.close()
        file_count += 1
        if allele_count < sum_alleles:
            print("* Warning: not all alleles of cluster %s are found in the allele database." % cid)
    print("* Totally %i FASTA files are generated." % file_count)
    
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
