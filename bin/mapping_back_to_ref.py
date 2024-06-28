#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import os
import subprocess
from functools import reduce
import glob
from subprocess import run, PIPE

def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast results")

    # All the required arguments #
    parser.add_argument("--results", type=str)
    parser.add_argument("--fastq", type=str)
    parser.add_argument("--sample", type=str)
    parser.add_argument("--blastdbpath", type=str)
    args = parser.parse_args()
    
    results_path = args.results
    sample = args.sample
    fastq = args.fastq
    blastdbpath = args.blastdbpath

    for blast_results in glob.glob(results_path):
        blast_df = pd.read_csv(blast_results, header=0, sep="\t",index_col=None)
        print(blast_df)
        target_dict = {}
        target_dict = pd.Series(blast_df.species.values,index=blast_df.sacc).to_dict()
        print (target_dict)

        for refid, refspname in target_dict.items():
            print (refid)
            print (refspname)
            combinedid = str(refid + " " + refspname).replace("sp.","sp").replace(" ","_")

            print("Extract sequence from blast database")
            fastafile = (sample + "_" + combinedid + ".fa").replace(" ","_")

            single_fasta_entry = open(fastafile, "w")
            command_line = ["blastdbcmd","-db", blastdbpath, "-entry", refid, \
                            "-outfmt","'%f'"]
            subprocess.call(command_line, stdout=single_fasta_entry)
            single_fasta_entry.close()

            print("Aligning original reads")
            index=(sample + "_" + combinedid).replace(" ","_")
            minimap2_output = str(index + ".sam")
            aligning = ["minimap2", "-ax", "map-ont", "--MD", "--sam-hit-only", fastafile, fastq]
            subprocess.call(aligning, stdout=open(minimap2_output,"w"))

            print("Derive a bam file")
            bamoutput = str(index + ".bam")
            #derivebam = ["samtools", "view", "-F", "-@", cpus, "-bS", samoutput]
            derivebam = ["samtools", "view", "-F", "4", "-bS", minimap2_output]
            subprocess.call(derivebam, stdout=open(bamoutput,"w"))

            print("Sorting bam file")
            sortedbamoutput = str(index + ".sorted.bam")
            #sorting = ["samtools", "sort", "-@", cpus, bamoutput, "-o", sortedbamoutput]
            sorting = ["samtools", "sort", bamoutput, "-o", sortedbamoutput]
            subprocess.call(sorting)

            print("Indexing bam file")
            bamindex = str(index + ".sorted.bam.bai")
            indexing = ["samtools", "index", sortedbamoutput]
            subprocess.call(indexing, stdout=open(bamindex,"w"))

            print("Derive coverage")
            derivecov = ["mosdepth", "--thresholds", "1,10,20,30", sample + "_" + combinedid, sortedbamoutput]
            subprocess.call(derivecov)

            #samtools coverage ${sampleid}_aln.sorted.bam > ${sampleid}_histogram.txt  > ${sampleid}_coverage.txt
            #samtools coverage -A -w 50 ${sampleid}_aln.sorted.bam > ${sampleid}_histogram



if __name__ == "__main__":
    main()
