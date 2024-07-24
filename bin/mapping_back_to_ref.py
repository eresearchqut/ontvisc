#!/usr/bin/env python
import argparse
from pathlib import Path
import subprocess
from functools import reduce
import glob
from subprocess import run, PIPE

def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast results")

    # All the required arguments #
    parser.add_argument("--fastq", type=str)
    args = parser.parse_args()
    fastq = args.fastq

    for reference in glob.glob("*.fasta"):
        #file_name = reference.removesuffix('.fasta')
        file_name = str(Path(reference).with_suffix(""))
        print(file_name)

        print("Aligning original reads")
        minimap2_output = str(file_name + ".sam")
        aligning = ["minimap2", "-ax", "map-ont", "--sam-hit-only", "-L", reference, fastq]
        subprocess.call(aligning, stdout=open(minimap2_output,"w"))

        print("Derive a bam file")
        bamoutput = str(file_name + ".bam")
        derivebam = ["samtools", "view", "-F", "4", "-bS", minimap2_output]
        subprocess.call(derivebam, stdout=open(bamoutput,"w"))

        print("Confirm we recover mapped reads")
        flagstatoutput = str(file_name + ".flagstat")
        flagstat = ["samtools", "flagstat", bamoutput]
        subprocess.call(flagstat, stdout=open(flagstatoutput,"w"))
        with open(flagstatoutput) as myfile:
            if '0 + 0 mapped' in myfile.read():
                print('No reads mapping!')
                subprocess.call(["rm","-r", minimap2_output])
                subprocess.call(["rm","-r", bamoutput])
                subprocess.call(["rm","-r", flagstatoutput])
            else:
                print("Sorting bam file")
                sortedbamoutput = str(file_name + ".sorted.bam")
                sorting = ["samtools", "sort", bamoutput, "-o", sortedbamoutput]
                subprocess.call(sorting)

                print("Indexing bam file")
                bamindex = str(file_name + ".sorted.bam.bai")
                indexing = ["samtools", "index", sortedbamoutput]
                subprocess.call(indexing, stdout=open(bamindex,"w"))

if __name__ == "__main__":
    main()