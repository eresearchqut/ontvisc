#!/usr/bin/env python
import pandas as pd
import argparse
import numpy as np



def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast results")

    # All the required arguments #
    parser.add_argument("--blastn_results", type=str)
    parser.add_argument("--sample_name", type=str)
    args = parser.parse_args()
    
    blastn_results_path = args.blastn_results
    sample_name = args.sample_name

    blastn_results = pd.read_csv(blastn_results_path, sep="\t", index_col=False, names=["qseqid", "sgi", "sacc", "length", "pident", "mismatch", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe", "species"], dtype={"stitle": 'str', "staxids": 'str', "species": 'str'})
    
    #remove synthetic construct hits
    blastn_results = blastn_results[~blastn_results["species"].str.contains("synthetic construct")]
    blastn_top_hit = blastn_results.drop_duplicates(subset=["qseqid"], keep="first")
    

    
    blastn_top_hit.to_csv(sample_name +  "_blastn_vs_NT_top_hits.txt", index=False, sep="\t")
    #only retain intagrated and episomal viral hits
    #l=['virus', 'viroid']
    
    #blastn_viral_top_hit = blastn_top_hit[blastn_top_hit['stitle'].str.contains('|'.join(l))]
    blastn_viral_top_hit = blastn_top_hit[blastn_top_hit["species"].str.contains('virus|viroid')]
    #replace space with underscore
    #blastn_blastn_top_evalue.replace('Elephantopus_scaber_closterovirus', 'Citrus_tristeza_virus',regex=True,inplace=True)
    #blastn_blastn_top_evalue.replace('Hop_stunt_viroid_-_cucumber', 'Hop_stunt_viroid',regex=True,inplace=True)
    #extract all seqeunces showing top blast homology to virus or virois hits.
    blastn_viral_top_hit.to_csv(sample_name +  "_blastn_vs_NT_top_viral_hits.txt", index=False, sep="\t")
    #just retain longest contig for each virus/viroid species
    blastn_viral_top_hit_spp= blastn_viral_top_hit.sort_values(["species", "qlen"], ascending=[True, False]).groupby("species").first().reset_index()
    print(blastn_viral_top_hit.head())
    blastn_viral_top_hit_spp.to_csv(sample_name +  "_blastn_vs_NT_top_viral_spp_hits.txt", index=False, sep="\t")

                          
if __name__ == "__main__":
    main()   