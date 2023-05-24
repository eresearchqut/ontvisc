#!/usr/bin/env python
import pandas as pd
import argparse
import numpy as np



def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast results")

    # All the required arguments #
    parser.add_argument("--megablast_results", type=str)
    parser.add_argument("--sample_name", type=str)
    args = parser.parse_args()
    
    megablast_results_path = args.megablast_results
    sample_name = args.sample_name

    megablast_results = pd.read_csv(megablast_results_path, sep="\t", index_col=False, names=["qseqid", "sgi", "sacc", "length", "pident", "mismatch", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe", "species"], dtype={"stitle": 'str', "staxids": 'str', "species": 'str'})
    
    megablast_results = megablast_results[~megablast_results["species"].str.contains("synthetic construct")]
    megablast_top_hit = megablast_results.drop_duplicates(subset=["qseqid"], keep="first")
    

    
    megablast_top_hit.to_csv(sample_name +  "_blastn_vs_NT_top_hits.txt", index=False, sep="\t")
    #only retain intagrated and episomal viral hits
    #l=['virus', 'viroid']
    
    #megablast_viral_top_hit = megablast_top_hit[megablast_top_hit['stitle'].str.contains('|'.join(l))]
    megablast_viral_top_hit = megablast_top_hit[megablast_top_hit["species"].str.contains('virus|viroid')]
    #replace space with underscore
    #megablast_blastn_top_evalue.replace('Elephantopus_scaber_closterovirus', 'Citrus_tristeza_virus',regex=True,inplace=True)
    #megablast_blastn_top_evalue.replace('Hop_stunt_viroid_-_cucumber', 'Hop_stunt_viroid',regex=True,inplace=True)
    
    megablast_viral_top_hit.to_csv(sample_name +  "_blastn_vs_NT_top_viral_hits.txt", index=False, sep="\t")
    megablast_viral_top_hit_spp= megablast_viral_top_hit.sort_values(["species", "length"], ascending=[True, False]).groupby("species").first().reset_index()
    print(megablast_viral_top_hit.head())
    megablast_viral_top_hit_spp.to_csv(sample_name +  "_blastn_vs_NT_top_viral_spp_hits.txt", index=False, sep="\t")

                          
if __name__ == "__main__":
    main()   