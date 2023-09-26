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
    parser.add_argument("--mode", type=str)
    args = parser.parse_args()
    
    blastn_results_path = args.blastn_results
    sample_name = args.sample_name
    mode = args.mode

    if mode == "ncbi":
        blastn_results = pd.read_csv(blastn_results_path, sep="\t", index_col=False, names=["qseqid", "sgi", "sacc", "length", "pident", "mismatch", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe", "species"], dtype={"stitle": 'str', "staxids": 'str', "species": 'str'})
        
        #remove synthetic construct hits
        blastn_results = blastn_results[~blastn_results["species"].str.contains("synthetic construct")]
        

    elif mode == "localdb":
        #retrieve spp name and accession from local db fasta header
        #rearrange column so it matches the one for NCBI
        blastn_results = pd.read_csv(blastn_results_path, sep="\t", index_col=False, names=["qseqid", "sgi", "seq_desc", "length", "pident", "mismatch", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe"], dtype={"stitle": 'str', "staxids": 'str'})
        blastn_results['sacc'] = blastn_results['seq_desc'].str.split('|').str[0]
        blastn_results['species'] = blastn_results['seq_desc'].str.split('|').str[1]
        blastn_results['species'] = blastn_results['species'].str.replace("Species:","")
        blastn_results = blastn_results[["qseqid", "sgi", "sacc", "length", "pident", "mismatch", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe", "species"]]

    blastn_top_hit = blastn_results.drop_duplicates(subset=["qseqid"], keep="first")    
    blastn_top_hit.to_csv(sample_name +  "_blastn_top_hits.txt", index=False, sep="\t",float_format="%.2f")
    #only retain intagrated and episomal viral hits
    #l=['virus', 'viroid']
    
    #blastn_viral_top_hit = blastn_top_hit[blastn_top_hit['stitle'].str.contains('|'.join(l))]
    #extract all sequences showing top blast homology to virus or viroid hits.
    blastn_viral_top_hit = blastn_top_hit[blastn_top_hit["species"].str.contains('virus|viroid')]
    blastn_viral_top_hit.drop(blastn_viral_top_hit[blastn_viral_top_hit["qcovs"] < 90].index, inplace = True)
    spp = blastn_viral_top_hit[['sacc','species','qseqid']]
    #spp = blastn_viral_top_hit[["sacc", "species", "qseqid"]]
    #collapse all contigs to given accession number and species
    f = lambda x: x.tolist() if len(x) > 1 else x
    spp = spp.groupby(['species','sacc'])['qseqid'].agg(f).reset_index().reindex(spp.columns, axis=1)
    #remove brackets from list in qseqid column
    spp['qseqid'] = spp['qseqid'].str.join(',')
    #reorder columns before saving
    spp = spp[["species", "sacc", "qseqid"]]
    spp.to_csv(sample_name +  "_contig_list_with_viral_match.txt", index=False, sep="\t",float_format="%.2f")
    #replace space with underscore
    #blastn_blastn_top_evalue.replace('Elephantopus_scaber_closterovirus', 'Citrus_tristeza_virus',regex=True,inplace=True)
    #blastn_blastn_top_evalue.replace('Hop_stunt_viroid_-_cucumber', 'Hop_stunt_viroid',regex=True,inplace=True)
    
    #remove hits that cover less than 50% of query
    
    blastn_viral_top_hit.to_csv(sample_name +  "_blastn_top_viral_hits.txt", index=False, sep="\t",float_format="%.2f")
        #just retain longest contig for each virus/viroid species
    blastn_viral_top_hit_spp= blastn_viral_top_hit.sort_values(["evalue", "length"], ascending=[True, False]).groupby("species").first().reset_index()
    #print(blastn_viral_top_hit.head())
    blastn_viral_top_hit_spp.to_csv(sample_name + "_blastn_top_viral_spp_hits.txt", index=False, sep="\t",float_format="%.2f")
    #print(blastn_viral_top_hit_spp.dtypes)

    
    
# elif mode == "localdb":
#     blastn_results = pd.read_csv(blastn_results_path, sep="\t", index_col=False, names=["qseqid", "sgi", "seq_desc", "length", "pident", "mismatch", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe"], dtype={"stitle": 'str', "staxids": 'str'})
#     blastn_results['sacc'] = blastn_results['seq_desc'].str.split('|').str[0]
#     blastn_results['species'] = blastn_results['seq_desc'].str.split('|').str[1]
#     blastn_results['species'] = blastn_results['species'].str.replace("Species:","")
#     blastn_results = blastn_results[["qseqid", "sgi", "sacc", "length", "pident", "mismatch", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe", "species"]]
    
#     blastn_viral_top_hit = blastn_results[blastn_results["species"].str.contains('virus|viroid')]
#     blastn_viral_top_hit.drop(blastn_viral_top_hit[blastn_viral_top_hit["qcovs"] < 50].index, inplace = True)
#     print(blastn_viral_top_hit)

                          
if __name__ == "__main__":
    main()   