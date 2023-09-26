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
        blastn_results = blastn_results[~blastn_results["species"].str.contains("synthetic construct", na=False)]
        

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
    blastn_viral_top_hit = blastn_top_hit[blastn_top_hit["species"].str.contains('virus|viroid', na=False)]
    #only retain blast hits with qcovs < 90
    blastn_viral_top_hit.drop(blastn_viral_top_hit[blastn_viral_top_hit["qcovs"] < 90].index, inplace = True)
    #blastn_viral_top_hit['count'] = blastn_viral_top_hit.groupby('species')['species'].transform('count')
    #derive read/contig count per viral spp
    summary_per_spp = blastn_viral_top_hit['species'].value_counts()
    summary_per_spp.index.name = 'species'
    summary_per_spp.to_csv(sample_name + "_viral_spp_abundance.txt", index=True, sep="\t", header=["Count"])

    #summary_per_qseqid = blastn_viral_top_hit[['species','sacc']].value_counts()
    #summary_per_qseqid.to_csv(sample_name + "_viral_spp_by_ids_abundance.txt", index=True, sep="\t", header=["Count"])

    spp = blastn_viral_top_hit[['sacc','species','qseqid']]
    #spp = spp.assign(count=spp.groupby(['species', 'sacc'])['qseqid'].transform('size'))
    #spp['count'] = spp.groupby(['species', 'sacc']).transform('count')
    #grouped = df.groupby(['date', 'store']).size().reset_index(name='count')
    spp['count'] = spp.groupby(['species', 'sacc'])['qseqid'].transform('size')
    #spp = blastn_viral_top_hit[["sacc", "species", "qseqid"]]
    #collapse all contigs to given accession number and species
    f = lambda x: x.tolist() if len(x) > 1 else x
    spp = spp.groupby(['species','sacc', 'count'])['qseqid'].agg(f).reset_index().reindex(spp.columns, axis=1)
    
    #reorder columns before saving
    spp = spp[["species", "sacc", "count", "qseqid"]].sort_values(["count"], ascending=[False])
    
    print(spp)
    spp.to_csv(sample_name +  "_contigs_list_with_viral_match.txt", index=False, sep="\t",float_format="%.2f")
    #replace space with underscore
    #blastn_blastn_top_evalue.replace('Elephantopus_scaber_closterovirus', 'Citrus_tristeza_virus',regex=True,inplace=True)
    #blastn_blastn_top_evalue.replace('Hop_stunt_viroid_-_cucumber', 'Hop_stunt_viroid',regex=True,inplace=True)
    
    blastn_viral_top_hit.to_csv(sample_name +  "_blastn_top_viral_hits.txt", index=False, sep="\t",float_format="%.2f")
    #just retain longest contig for each virus/viroid species
    blastn_viral_top_hit_spp= blastn_viral_top_hit.sort_values(["evalue", "length"], ascending=[True, False])
    #print(blastn_viral_top_hit.head())
    blastn_viral_top_hit_spp.to_csv(sample_name + "_blastn_top_viral_spp_hits.txt", index=False, sep="\t",float_format="%.2f")
    #print(blastn_viral_top_hit_spp.dtypes)

    




                          
if __name__ == "__main__":
    main()   