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
    parser.add_argument("--analysis_method", type=str)
    args = parser.parse_args()
    
    blastn_results_path = args.blastn_results
    sample_name = args.sample_name
    mode = args.mode
    analysis_method = args.analysis_method

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
    #print(blastn_results.dtypes)
    blastn_top_hit = blastn_results.drop_duplicates(subset=["qseqid"], keep="first").copy()
    blastn_top_hit.to_csv(sample_name + "_" + analysis_method + "_blastn_top_hits.txt", index=False, sep="\t")
    #extract all sequences showing top blast homology to virus or viroid hits.
    #it is important to include copy() here as otherwise it will complain further down with the error
    #A value is trying to be set on a copy of a slice from a DataFrame.

    blastn_viral_top_hit = blastn_top_hit[blastn_top_hit["species"].str.contains('virus|viroid', na=False)].copy()
    #only retain blast hits with qcovs < 90
    if mode == "ncbi":
        #blastn_viral_top_hit['FLAG'] = blastn_viral_top_hit['qcovs'].apply(lambda x: 'low_confidence' if x < 90 else '')
        #blastn_viral_top_hit['FLAG'] = np.where((blastn_viral_top_hit['qcovs'] < 90) ^ (blastn_viral_top_hit['pident'] < 85) , "FLAG", "")
        blastn_viral_top_hit_high_conf = blastn_viral_top_hit.drop(blastn_viral_top_hit[blastn_viral_top_hit["qcovs"] < 90].index)
        #blastn_viral_top_hit_high_conf = blastn_viral_top_hit.drop(blastn_viral_top_hit[blastn_viral_top_hit["qcovs"] > 90].index).sort_values(by=['evalue'], ascending=True)
    elif mode == "localdb":
        #blastn_viral_top_hit['FLAG'] = np.where((blastn_viral_top_hit['qcovs'] < 95) ^ (blastn_viral_top_hit['pident'] < 85) , "FLAG", "")
        #blastn_viral_top_hit_high_conf = blastn_viral_top_hit.drop(blastn_viral_top_hit[blastn_viral_top_hit["qcovs"] > 95].index).sort_values(by=['evalue'], ascending=True)
        blastn_viral_top_hit_high_conf = blastn_viral_top_hit.drop(blastn_viral_top_hit[blastn_viral_top_hit["qcovs"] < 95].index)
    #derive read/contig count per viral spp
    summary_per_spp = blastn_viral_top_hit['species'].value_counts().rename_axis('species').reset_index(name='count')
    summary_per_spp_high_conf = blastn_viral_top_hit_high_conf['species'].value_counts().rename_axis('species').reset_index(name='count')
    #summary_per_spp_high_conf = blastn_viral_top_hit_high_conf['species'].value_counts().rename_axis('species').reset_index(name='counts')
    #summary_per_spp.index.name = 'species'
    #summary_per_spp_high_conf.index.name = 'species'
    summary_per_spp.to_csv(sample_name + "_" + analysis_method + "_viral_spp_abundance.txt", index=False, sep="\t")

    spp = blastn_viral_top_hit[['sacc','species','qseqid']].copy()
    spp['count'] = spp.groupby(['species', 'sacc'])['qseqid'].transform('size')
    #collapse all contigs to given accession number and species
    f = lambda x: x.tolist() if len(x) > 1 else x
    spp = spp.groupby(['species','sacc', 'count'])['qseqid'].agg(f).reset_index().reindex(spp.columns, axis=1)
    
    #reorder columns before saving
    spp = spp[["species", "sacc", "count", "qseqid"]].sort_values(["count"], ascending=[False])

    spp.to_csv(sample_name + "_" + analysis_method + "_queryid_list_with_viral_match.txt", index=False, sep="\t",float_format="%.2f")
    #replace space with underscore
    #blastn_blastn_top_evalue.replace('Elephantopus_scaber_closterovirus', 'Citrus_tristeza_virus',regex=True,inplace=True)
    #blastn_blastn_top_evalue.replace('Hop_stunt_viroid_-_cucumber', 'Hop_stunt_viroid',regex=True,inplace=True)
    
    blastn_viral_top_hit.to_csv(sample_name + "_" + analysis_method + "_blastn_top_viral_hits.txt", index=False, sep="\t",float_format="%.2f")
    #just retain longest contig for each virus/viroid species
    blastn_viral_top_hit_spp= blastn_viral_top_hit.sort_values(["qlen"], ascending=[True]).groupby("species", as_index=False).first().copy()

    blastn_viral_top_hit_spp.to_csv(sample_name + "_" + analysis_method +  "_blastn_top_viral_spp_hits.txt", index=False, sep="\t",float_format="%.2f")

    blastn_viral_top_hit_f = blastn_viral_top_hit[["qseqid", "qlen", "species", "sacc", "stitle", "slen", "pident", "sstrand", "evalue", "bitscore", "qcovs"]]
    blastn_viral_top_hit_spp_evalue_based = blastn_viral_top_hit_f.sort_values(["evalue"], ascending=[True]).groupby("species", as_index=False).first().copy()
    blastn_viral_top_hit_spp_pident_based = blastn_viral_top_hit_f.sort_values(["pident"], ascending=[False]).groupby("species", as_index=False).first().copy()
    blastn_viral_top_hit_spp_length_based = blastn_viral_top_hit_f.sort_values(["qlen"], ascending=[False]).groupby("species", as_index=False).first().copy()
    blastn_viral_top_hit_spp_bitscore_based = blastn_viral_top_hit_f.sort_values(["bitscore"], ascending=[False]).groupby("species", as_index=False).first().copy()

    summary_per_spp = summary_per_spp.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')
    summary_per_spp_high_conf = summary_per_spp_high_conf.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">') # use bootstrap styling
    blastn_viral_top_hit_spp = blastn_viral_top_hit_spp.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')
    blastn_viral_top_hit_spp_evalue_based = blastn_viral_top_hit_spp_evalue_based.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')
    blastn_viral_top_hit_spp_pident_based = blastn_viral_top_hit_spp_pident_based.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')
    blastn_viral_top_hit_spp_length_based = blastn_viral_top_hit_spp_length_based.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')
    blastn_viral_top_hit_spp_bitscore_based = blastn_viral_top_hit_spp_bitscore_based.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')
    spp = spp.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')
    
    html_string = '''
    <html>
        <head>
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
            <style>body{ margin:0 100; background:whitesmoke; }
            .collapsible {
            background-color: #777;
            color: white;
            cursor: pointer;
            padding: 18px;
            width: 100%;
            border: none;
            text-align: left;
            outline: none;
            font-size: 15px;
            }
            .active, .collapsible:hover {
            background-color: #555;
            }
            .content {
            padding: 0 18px;
            display: none;
            overflow: hidden;
            background-color: #f1f1f1;
            }
            </style>
        </head>

        <body>
            <h1>Homology blast results</h1>

            <button type="button" class="collapsible"> Total number of match(es) to viral species</button>
            <div class="content">
                ''' + summary_per_spp + '''
            </div>

            <button type="button" class="collapsible"> Total number of match(es) to viral species (filtered) </button>
            <div class="content">
                <p>Only blast viral matches which show >90% query coverage for NCBI and >95% query coverage for local viral database were considered here.</p>
                ''' + summary_per_spp_high_conf + '''
            </div>

            <button type="button" class="collapsible"> Total number of match(es) to specific viral accession number </h2></button>
            <div class="content">
                ''' + spp + '''
            </div>

            <button type="button" class="collapsible"> Top viral species match(es) based on evalue </h2></button>
            <div class="content">
                ''' + blastn_viral_top_hit_spp_evalue_based + '''
            </div>

            <button type="button" class="collapsible"> Top viral species match(es) based on % identity (pident) </h2></button>
            <div class="content">
                ''' + blastn_viral_top_hit_spp_pident_based + '''
            </div>

            <button type="button" class="collapsible"> Top viral species match(es) based on query length (qlen) </h2></button>
            <div class="content">
                ''' + blastn_viral_top_hit_spp_length_based + '''
            </div>

            <button type="button" class="collapsible"> Top viral species match(es) based on bitscore </h2></button>
            <div class="content">
                ''' + blastn_viral_top_hit_spp_bitscore_based + '''
            </div>

        <script>
        var coll = document.getElementsByClassName("collapsible");
        var i;

        for (i = 0; i < coll.length; i++) {
        coll[i].addEventListener("click", function() {
            this.classList.toggle("active");
            var content = this.nextElementSibling;
            if (content.style.display === "block") {
            content.style.display = "none";
            } else {
            content.style.display = "block";
            }
        });
        }
        </script>

        </body>
    </html>'''

    report = open(sample_name + "_blast_report.html", "w")
    report.write(html_string)
    report.close()

if __name__ == "__main__":
    main()
