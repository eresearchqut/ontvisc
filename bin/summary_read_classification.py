#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import os
import subprocess
from functools import reduce
from glob import glob
from subprocess import run, PIPE

def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load coverage stats files")

    # All the required arguments #
    parser.add_argument("--sample", type=str)
    
    
        
#    parser.add_argument("--fasta", type=str)
#    parser.add_argument("--blast_results", type=str)
#    parser.add_argument("--mosdepth_results", type=str)
#    parser.add_argument("--coverm-results", type=str)
    

    args = parser.parse_args()
    
#    blast = args.blast_results
#    mosdepth = args.mosdepth_results
#    coverm = args.coverm-results
    sample_name = args.sample
    print(sample_name)
    PCT_10X_all = ()
    coverm_all = pd.DataFrame()
    blast_df = pd.DataFrame()
    for blast_results in glob("*_blastn_top_viral_spp_hits.txt"):
        blastn_results = pd.read_csv(blast_results, sep="\t", index_col=False)
        blast_df = blastn_results[["species", "stitle", "qseqid", "sacc", "length", "pident", "sstrand", "evalue", "bitscore", "qcovs"]]
        sacc_list = blast_df["sacc"].tolist()
        for sacc in sacc_list:
            for coverm_results in glob("*_coverm_summary.txt"):
                if sacc in str(coverm_results):
                    coverm_results = pd.read_csv(coverm_results, sep="\t", index_col=False)
                    coverm_results.columns = ["genome", "read_counts", "mean_cov", "variance", "RPKM", "%_bases_cov", "length"]
                    coverm_results["sacc"] = sacc
                    coverm_df=coverm_results[["sacc", "read_counts", "mean_cov", "RPKM"]]
                    #print(coverm_df)
                    coverm_all = coverm_all.append(coverm_df)
            for mosdepth_results in glob("*mosdepth.global.dist.txt"):
                if sacc in str(mosdepth_results):
                    mosdepth_results = pd.read_csv(mosdepth_results, sep="\t", index_col=False)
                    mosdepth_results.columns = ["genome", "pc_coverage", "depth"]
                    PCT_10X = mosdepth_results[(mosdepth_results['pc_coverage']==10) & (mosdepth_results['genome']==sacc)]
                    print(PCT_10X)
                    #PCT_10X_all = PCT_10X_all.append(PCT_10X)



        blast_df = pd.merge(blast_df, coverm_all, on='sacc')
        print(blast_df)

    #def cov_stats(blastdbpath, cpus, dedup, fastqfiltbysize, final_data, rawfastq, read_size, sample, target_dict, mode):
    #print("Align reads and derive coverage and depth for best hit")
    #rawfastq_read_counts = (len(open(rawfastq).readlines(  ))/4)
        blast_df.to_csv(sample_name +  "_top_scoring_targets_with_cov_stats.txt", index=None, sep="\t",float_format="%.2f")

    #cov_dict = {}
    #dedup_read_counts_dict = {}
    #dup_pc_dict = {}
    #fpkm_dict = {}
    #PCT_1X_dict = {}
    #PCT_5X_dict = {}
    #PCT_10X_dict = {}
    #PCT_20X_dict = {}
    ##read_counts_dict = {}
    #rpm_dict = {}
    #consensus_dict = {}

    #read_counts_dedup_df = pd.DataFrame()
    #dup_pc_df = pd.DataFrame()
    #cov_df = pd.DataFrame()
    #read_counts_df = pd.DataFrame()
    #rpm_df = pd.DataFrame()
    #fpkm_df = pd.DataFrame()
    #PCT_1X_df = pd.DataFrame()
    #PCT_5X_df = pd.DataFrame()
    #PCT_10X_df = pd.DataFrame()
    #PCT_20X_df = pd.DataFrame()
    #consensus_df = pd.DataFrame()


    #reflen = ()
    #cov = ()
    #PCT_1X = ()
    #PCT_5X = ()
    #PCT_10X = ()
    #PCT_20X = ()

    #with open(picard_output) as f:
#a = " "
#while(a):
#    a = f.readline()
#    l = a.find("MEAN_COVERAGE") #Gives a non-negative value when there is a match
#    if ( l >= 0 ):
#        line = f.readline()
#        elements = line.split("\t")
#        reflen, cov, PCT_1X, PCT_5X, PCT_10X, PCT_20X = elements[0], elements[1], elements[13], elements[14], elements[15],elements[17]
#f.close()
#cov_dict[refspname] = cov
#PCT_1X_dict[refspname] = PCT_1X
#PCT_5X_dict[refspname] = PCT_5X
#PCT_10X_dict[refspname] = PCT_10X
#PCT_20X_dict[refspname] = PCT_20X

#fpkm = round(int(final_read_counts)/(int(reflen)/1000*int(rawfastq_read_counts)/1000000))
#rpm = round(int(final_read_counts)*1000000/int(rawfastq_read_counts))

#rpm_dict[refspname] = rpm
#fpkm_dict[refspname] = fpkm

#cov_df = pd.DataFrame(cov_dict.items(),columns=["Species_updated", "mean_read_depth"])
#read_counts_df = pd.DataFrame(read_counts_dict.items(),columns=["Species_updated", "read_count"])
#rpm_df = pd.DataFrame(rpm_dict.items(),columns=["Species_updated", "RPM"])
#fpkm_df = pd.DataFrame(fpkm_dict.items(),columns=["Species_updated", "FPKM"])
#PCT_1X_df = pd.DataFrame(PCT_1X_dict.items(),columns=["Species_updated", "PCT_1X"])
#PCT_5X_df = pd.DataFrame(PCT_5X_dict.items(),columns=["Species_updated", "PCT_5X"])
#PCT_10X_df = pd.DataFrame(PCT_10X_dict.items(),columns=["Species_updated", "PCT_10X"])
#PCT_20X_df = pd.DataFrame(PCT_20X_dict.items(),columns=["Species_updated", "PCT_20X"])
#consensus_df = pd.DataFrame(consensus_dict.items(),columns=["Species_updated", "consensus_fasta"])

#project_files = glob(index + "*ebwt") + glob(index + ".vcf.gz*")
#    for fl in project_files:
#        subprocess.call(["rm","-r", fl])
#except OSError as err:
#    print("OS error: {0}".format(err))

#print("Deriving summary table with coverage statistics")

#if read_counts_dedup_df.empty:
#dfs = [final_data, cov_df, read_counts_df, rpm_df, fpkm_df, PCT_1X_df, PCT_5X_df, PCT_10X_df, PCT_20X_df, consensus_df]
#else:
#dfs = [final_data, cov_df, read_counts_df, read_counts_dedup_df, dup_pc_df, rpm_df, fpkm_df, PCT_1X_df, PCT_5X_df, PCT_10X_df, PCT_20X_df, consensus_df]

#full_table = reduce(lambda left,right: pd.merge(left,right,on=["Species_updated"],how='outer'), dfs)

#full_table["mean_read_depth"] = full_table["mean_read_depth"].astype(float)
#full_table["PCT_1X"] = full_table["PCT_1X"].astype(float)
#full_table["PCT_5X"] = full_table["PCT_5X"].astype(float)
#full_table["PCT_10X"] = full_table["PCT_10X"].astype(float)
#full_table["PCT_20X"] = full_table["PCT_20X"].astype(float)
#if "duplication_rate" in full_table.columns:
#full_table["duplication_rate"] = full_table["duplication_rate"].astype(float)
#full_table.insert(0, "Sample", sample)

#if mode == 'ncbi':
#full_table = full_table.drop(["Species"], axis=1)
#full_table = full_table.rename(columns={"Species_updated": "Species"})
#full_table.to_csv(sample + "_" + read_size + "_top_scoring_targets_with_cov_stats.txt", index=None, sep="\t",float_format="%.2f")

#elif mode == 'viral_db':
#full_table = full_table.rename(columns={"Species_updated": "Species"})
#full_table.to_csv(sample + "_" + read_size + "_top_scoring_targets_with_cov_stats_viral_db.txt", index=None, sep="\t",float_format="%.2f")
if __name__ == "__main__":
    main()