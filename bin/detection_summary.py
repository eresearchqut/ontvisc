#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from functools import reduce
import glob
import time


def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load VirReport pipeline results")
    # All the required arguments #
    parser.add_argument("--threshold", type=float)
    
    args = parser.parse_args()
    threshold = args.threshold
    
#     timestr = time.strftime("%Y%m%d-%H%M%S")

    run_data = pd.DataFrame()
    for fl in glob.glob("*_top_blast_with_cov_stats.txt"):
        sample_data = pd.read_csv(fl, header=0, sep="\t",index_col=None)
        
        run_data = pd.concat([run_data, sample_data], axis = 0)
    print(run_data)
    
    
    run_data = run_data[["Sample","species","stitle","qseqid","sacc","length","pident","sstrand","evalue","bitscore","bitscore","qcovs","read_count","mean_cov","RPKM","PCT_5X","PCT_10X","PCT_20X"]]
    contamination_flag(run_data,threshold)
    
    run_data.to_csv("summary_detection.txt", index=None, sep="\t")


def contamination_flag(df, threshold):
    df["RPKM"] = df["RPKM"].astype(float)
    df["RPKM_max"] = df.groupby(["species"])["RPKM"].transform(max)
    df["RPKM_max"] =  df["RPKM_max"].round(1)
    df["threshold_value"]=df["RPKM_max"]*threshold
    df["threshold_value"] = df["threshold_value"].round(1)
    df["contamination_flag"] = np.where(df["RPKM"] <= df["threshold_value"], True, False)
    df["contamination_flag"] = np.where(df["RPKM_max"] <= 10, "NA", df["contamination_flag"])
    df = df.sort_values(["Sample", "stitle"], ascending = (True, True))

if __name__ == "__main__":
    main()