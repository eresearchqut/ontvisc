#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from functools import reduce
import glob
import re
import os
import time

def main():
    parser = argparse.ArgumentParser(description="Derive a qc report")
    parser.add_argument("--host_filtering", type=str)
    parser.add_argument("--adapter_trimming", type=str)
    parser.add_argument("--quality_trimming", type=str)
    args = parser.parse_args()
    host_filtering = args.host_filtering
    adapter_trimming = args.adapter_trimming
    quality_trimming = args.quality_trimming

    


    timestr = time.strftime("%Y%m%d-%H%M%S")

    summary_dict = {}

    #if os.path.isfile("*raw_NanoStats.txt"):
    for raw_read_out in glob.glob("*raw_NanoStats.txt"):
        raw_reads = ()
        line_number = 0
        sample = (os.path.basename(raw_read_out).replace('_raw_NanoStats.txt', ''))
        with open(raw_read_out, 'r') as f:
            for line in f:
                string_to_search = ("number_of_reads")
                line_number += 1
                if string_to_search in line:
                    elements = line.rsplit('\t')
                    raw_reads = int(elements[1].strip())
                    print(raw_reads)

            #first_line = next(f)
            #raw_reads = int(first_line[0].strip())
        summary_dict[sample] = [raw_reads]
        print(summary_dict)
        f.close()

    
    for qt_read_out in glob.glob("*filtered_NanoStats.txt"):
        qt_reads = ()
        line_number = 0 
        sample = (os.path.basename(qt_read_out).replace('_filtered_NanoStats.txt', ''))
        with open(qt_read_out, 'r') as f:
            for line in f:
                string_to_search = ("number_of_reads")
                line_number += 1
                if string_to_search in line:
                    elements = line.rsplit('\t')
                    qt_reads = int(elements[1].strip())
                    print(qt_reads)
            #first_line = next(f)
            #qt_reads = int(first_line[0].strip())
        f.close()
        summary_dict[sample].append(qt_reads)
    
    if host_filtering == "true":
        for host_filt_read_out in glob.glob("*unaligned_reads_count.txt"):
            hf_reads = ()
            sample = (os.path.basename(host_filt_read_out).replace('_unaligned_reads_count.txt', ''))
            with open(host_filt_read_out, 'r') as f:
                first_line = next(f)
                hf_reads = int(first_line.strip())
                print(hf_reads)
            f.close()
            summary_dict[sample].append(hf_reads)
            print(summary_dict)

    if (adapter_trimming == "true" or quality_trimming == "true") and host_filtering == "true":
        run_data_df = pd.DataFrame([([k] + v) for k, v in summary_dict.items()], columns=['Sample','raw_reads','quality_filtered_reads', 'host_filtered_reads'])
        run_data_df['percent_quality_filtered'] = run_data_df['quality_filtered_reads'] / run_data_df['raw_reads'] * 100
        run_data_df['percent_host_filtered'] = run_data_df['host_filtered_reads'] / run_data_df['raw_reads'] * 100
        run_data_df['percent_quality_filtered'] = run_data_df['percent_quality_filtered'].apply(lambda x: float("{:.2f}".format(x)))
        run_data_df['percent_host_filtered'] = run_data_df['percent_host_filtered'].apply(lambda x: float("{:.2f}".format(x)))
        run_data_df = run_data_df.sort_values("Sample")
        run_data_df.to_csv("run_qc_report_" + timestr + ".txt", index = None, sep="\t")
    else:
        if adapter_trimming == "true" or quality_trimming == "true":
            run_data_df = pd.DataFrame([([k] + v) for k, v in summary_dict.items()], columns=['Sample','raw_reads','quality_filtered_reads'])
            run_data_df['percent_quality_filtered'] = run_data_df['quality_filtered_reads'] / run_data_df['raw_reads'] * 100
            run_data_df['percent_quality_filtered'] = run_data_df['percent_quality_filtered'].apply(lambda x: float("{:.2f}".format(x)))
            run_data_df = run_data_df.sort_values("Sample")
            run_data_df.to_csv("run_qc_report_" + timestr + ".txt", index = None, sep="\t")
        elif host_filtering == "true":
            run_data_df = pd.DataFrame([([k] + v) for k, v in summary_dict.items()], columns=['Sample','raw_reads','host_filtered_reads'])
            run_data_df['percent_host_filtered'] = run_data_df['host_filtered_reads'] / run_data_df['raw_reads'] * 100
            run_data_df['percent_host_filtered'] = run_data_df['percent_host_filtered'].apply(lambda x: float("{:.2f}".format(x)))
            run_data_df = run_data_df.sort_values("Sample")
            run_data_df.to_csv("run_qc_report_" + timestr + ".txt", index = None, sep="\t")

if __name__ == '__main__':
    main()
