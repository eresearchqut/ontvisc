#!/usr/bin/env python
import pandas as pd
import argparse
import glob


def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast results")
    # All the required arguments #
    parser.add_argument("--sample", type=str)
    args = parser.parse_args()
    sample_name = args.sample

    for bracken_file in glob.glob('*_bracken_report_viral.txt'):
        bracken_df = pd.read_csv(bracken_file, sep="\t", index_col=False)
        bracken_df["fraction_total_reads"] = pd.to_numeric(bracken_df["fraction_total_reads"], errors='coerce', downcast="float")
        bracken_df_filtered = bracken_df.drop(bracken_df[bracken_df["fraction_total_reads"] < 0.0001].index).sort_values(by=['kraken_assigned_reads'], ascending=False)
        bracken_df_html = bracken_df.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">') # use bootstrap styling
        bracken_df_filtered_html = bracken_df_filtered.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">') # use bootstrap styling
        
    for kaiju_file in glob.glob('*_kaiju_summary_viral.tsv'):
        kaiju_df = pd.read_csv(kaiju_file, sep="\t", index_col=False)
        kaiju_df["percent"] = pd.to_numeric(kaiju_df["percent"], errors='coerce', downcast="float")
        kaiju_df_filtered = kaiju_df.drop(kaiju_df[kaiju_df["percent"] < 0.05].index).sort_values(by=['reads'], ascending=False)
        kaiju_df_html = kaiju_df.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">') # use bootstrap styling
        kaiju_df_filtered_html = kaiju_df_filtered.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">') # use bootstrap styling 

    for blastn_file in glob.glob('*_viral_spp_abundance.txt'):
        for blastn_file_filt in glob.glob('*_viral_spp_abundance_filtered.txt'):
            blastn_df = pd.read_csv(blastn_file, sep="\t", index_col=False)
            megablast_summary_per_spp = blastn_df.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')
            blastn_df_high_conf = pd.read_csv(blastn_file_filt, sep="\t", index_col=False)
            megablast_summary_per_spp_high_conf = blastn_df_high_conf.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')

    #consider all options
    #A-B-C kaiju braken megablast
    #A-B kaiju braken done
    #B-C kaiju megablast
    #A-C
    #A kaiju done
    #B braken done
    #C homology search
    if glob.glob("*_bracken_report_viral.txt") and glob.glob("*_kaiju_summary_viral.tsv") and not glob.glob("*_viral_spp_abundance*.txt"):
        html_string = '''
            <html>
                <head>
                    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
                    <style>
                    body { margin:0 100; background:whitesmoke; 
                    }
                    container {
                        border: 1px solid grey;
                        margin: 1rem;
                    }
                    [data-tab-info] {
                        display: none;
                    }
 
                    .active[data-tab-info] {
                        display: block;
                    }
 
                    .tab-content {
                        margin-top: 1rem;
                        padding-left: 1rem;
                        font-size: 20px;
                        font-family: sans-serif;
                        font-weight: bold;
                        color: rgb(0, 0, 0);
                    }
 
                    .tabs {
                        border-bottom: 1px solid grey;
                        background-color: rgb(170, 245, 144);
                        font-size: 25px;
                        color: rgb(0, 0, 0);
                        display: flex;
                        margin: 0;
                    }
 
                    .tabs span {
                        background: rgb(170, 245, 144);
                        padding: 10px;
                        border: 1px solid rgb(255, 255, 255);
                    }
            
                    .tabs span:hover {
                        background: rgb(55, 219, 46);
                        cursor: pointer;
                        color: black;
                    }
                    </style>
                </head>
                <body>
                    <!-- Body Container -->
                    <div class="container">
                    
                        <!-- Tabs Detail -->
                        <div class="tabs">
                            <span data-tab-value="#tab_1">Kraken</span>
                            <span data-tab-value="#tab_2">Kaiju</span>
                        </div>

                        <!-- Tab content -->
                        <div class="tab-content">
                            <div class="tabs__tab active" id="tab_1" data-tab-info>
                                <h1>Direct read taxonomic classification (nucleotide-based)</h1>
                                <!-- *** Section 1 *** --->
                                <h2>Section 1: Kraken2 and Bracken viral results filtered </h2>
                                <p>Only viral matches which represented >=0.0001 of the total read fraction are retained </p>
                                ''' + bracken_df_filtered_html + '''

                                <!-- *** Section 2 *** --->
                                <h2>Section 2: All Kraken2 and Bracken viral results </h2>
                                ''' + bracken_df_html + '''
                            </div>
                            <div class="tabs__tab" id="tab_2" data-tab-info>
                                <h1>Direct read taxonomic classification (protein-based)</h1>
                                <!-- *** Section 1 *** --->
                                <h2>Section 1: Kaiju results filtered </h2>
                                <p>Only viral matches which represented >=0.05 of the total read fraction are retained </p>
                                ''' + kaiju_df_filtered_html + '''

                                <!-- *** Section 2 *** --->
                                <h2>Section 2: All Kaiju results </h2>
                                ''' + kaiju_df_html + '''
                            </div>
                        </div>
                    </div>
                     <script type="text/javascript">
       
                        // function to get each tab details
                        const tabs = document.querySelectorAll('[data-tab-value]')
                        const tabInfos = document.querySelectorAll('[data-tab-info]')
                
                        tabs.forEach(tab => {
                            tab.addEventListener('click', () => {
                                const target = document
                                    .querySelector(tab.dataset.tabValue);
                                tabInfos.forEach(tabInfo => {
                                    tabInfo.classList.remove('active')
                                })
                                target.classList.add('active');
                            })
                        })
                    </script>
                </body>
            </html>'''

        report = open(sample_name + "_read_classification_report.html", "w")
        report.write(html_string)
        report.close()

    elif glob.glob("*_bracken_report_viral.txt") and glob.glob("*_kaiju_summary_viral.tsv") and glob.glob("*_viral_spp_abundance*.txt"):
        html_string = '''
            <html>
                <head>
                    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
                    <style>
                    body { margin:0 100; background:whitesmoke; 
                    }
                    container {
                        border: 1px solid grey;
                        margin: 1rem;
                    }
                    [data-tab-info] {
                        display: none;
                    }
 
                    .active[data-tab-info] {
                        display: block;
                    }
 
                    .tab-content {
                        margin-top: 1rem;
                        padding-left: 1rem;
                        font-size: 20px;
                        font-family: sans-serif;
                        font-weight: bold;
                        color: rgb(0, 0, 0);
                    }
 
                    .tabs {
                        border-bottom: 1px solid grey;
                        background-color: rgb(170, 245, 144);
                        font-size: 25px;
                        color: rgb(0, 0, 0);
                        display: flex;
                        margin: 0;
                    }
 
                    .tabs span {
                        background: rgb(170, 245, 144);
                        padding: 10px;
                        border: 1px solid rgb(255, 255, 255);
                    }
            
                    .tabs span:hover {
                        background: rgb(55, 219, 46);
                        cursor: pointer;
                        color: black;
                    }
                    </style>
                </head>
                <body>
                    <!-- Body Container -->
                    <div class="container">
                        <!-- Tabs Detail -->
                        <div class="tabs">
                            <span data-tab-value="#tab_1">Megablast</span>
                            <span data-tab-value="#tab_2">Kraken</span>
                            <span data-tab-value="#tab_3">Kaiju</span>
                            
                        </div>
                        <!-- Tab content -->
                        <div class="tab-content">
                            <div class="tabs__tab active" id="tab_1" data-tab-info>
                                <h1>Direct read homology search (megablast)</h1>
                                <!-- *** Section 1 *** --->
                                <h2>Section 1: Number of reads matching to viral species filtered </h2>
                                <p>Only blast viral matches which show >90% query coverage for NCBI and >95% query coverage for local viral database were considered here.</p>
                                ''' + megablast_summary_per_spp_high_conf + '''
                                <h2>Section 2: Number of reads matching to viral species </h2>
                                <p>All blast viral matches were considered here.</p>
                                ''' + megablast_summary_per_spp + '''
                            </div>
                            <div class="tabs__tab" id="tab_2" data-tab-info>
                                <h1>Direct read taxonomic classification (nucleotide-based)</h1>
                                <!-- *** Section 1 *** --->
                                <h2>Section 1: Kraken2 and Bracken viral results filtered </h2>
                                <p>Only viral matches which represented >=0.0001 of the total read fraction are retained </p>
                                ''' + bracken_df_filtered_html + '''
                                <!-- *** Section 2 *** --->
                                <h2>Section 2: All Kraken2 and Bracken viral results </h2>
                                ''' + bracken_df_html + '''
                            </div>
                            <div class="tabs__tab" id="tab_3" data-tab-info>
                                <h1>Direct read taxonomic classification (protein-based)</h1>
                                <!-- *** Section 1 *** --->
                                <h2>Section 1: Kaiju results filtered </h2>
                                <p>Only viral matches which represented >=0.05 of the total read fraction are retained </p>
                                ''' + kaiju_df_filtered_html + '''
                                <!-- *** Section 2 *** --->
                                <h2>Section 2: All Kaiju results </h2>
                                ''' + kaiju_df_html + '''
                            </div>
                        </div>
                    </div>
                     <script type="text/javascript">
       
                        // function to get each tab details
                        const tabs = document.querySelectorAll('[data-tab-value]')
                        const tabInfos = document.querySelectorAll('[data-tab-info]')
                
                        tabs.forEach(tab => {
                            tab.addEventListener('click', () => {
                                const target = document
                                    .querySelector(tab.dataset.tabValue);
                                tabInfos.forEach(tabInfo => {
                                    tabInfo.classList.remove('active')
                                })
                                target.classList.add('active');
                            })
                        })
                    </script>
                </body>
            </html>'''
        report = open(sample_name + "_read_classification_report.html", "w")
        report.write(html_string)
        report.close()

    elif glob.glob("*_bracken_report_viral.txt") and not glob.glob("*_kaiju_summary_viral.tsv") and not glob.glob("*_viral_spp_abundance*.txt"):
        html_string = '''
        <html>
            <head>
                <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
                <style>body{ margin:0 100; background:whitesmoke; }</style>
            </head>
            <body>
                <h1>Direct read taxonomic classification (nucleotide-based)</h1>
                <!-- *** Section 1 *** --->
                <h2>Section 1: Kraken2 and Bracken viral results filtered </h2>
                <p>Only viral matches which represented >=0.0001 of the total read fraction are retained </p>
                ''' + bracken_df_filtered_html + '''

                <!-- *** Section 2 *** --->
                <h2>Section 2: All Kraken2 and Bracken viral results </h2>
                ''' + bracken_df_html + '''
            </body>
        </html>'''

        report = open(sample_name + "_read_classification_report.html", "w")
        report.write(html_string)
        report.close()

    elif glob.glob("*_kaiju_summary_viral.tsv") and not glob.glob("*_bracken_report_viral.txt") and not glob.glob("*_viral_spp_abundance*.txt"):
        html_string = '''
        <html>
            <head>
                <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
                <style>body{ margin:0 100; background:whitesmoke; }</style>
            </head>
            <body>
                <h1>Direct read taxonomic classification (protein-based)</h1>
                <!-- *** Section 1 *** --->
                <h2>Section 1: Kaiju results filtered </h2>
                <p>Only viral matches which represented >=0.05 of the total read fraction are retained </p>
                ''' + kaiju_df_filtered_html + '''

                <!-- *** Section 2 *** --->
                <h2>Section 2: All Kaiju results </h2>
                ''' + kaiju_df_html + '''
            </body>
        </html>'''

        report = open(sample_name + "_read_classification_report.html", "w")
        report.write(html_string)
        report.close()


    elif glob.glob("*_viral_spp_abundance*.txt") and not glob.glob("*_bracken_report_viral.txt") and not glob.glob("*_kaiju_summary_viral.tsv"):
        html_string = '''
            <html>
            <head>
                <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
                <style>body{ margin:0 100; background:whitesmoke; }</style>
            </head>
            <body>
                <h1>Direct read homology search (megablast)</h1>
               <!-- *** Section 1 *** --->
                <h2>Section 1: Number of reads matching to viral species filtered </h2>
                <p>Only blast viral matches which show >90% query coverage for NCBI and >95% query coverage for local viral database were considered here.</p>
                ''' + megablast_summary_per_spp_high_conf + '''
                
                <h2>Section 2: Number of reads matching to viral species </h2>
                <p>All blast viral matches were considered here.</p>
                ''' + megablast_summary_per_spp + '''
            </body>
        </html>'''

        report = open(sample_name + "_read_classification_report.html", "w")
        report.write(html_string)
        report.close()
    
    elif glob.glob("*_bracken_report_viral.txt") and not glob.glob("*_kaiju_summary_viral.tsv") and glob.glob("*_viral_spp_abundance*.txt"):
        html_string = '''
            <html>
                <head>
                    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
                    <style>
                    body { margin:0 100; background:whitesmoke; 
                    }
                    container {
                        border: 1px solid grey;
                        margin: 1rem;
                    }
                    [data-tab-info] {
                        display: none;
                    }
 
                    .active[data-tab-info] {
                        display: block;
                    }
 
                    .tab-content {
                        margin-top: 1rem;
                        padding-left: 1rem;
                        font-size: 20px;
                        font-family: sans-serif;
                        font-weight: bold;
                        color: rgb(0, 0, 0);
                    }
 
                    .tabs {
                        border-bottom: 1px solid grey;
                        background-color: rgb(170, 245, 144);
                        font-size: 25px;
                        color: rgb(0, 0, 0);
                        display: flex;
                        margin: 0;
                    }
 
                    .tabs span {
                        background: rgb(170, 245, 144);
                        padding: 10px;
                        border: 1px solid rgb(255, 255, 255);
                    }
            
                    .tabs span:hover {
                        background: rgb(55, 219, 46);
                        cursor: pointer;
                        color: black;
                    }
                    </style>
                </head>
                <body>
                    <!-- Body Container -->
                    <div class="container">
                    
                        <!-- Tabs Detail -->
                        <div class="tabs">
                            <span data-tab-value="#tab_1">Megablast</span>
                            <span data-tab-value="#tab_2">Kraken</span>
                        </div>

                        <!-- Tab content -->
                        <div class="tab-content">
                            <div class="tabs__tab active" id="tab_1" data-tab-info>
                                <h1>Direct read homology search (megablast)</h1>
                                <!-- *** Section 1 *** --->
                                <h2>Section 1: Number of reads matching to viral species filtered </h2>
                                <p>Only blast viral matches which show >90% query coverage for NCBI and >95% query coverage for local viral database were considered here.</p>
                                ''' + megablast_summary_per_spp_high_conf + '''
                                
                                <h2>Section 2: Number of reads matching to viral species </h2>
                                <p>All blast viral matches were considered here.</p>
                                ''' + megablast_summary_per_spp + '''
                            </div>
                            <div class="tabs__tab" id="tab_2" data-tab-info>
                                <h1>Direct read taxonomic classification (nucleotide-based)</h1>
                                <!-- *** Section 1 *** --->
                                <h2>Section 1: Kraken2 and Bracken viral results filtered </h2>
                                <p>Only viral matches which represented >=0.0001 of the total read fraction are retained </p>
                                ''' + bracken_df_filtered_html + '''
                                <!-- *** Section 2 *** --->
                                <h2>Section 2: All Kraken2 and Bracken viral results </h2>
                                ''' + bracken_df_html + '''
                            </div>
                        </div>
                    </div>
                     <script type="text/javascript">
       
                        // function to get each tab details
                        const tabs = document.querySelectorAll('[data-tab-value]')
                        const tabInfos = document.querySelectorAll('[data-tab-info]')
                
                        tabs.forEach(tab => {
                            tab.addEventListener('click', () => {
                                const target = document
                                    .querySelector(tab.dataset.tabValue);
                                tabInfos.forEach(tabInfo => {
                                    tabInfo.classList.remove('active')
                                })
                                target.classList.add('active');
                            })
                        })
                    </script>                    
                </body>
            </html>'''

        report = open(sample_name + "_read_classification_report.html", "w")
        report.write(html_string)
        report.close()

    elif not glob.glob("*_bracken_report_viral.txt") and glob.glob("*_kaiju_summary_viral.tsv") and glob.glob("*_viral_spp_abundance*.txt"):
        html_string = '''
            <html>
                <head>
                    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
                    <style>
                    body { margin:0 100; background:whitesmoke; 
                    }
                    container {
                        border: 1px solid grey;
                        margin: 1rem;
                    }
                    [data-tab-info] {
                        display: none;
                    }
 
                    .active[data-tab-info] {
                        display: block;
                    }
 
                    .tab-content {
                        margin-top: 1rem;
                        padding-left: 1rem;
                        font-size: 20px;
                        font-family: sans-serif;
                        font-weight: bold;
                        color: rgb(0, 0, 0);
                    }
 
                    .tabs {
                        border-bottom: 1px solid grey;
                        background-color: rgb(170, 245, 144);
                        font-size: 25px;
                        color: rgb(0, 0, 0);
                        display: flex;
                        margin: 0;
                    }
 
                    .tabs span {
                        background: rgb(170, 245, 144);
                        padding: 10px;
                        border: 1px solid rgb(255, 255, 255);
                    }
            
                    .tabs span:hover {
                        background: rgb(55, 219, 46);
                        cursor: pointer;
                        color: black;
                    }
                    </style>
                </head>
                <body>
                    <!-- Body Container -->
                    <div class="container">
                    
                        <!-- Tabs Detail -->
                        <div class="tabs">
                            <span data-tab-value="#tab_1">Megablast</span>
                            <span data-tab-value="#tab_2">Kaiju</span>
                        </div>

                        <!-- Tab content -->
                        <div class="tab-content">
                            <div class="tabs__tab active" id="tab_1" data-tab-info>
                                <h1>Direct read homology search (megablast)</h1>
                                <!-- *** Section 1 *** --->
                                <h2>Section 1: Number of reads matching to viral species filtered </h2>
                                <p>Only blast viral matches which show >90% query coverage for NCBI and >95% query coverage for local viral database were considered here.</p>
                                ''' + megablast_summary_per_spp_high_conf + '''
                                
                                <h2>Section 2: Number of reads matching to viral species </h2>
                                <p>All blast viral matches were considered here.</p>
                                ''' + megablast_summary_per_spp + '''
                            </div>
                            <div class="tabs__tab" id="tab_2" data-tab-info>
                                <h1>Direct read taxonomic classification (protein-based)</h1>
                                <!-- *** Section 1 *** --->
                                <h2>Section 1: Kaiju results filtered </h2>
                                <p>Only viral matches which represented >=0.05 of the total read fraction are retained </p>
                                ''' + kaiju_df_filtered_html + '''

                                <!-- *** Section 2 *** --->
                                <h2>Section 2: All Kaiju results </h2>
                                ''' + kaiju_df_html + '''
                            </div>
                        </div>
                    </div>
                     <script type="text/javascript">
       
                        // function to get each tab details
                        const tabs = document.querySelectorAll('[data-tab-value]')
                        const tabInfos = document.querySelectorAll('[data-tab-info]')
                
                        tabs.forEach(tab => {
                            tab.addEventListener('click', () => {
                                const target = document
                                    .querySelector(tab.dataset.tabValue);
                                tabInfos.forEach(tabInfo => {
                                    tabInfo.classList.remove('active')
                                })
                                target.classList.add('active');
                            })
                        })
                    </script>                    
                </body>
            </html>'''

        report = open(sample_name + "_read_classification_report.html", "w")
        report.write(html_string)
        report.close()

if __name__ == "__main__":
    main()  
