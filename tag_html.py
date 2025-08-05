#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os 
import pandas as pd
import argparse

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    # important parameters 
    parser.add_argument("-sample", "--sample", help="sample prefix used to find fastq file")
    parser.add_argument("-outdir", "--outdir", default = ".", help="output directory name")
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    sample = args.sample
    output_dir = args.outdir
    # read GFP results 
    print("Read GFP results ...")
    gfp_sublibrary_summary_df = pd.read_excel(os.path.join(output_dir, sample + "GFP.xlsx"), sheet_name  = "summary")
    gfp_tag_gene_df = pd.read_excel(os.path.join(output_dir, sample + "GFP.xlsx"), sheet_name  = "flag_gene")
    gfp_tag_gene_df.columns = [c.strip("_") for c in gfp_tag_gene_df.columns]
    # gfp_tag_gene_df["index"] = gfp_tag_gene_df["barcode"].str.slice(0,15)
    # gfp_tag_gene_df = pd.merge(gfp_sublibrary_summary_df[["library", "index"]], gfp_tag_gene_df, how = "right", on = "index")
    # gfp_tag_gene_df.sort_values(by = ['library', "index"], inplace=True, ignore_index=True)
    # read ChiC summary 
    print("Read ChiC results ...")
    chic_sublibrary_summary = pd.read_table(os.path.join(output_dir, sample + ".pre.txt"), sep = "\t", header = 0)
    total_GFP_read = gfp_sublibrary_summary_df["total_reads"].sum()
    primer_yield = gfp_sublibrary_summary_df["primer_yield"].mean()
    barcode_yield = gfp_sublibrary_summary_df["barcode_yield"].mean()
    gfp_yiled = gfp_sublibrary_summary_df["final_yield"].mean()
    cds_total = gfp_sublibrary_summary_df["CDS_align"].sum()
    frame_total = gfp_sublibrary_summary_df["frame_count"].sum()
    target_frame = round(frame_total/cds_total*100, 1)
    gene_freq = gfp_tag_gene_df.groupby("symbol").size().reset_index()
    gene_freq.columns = ["symbol", "freq"]
    gene_freq.sort_values(by = ["freq"], ascending=False, inplace = True, ignore_index=True)
    gene1 = gene_freq.iloc[0,0]; gene2 = gene_freq.iloc[1,0]; gene3 = gene_freq.iloc[2,0]
    freq1 = gene_freq.iloc[0,1]; freq2 = gene_freq.iloc[1,1]; freq3 = gene_freq.iloc[2,1]
    # read ChiC read per cell with GFP data
    if os.path.exists(os.path.join(output_dir, sample + "_clean_barcode.stat")):
        GFP_ChiC_df = pd.read_table(os.path.join(output_dir, sample + "_clean_barcode.stat"), sep = "\t", header = 0)
    else:
        gfp_tag_gene_df
    ### 
    gfp_df = pd.DataFrame({"Total reads":[round(total_GFP_read)], 
                  "Valid primer (%)":[round(primer_yield,1)], 
                  "Valid barcodes (%)":[round(barcode_yield, 1)], 
                  "Valid reads (%)":[round(gfp_yiled, 1)]}).transpose()
    gfp_gene_df = pd.DataFrame({"Valid frame (%)" : [target_frame], 
                  "Total unique genes":[len(gene_freq)],
                  "Top target genes":[", ".join([gene1, gene2, gene3])],
                  "Top target freq":[", ".join([str(freq1), str(freq2), str(freq3)])]}).transpose()
    #print(gfp_df)
    # read chic summary 
    total_ChiC_read = sum(chic_sublibrary_summary["total_reads"])
    ChiC_barcode_yield = chic_sublibrary_summary["barcode_yield"].mean()
    from utilities.parse_log import markdup_log_parser
    chic_dedup_df = pd.read_table(os.path.join(output_dir, sample + ".dedup.stat"), sep = "\t", header = 0, index_col = 0)
    duplication_ratio = chic_dedup_df.loc["PERCENT_DUPLICATION"].item()
    chic_df = pd.DataFrame({"Total reads":[round(total_ChiC_read)], 
                            "Valid barcodes (%)":[round(ChiC_barcode_yield, 1)],
                            "Duplication (%)":[round(float(duplication_ratio)*100, 1)]}).transpose()
    #print(chic_df)
    # read chic png 
    chic_barcode_distribution_png = os.path.join(output_dir, sample + "_barcode.png")
    import base64
    with open(chic_barcode_distribution_png, "rb") as img_file:
        encoded_string = base64.b64encode(img_file.read()).decode('utf-8')
    img_html = f'<img src="data:image/png;base64,{encoded_string}" width="800"/>'
    # 
    with open(os.path.join(output_dir, sample + "_clean_barcode.png"), "rb") as img_file:
        encoded_string = base64.b64encode(img_file.read()).decode('utf-8')
    clean_img_html = f'<img src="data:image/png;base64,{encoded_string}" width="800"/>'
    # with open(os.path.join(output_dir, sample + "pair.png"), "rb") as img_file:
    #     encoded_string = base64.b64encode(img_file.read()).decode('utf-8')
    # pair_html = f'<img src="data:image/png;base64,{encoded_string}" width="800"/>'
    # read chic barcode stat 
    chic_barcode_stat = pd.read_table(os.path.join(output_dir, sample + "_barcode.stat"), sep = "\t", header = 0)
    chic_max = chic_barcode_stat["freq"].max(); chic_mean = chic_barcode_stat["freq"].mean(); chic_median = chic_barcode_stat["freq"].median()
    chic_img_stat = pd.DataFrame({"Total cell (no filtering)":[len(chic_barcode_stat)], 
                                  "Max reads per cell":[chic_max], 
                                  "Mean reads per cell":[chic_mean], 
                                  "Median reads per cell":[chic_median],
                                  }).transpose()
    chic_clean_img_stat = pd.DataFrame({"Total cells (filtered by GFP)":[len(GFP_ChiC_df)],
                                        "Valid cells (filtered by GFP & ChiC readcount)":[len(GFP_ChiC_df.query("label == 1"))],
                                        "Max reads per cell":[GFP_ChiC_df["ChiC_read_count"].max()], 
                                        "Mean reads per cell":[GFP_ChiC_df["ChiC_read_count"].mean()],
                                        "Median reads per cell":[GFP_ChiC_df["ChiC_read_count"].median()]}).transpose()
    #print(chic_img_stat)
    # tab pages
    gfp_summary = "<h3>GFP sequencing</h3>" + gfp_df.reset_index().to_html(index = False, header = False,border=0)
    gfp_gene_summary = "<h3>GFP gene</h3>" + gfp_gene_df.reset_index().to_html(index = False, header = False,border=0)
    chic_summary = "<h3>ChiC sequencing</h3>" + chic_df.reset_index().to_html(index = False, header = False,border=0)
    img_summary = "<h3>ChiC reads per cell</h3>" + img_html
    clean_img_summary = "<h3>ChiC reads per cell filtered by GFP</h3>" + clean_img_html
    chic_img_stat_html = chic_img_stat.reset_index().to_html(index = False, header = False,border=0)
    chic_clean_img_stat_html = chic_clean_img_stat.reset_index().to_html(index = False, header = False, border = 0)
    # Build buttons/tabs (label)
    buttons_html = f"""
    <button class="tablinks" onclick="openTab(event, 'tab1')" id="defaultOpen1">Summary</button>
    <button class="tablinks" onclick="openTab(event, 'tab2')">GFP sublibrary summary</button>
    <button class="tablinks" onclick="openTab(event, 'tab3')">GFP tag genes</button>
    <button class="tablinks" onclick="openTab(event, 'tab4')">GFP tag genes frequency</button>
    <button class="tablinks" onclick="openTab(event, 'tab5')">GFP tag gene & ChiC read count</button>
    """
    gfp_sublibrary_summary_df.to_csv(os.path.join(output_dir, sample + "_gfp.csv"), sep = "\t", header = True, index = False)
    download_file = os.path.join(sample + "_gfp.csv")
    download_gfp = f"""
    <a href="{download_file}" download>
    <button style="margin: 10px 0;">Download CSV</button>
    </a>
    """
    gfp_tag_gene_df.to_csv(os.path.join(output_dir, sample + "_gene.csv") ,sep = "\t", header = True, index = False)
    download_file = os.path.join(sample + "_gene.csv")
    download_gene = f"""
    <a href="{download_file}" download>
    <button style="margin: 10px 0;">Download CSV</button>
    </a>
    """
    # Build contents html (label content)
    contents_html = f"""
    <div id="tab1" class="tabcontent">
    {gfp_summary}
    {gfp_gene_summary}
    {chic_summary}
    {img_summary}
    {chic_img_stat_html}
    {clean_img_summary}
    {chic_clean_img_stat_html}
    <br><br></div>
    <div id="tab2" class="tabcontent">{download_gfp}{gfp_sublibrary_summary_df.to_html(index = False,border=0, table_id = "gfp")}</div>
    <div id="tab3" class="tabcontent">{download_gene}{gfp_tag_gene_df.to_html(index = False,border=0, table_id = "gene")}</div>
    <div id="tab4" class="tabcontent">{gene_freq.to_html(index = False,border=0, table_id = "freq")}</div>
    <div id="tab5" class="tabcontent">{GFP_ChiC_df.to_html(index = False, border = 0, table_id="combined")}</div>
    """
    # html template 
    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
    <meta charset="UTF-8" />
    <title>{sample}</title>
    <!-- DataTables CSS -->
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css" />

    <!-- jQuery -->
    <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>

    <!-- DataTables JS -->
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

    <style>
        .tab {{
        overflow: hidden;
        border-bottom: 1px solid #ccc;
        }}
        .tab button {{
        background-color: inherit;
        border: none;
        outline: none;
        padding: 10px 20px;
        cursor: pointer;
        transition: background-color 0.3s;
        font-size: 16px;
        }}
        .tab button:hover {{
        background-color: #ddd;
        }}
        .tab button.active {{
        background-color: #ccc;
        }}
        .tabcontent {{
        display: none;
        padding: 20px;
        border: 1px solid #ccc;
        border-top: none;
        }}
        table {{
        border-collapse: collapse;
        width: 50%;
        }}
        th, td {{
        border-bottom: 1px solid #ccc;  /* Row lines */
        padding: 8px;
        text-align: left;
        }}
        tr:last-child td {{
        border-bottom: none;  /* Optional: remove bottom border on last row */
        }}
    </style>
    </head>
    <body>

    <h2>Summary report for {sample}</h2>

    <div class="tab">
    {buttons_html}
    </div>

    {contents_html}
    
    <script>
    function openTab(evt, tabName) {{
    var i, tabcontent, tablinks;
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {{
        tabcontent[i].style.display = "none";
    }}
    tablinks = document.getElementsByClassName("tablinks");
    for (i = 0; i < tablinks.length; i++) {{
        tablinks[i].className = tablinks[i].className.replace(" active", "");
    }}
    document.getElementById(tabName).style.display = "block";
    evt.currentTarget.className += " active";
    if (tabName === 'tab2') {{
        if (!$.fn.DataTable.isDataTable('#gfp')) {{
        $('#gfp').DataTable();
        }}
    }}
    if (tabName === 'tab3') {{
        if (!$.fn.DataTable.isDataTable('#gene')) {{
        $('#gene').DataTable();
        }}
    }}
    if (tabName === 'tab4') {{
        if (!$.fn.DataTable.isDataTable('#freq')) {{
        $('#freq').DataTable();
        }}
    }}
    if (tabName === 'tab5') {{
        if (!$.fn.DataTable.isDataTable('#combined')) {{
        $('#combined').DataTable();
        }}
    }}
    }}
    
    // Open the first tab by default
    document.getElementById("defaultOpen1").click();
    </script>

    </body>
    </html>
    """
    # write into html file 
    with open(os.path.join(output_dir, sample + ".html"), "w") as f:
        f.write(html_template)

if __name__ == "__main__":
    main()

