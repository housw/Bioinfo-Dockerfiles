#!/usr/bin/env python


import os
import sys
import textwrap
import glob
from Bio import SeqIO


def get_GLASSgo_html(GLASSgo_folder, ID):
    return GLASSgo_folder +"/"+ "GLASSgo_output_" +ID+ ".html"


def summarize_GLASSgo_results(GLASSgo_folder, ID):
    GLASSgo_result = GLASSgo_folder +"/"+ "GLASSgo_output_" +ID+ ".fa"
    hits = SeqIO.parse(GLASSgo_result, "fasta")
    Nr_of_hits = 0 
    Identities = []
    # passover the query
    next(hits)
    for hit in hits:
        Nr_of_hits += 1
        identity = float(hit.description.split("VAL:")[-1].split("%", 1)[0])
        Identities.append(identity)
    mean_identity = sum(Identities) / Nr_of_hits
    ret = "hits:%d, mean_ident:%.2f"%(Nr_of_hits, mean_identity)
    return ret

def get_GLASSgo2CopraRNA_results(GLASSgo2CopraRNA_folder, ID, file_type="fasta"):
    # file_type can be "fasta" or "pdf"
    file_type = file_type
    result_folder = GLASSgo2CopraRNA_folder +"/"+ ID
    for f in os.listdir(result_folder):
        if f.endswith(file_type):
            return result_folder +"/"+ f    


def get_simplified_CopraRNA2_tab(CopraRNA2_simp_folder, ID):
    return CopraRNA2_simp_folder +"/"+ ID+"_CopraRNA_result_simp.tab"


def get_Enrichment_results(CopraRNA2_simp_folder, ID, suffix):
    # suffix can be _BP.html, _CC.html, _MF.html, _GO_enrichment.pdf 
    for f in os.listdir(CopraRNA2_simp_folder):
        if f.startswith(ID) and f.endswith(suffix):
            return CopraRNA2_simp_folder +"/"+ f

def get_alignment_results(Alignment_folder, ID, file_type="fasta"):
    # file_type can be "fasta" or "html" or "nwk"    
    file_type = file_type
    result_folder = Alignment_folder
    for f in os.listdir(result_folder):
        if f.startswith(ID) and f.endswith(file_type):
            return result_folder +"/"+ f 

def get_RNAalifold_results(RNAalifold_folder, ID, suffix="_aln.pdf"):
    # suffix can be "_aln.pdf" or "_alirna.pdf"
    result_folder = RNAalifold_folder
    for f in os.listdir(result_folder):
        if f.startswith(ID) and f.endswith(suffix):
            return result_folder +"/"+ f

def get_Rscape_results(Rscape_folder, ID, suffix="_output.txt.sorted"):
    # suffix can be "_output.txt.sorted" or ".R2R.sto.pdf"
    result_folder = Rscape_folder
    for f in os.listdir(result_folder):
        if f.startswith(ID) and f.endswith(suffix):
            return result_folder +"/"+ f

def get_RNAcode_results(RNAcode_folder, ID):
    for folder in os.listdir(RNAcode_folder):
        if folder.startswith(ID):
            return """<a href="{RNAcode_results}">detected</a>""".format(RNAcode_results = RNAcode_folder+"/"+folder)
    else:
        return None
    

def get_Synteny_pdf(Synteny_folder, ID, suffix="_synteny.pdf"):
    result_folder = Synteny_folder
    for f in os.listdir(result_folder):
        if f.startswith(ID) and f.endswith(suffix):
            return result_folder +"/"+ f

def get_MEME_html(Motif_folder, ID):
    return Motif_folder +"/"+ ID + "/MEME_results/meme.html"

def main():

    # global parameters 
    GLASSgo_folder = "01_GLASSgo_Results"
    GLASSgo2CopraRNA_folder = "02_GLASSgo2CopraRNA"
    CopraRNA2_simp_folder = "03_CopraRNA_Results_v0_simplified"
    Alignment_folder = "04_Mafft_Results"
    RNAalifold_folder = "05_RNAalifold_Results"
    Rscape_folder = "06_Rscape_Results"
    RNAcode_folder = "07_RNAcode_Results"
    Synteny_folder = "08_Synteny_Results"
    Motif_folder = "09_Promoter_Results"

    # parse fasta to get all sRNA IDs
    sRNA_IDs = []
    sRNAs = SeqIO.parse("Sa/Sa_HG001_sRNAs.fa", "fasta")
    for sRNA in sRNAs:
        header = sRNA.name
        sRNA_IDs.append(header)

    # write html file
    with open("summarized_information.html", "w") as oh:
        # write metadata
        metadata = textwrap.dedent(
        """
        <html>
            <head>
                <title>Staphylococcus aureus ncRNAs</title>
            </head>
            <body>
                <table border="1">
        """)
        oh.write(metadata+"\n")

        # write header
        table_header = "<tr><th>ID</th><th>GLASSgo</th><th>GLASSgo2CopraRNA</th><th>CopraRNA2</th><th>Enrichment</th><th>Alignment</th><th>RNAalifold</th><th>Rscape</th><th>RNAcode</th><th>Synteny</th><th>PromoterMotif</th></tr>"
        oh.write(table_header+"\n")

        # write content
        for ID in sRNA_IDs:
            row = textwrap.dedent(
            """
            <tr>
                    <td>{sRNA_ID}</td>
                    <td><a href="{GLASSgo_html}">{GLASSgo_summary}</a></td>
                    <td><a href="{GLASSgo2CopraRNA_fasta}">fasta</a>, <a href="{GLASSgo2CopraRNA_pdf}">pdf</a></td>
                    <td><a href="{CopraRNA2_tab}">simplified_tab</a></td>
                    <td><a href="{BP_html}">BP</a>, <a href="{CC_html}">CC</a>, <a href="{MF_html}">MF</a>, <a href="{revigo_pdf}">revigo</a></td>     
                    <td><a href="{mafft_fasta}">fasta</a>, <a href="{mview_html}">html</a>, <a href="{fasttree_nwk}">tree</a></td>
                    <td><a href="{RNAalifold_aln}">aln</a>, <a href="{RNAalifold_alirna}">alirna</a></td>
                    <td><a href="{Rscape_pdf}">pdf</a>, <a href="{Rscape_txt}">txt</a></td>
                    <td>{RNAcode_results_html_code}</td>
                    <td><a href="{Synteny_pdf}">pdf</a></td>
                    <td><a href="{meme_html}">html</a></td>                      
            </tr>
            """.format(
                sRNA_ID = ID, 
                GLASSgo_html= get_GLASSgo_html(GLASSgo_folder, ID),
                GLASSgo_summary= summarize_GLASSgo_results(GLASSgo_folder, ID), 
                GLASSgo2CopraRNA_fasta = get_GLASSgo2CopraRNA_results(GLASSgo2CopraRNA_folder, ID, file_type="fasta"),
                GLASSgo2CopraRNA_pdf = get_GLASSgo2CopraRNA_results(GLASSgo2CopraRNA_folder, ID, file_type="pdf"),
                CopraRNA2_tab = get_simplified_CopraRNA2_tab(CopraRNA2_simp_folder, ID),
                BP_html = get_Enrichment_results(CopraRNA2_simp_folder, ID, suffix="_BP.html"),
                CC_html = get_Enrichment_results(CopraRNA2_simp_folder, ID, suffix="_CC.html"),
                MF_html = get_Enrichment_results(CopraRNA2_simp_folder, ID, suffix="_MF.html"),
                revigo_pdf = get_Enrichment_results(CopraRNA2_simp_folder, ID, suffix="_GO_enrichment.pdf"),
                mafft_fasta = get_alignment_results(Alignment_folder, ID, file_type="fasta"),
                mview_html = get_alignment_results(Alignment_folder, ID, file_type="html"),
                fasttree_nwk = get_alignment_results(Alignment_folder, ID, file_type="nwk"),
                RNAalifold_aln = get_RNAalifold_results(RNAalifold_folder, ID, suffix="_aln.pdf"),
                RNAalifold_alirna = get_RNAalifold_results(RNAalifold_folder, ID, suffix="_alirna.pdf"),
                Rscape_pdf = get_Rscape_results(Rscape_folder, ID, suffix=".R2R.sto.pdf"),
                Rscape_txt = get_Rscape_results(Rscape_folder, ID, suffix="_output.txt.sorted"),
                RNAcode_results_html_code = get_RNAcode_results(RNAcode_folder, ID),
                Synteny_pdf = get_Synteny_pdf(Synteny_folder, ID, suffix="_synteny.pdf"),
                meme_html = get_MEME_html(Motif_folder, ID)
            ))

            oh.write(row +"\n")

        #get_Enrichment_results(CopraRNA2_simp_folder, ID),

        html_end = textwrap.dedent(
        """
                </table>
            </body>
        </html>
        """)
        oh.write(html_end+"\n")

if __name__ == "__main__":
    main()
