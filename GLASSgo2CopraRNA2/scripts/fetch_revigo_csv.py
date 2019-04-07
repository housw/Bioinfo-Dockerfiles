#!/usr/bin/env python


from robobrowser import RoboBrowser
import re
import sys
import os
import argparse


def GOstats2Revigo(input_GOstats_tsv, pvalue_cutoff=0.05, fdr_cutoff=1.0, output_column=3):
    """ parse GOstats tsv file, select GO terms meet pvalue and fdr cutoff, 
        return selected GO terms for revigo visualization

        output_column can be 2 or 3, when it's 2, only write GOID and Pvalue; 
        when it's 3, write GOID, Pvalue and Count
    """
    ret_str = ""
    with open(input_GOstats_tsv, "r") as ih:
        ih.readline() # header
        for line in ih:
            line = line.strip().split("\t")
            GO_ID = line[0]
            Pvalue = line[1]
            FDR = line[7]
            Count = line[4]
            if float(Pvalue) <= pvalue_cutoff and float(FDR) <= fdr_cutoff:
                if output_column == 3:
                    ret_str += GO_ID +"\t"+ Pvalue +"\t"+ Count +"\n"
                else:
                    ret_str += GO_ID +"\t"+ Pvalue +"\n" 
    return ret_str


# modified from https://gist.github.com/SamDM/b7e8a13a5529c24291e293ee6ebe2366
def scrape_revigo_csv(input_GOstats_tsv, out_file, pvalue_cutoff=0.05, fdr_cutoff=1.0):
    """ 
    """
    oh = open(out_file, "w")
    
    # get input goterms from GOstats result
    goterms = GOstats2Revigo(input_GOstats_tsv, pvalue_cutoff=pvalue_cutoff, fdr_cutoff=fdr_cutoff, output_column=3)
    if goterms:
        br = RoboBrowser(parser="lxml")
        br.open("http://revigo.irb.hr/")

        form = br.get_form()
        #print(form)
        form["goList"].value = goterms

        br.submit_form(form)

        download_rsc_link = br.find("a", href=re.compile("toR.jsp"))
        br.follow_link(download_rsc_link)
        #r_code = br.response.content.decode("utf-8")
        #print(r_code)

        br.back()

        download_csv_link = br.find("a", href=re.compile("export.jsp"))
        br.follow_link(download_csv_link)
        csv_content = br.response.content.decode("utf-8")
        oh.write(csv_content)
    else:
        oh.write("term_ID,description,frequency,plot_X,plot_Y,plot_size,log10 p-value,userVal_2,uniqueness,dispensability,representative,eliminated")

    oh.close()

def main():

    # main parser
    parser = argparse.ArgumentParser(description="fetch revigo summarized GO terms from GOstats output tsv file")
    parser.add_argument("input_GOstats_tsv", help="input GOstats enrichment result in tsv format")
    parser.add_argument("-p", "--pvalue_cutoff", type=float, default=0.05, help="p-value cutoff, default=0.05")
    parser.add_argument("-f", "--fdr_cutoff", type=float, default=1.0, help="fdr cutoff, default=1.0")
    parser.add_argument("-o", "--out_folder", help="output directory, default=./", default="./")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    if len(sys.argv) < 2:
        sys.stderr.write("\nERROR: Not enough parameters were provided, please refer to the usage.\n\n")
        sys.stderr.write(parser.format_help())
        sys.exit(1)

    args = parser.parse_args()

    # input and output handeling
    basename = os.path.basename(args.input_GOstats_tsv)
    prefix = basename.split(".")[0]
    out_file = os.path.join(args.out_folder, prefix+"_revigo.csv")

    if os.path.exists(out_file):
        sys.stdout.write("\nWarning: output file exists, will be overwriten!\n")

    # revigo
    scrape_revigo_csv(args.input_GOstats_tsv, out_file, args.pvalue_cutoff, args.fdr_cutoff)


if __name__ == "__main__":
    main()
