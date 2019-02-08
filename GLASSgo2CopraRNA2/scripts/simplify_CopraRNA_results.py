#!/usr/bin/env python 


import sys, os 
import subprocess


def get_refseq2genbank(input_refseq2genbank_tsv):
    refseq2genbank = {} # {refseq_ID:genbank_ID}
    with open(input_refseq2genbank_tsv, "r") as ih:
        for line in ih:
            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            refseq = line[0].strip()
            genbank = line[1].strip()
            if refseq not in refseq2genbank:
                refseq2genbank[refseq] = genbank
            else:
                sys.stdout.write("RefSeq ID %s was found multiple times"%refseq)
    return refseq2genbank


def main():
    if len(sys.argv) !=3:
        print("\nUsage: python %s [CopraRNA prediction parent folder] [RefSeq_ID to Genebank_ID tsv] \n"%sys.argv[0])
        sys.exit(0)


    # get refseq2genbank
    refseq2genbank = get_refseq2genbank(sys.argv[2])

    out_dir = sys.argv[1].strip('/') + "_simplified"
    os.system("mkdir -p %s"%out_dir)

    fdr_pos = 0
    pvalue_pos = 1
    geneID_pos = 2
    annotation_pos = -3
    fdr_cutoff = 0.5

    for folder in os.listdir(sys.argv[1]):
        ID = folder
        #if os.path.exists(os.path.join(folder, "CopraRNA_result.csv")):
        CopraRNA_result = os.path.join(sys.argv[1], folder, "CopraRNA_result.csv")
        if os.path.exists(CopraRNA_result):
            if os.stat(CopraRNA_result).st_size > 0:
                out_file = os.path.join(out_dir, ID+"_CopraRNA_result.tab")
                out_file_simp = os.path.join(out_dir, ID+"_CopraRNA_result_simp.tab")
                out_file_locus = os.path.join(out_dir, ID+"_CopraRNA_result_simp_locustag.txt")
                with open(CopraRNA_result, "r") as ih, open(out_file, "w") as oh1, open(out_file_simp, "w") as oh2, \
                    open(out_file_locus, "w") as oh3:
                    oh1.write(ih.readline().replace(",", "\t"))
                    oh2.write("CopraRNA_FDR\tCopraRNA_PValue\tLocus_Tag\tGene_Name\tEnergy[kcal/mol]\tIntaRNA_PValue\tPosition_mRNA\tPosition_sRNA\tannotation\n")
                    for line in ih:
                        oh1.write(line.replace(",", "\t"))
                        line = line.strip().split(",")
                        if line[fdr_pos] == "NA":
                            continue
                        fdr = float(line[fdr_pos])
                        pvalue = float(line[pvalue_pos])
                        locus_tag, other = line[geneID_pos].split("(")
                        locus_tag = locus_tag.upper()
                        other = other.strip(")").split("|")
                        Gene_Name = refseq2genbank[locus_tag] if refseq2genbank.get(locus_tag, None) else other[0]
                        Energy = other[1]
                        IntaRNA_PValue = other[2]
                        Position_mRNA = other[3] +" - "+ other[4]
                        Position_sRNA = other[5] +" - "+ other[6]
                        annotation = line[annotation_pos]
                        if fdr <= fdr_cutoff:
                            oh2.write(str(fdr) +"\t"+ str(pvalue) +"\t"+ locus_tag +"\t"+ Gene_Name +"\t"+ Energy +"\t"+ IntaRNA_PValue +"\t"+ Position_mRNA +"\t"+ Position_sRNA +"\t"+ annotation +"\n")
                            oh3.write(locus_tag+"\n")
            else:
                print("The result file is empty for %s"%ID)
        else:
            print("No result file returned for %s"%ID)

if __name__ == "__main__":
    main()
