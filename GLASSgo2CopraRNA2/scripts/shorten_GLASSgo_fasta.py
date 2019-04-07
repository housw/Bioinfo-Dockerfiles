#!/usr/bin/env python


import os, sys
import argparse

class Fasta(object):

    def __init__(self, header, seq):
        self.header = header
        self.seq = seq.upper()

    def get_gc_content(self):
        ret = 0
        if len(self.seq) == 0:
            return ret
        for c in self.seq:
            if c == "G" or c == "C":
                ret += 1
        return float(ret)/len(self.seq)

    def __str__(self):
        return ">"+self.header+"\n"+self.seq+"\n"


def parse_fasta(fasta_file):
    """
    :param fasta_file: input fasta file
    :return:           yield Fasta record as a generator
    """
    header = ""
    seq = []
    with open(fasta_file, "r") as ih:
        for line in ih:
            if line.startswith(">"):
                if header:
                    yield Fasta(header, "".join(seq))
                header = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        yield Fasta(header, "".join(seq))


def shorten_fasta(input_fasta, output_file, maximum_records=400, identity_cutoff=0.3):
    header2record_and_id = {} # fasta_record : identity
    known_headers = set()
    with open(output_file, "w") as oh:
        for rec in parse_fasta(input_fasta):
            # GLASSgo hits
            if " " in rec.header:
                if ">" in rec.header:
                    rec.header = rec.header.split(">")[0]
                short_header, description = rec.header.strip().split(' ', 1)
                short_header = short_header.replace(".", "_").replace(":", "_").replace("-", "_")
                rec.header = short_header
                try:
                    identity = float(description.split("VAL:")[-1].split("%-taxID")[0])/100
                except Exception as e:
                    print("No identity was found for %s, use 1.0"%rec.header)
                    identity = 1.0
            # query
            else:
                short_header = rec.header
                identity = 1.0
            
            # collect hits with minimum identity cutoff
            if short_header not in header2record_and_id:
                if identity >= identity_cutoff and short_header not in header2record_and_id:
                    header2record_and_id.update({short_header:[rec, identity]})

        # sort by identity
        from collections import OrderedDict
        sorted_header2record_and_id = OrderedDict(sorted(header2record_and_id.items(), key=lambda item: item[1][1], reverse=True))

        # apply maximum_records cutoff
        for i, v in enumerate(sorted_header2record_and_id.values()):
            if i <= maximum_records:
                oh.write(str(v[0]))

def main():

    # main parser
    parser = argparse.ArgumentParser(description="reformat the header of fasta alignment")
    parser.add_argument("input_file", help="input alignment in fasta format generated")
    parser.add_argument("-m", "--maximum_records", type=int, default=400, help="maximum sRNAs will be kept for RNAalifold, default=400")
    parser.add_argument("-i", "--identity", type=float, default=0.3, help="minimum sequence identity of sRNAs will be kept, default=0.3")
    parser.add_argument("-p", "--prefix", help="output prefix")
    parser.add_argument("-o", "--out_folder", help="output directory, default=./", default="./")
    parser.add_argument("-f", "--force", action="store_true", help="force to overwrite the output")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    # show usage
    if len(sys.argv) < 2:
        sys.stderr.write("\nError: Not enough parameters were provided, please refer to the usage.\n")
        sys.stderr.write(parser.format_help())
        sys.exit(1)

    # parse args
    args = parser.parse_args()

    # input and output handling
    if args.prefix:
        prefix = args.prefix
        out_file = os.path.join(args.out_folder, prefix+"_shorten.fa")
    else:
        basename = os.path.basename(args.input_file)
        prefix = os.path.splitext(basename)[0]
        out_file = os.path.join(args.out_folder, prefix+"_shorten.fa")

    if os.path.exists(out_file):
        if args.force:
            print("Warning: output file exists, will be overwritten!")
        else:
            sys.stderr.write("Error: output file detected, please remove it or use --force option to overwrite it")


    # shorten fasta
    shorten_fasta(args.input_file, out_file, args.maximum_records, args.identity)

if __name__ == "__main__":
    main()
