#!/usr/bin/env python

# Copyright (C) 2018  Shengwei Hou
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import print_function
import os
import sys
import argparse


def get_locustags(input_CopraRNA_result, out_file, pvalue_cutoff=1, fdr_cutoff=1, deliminator="\t"):
    """ get locus_tags from the top targets of CopraRNA result
    """
    with open(input_CopraRNA_result, "r") as ih, open(out_file, "w") as oh:
        ih.readline()
        for line in ih:
            line = line.strip().split(deliminator)
            fdr = float(line[0])
            pvalue = float(line[1])
            locus_tag = line[2].strip().split('(')[0]
            if fdr <= fdr_cutoff and pvalue <= pvalue_cutoff:
                oh.write(locus_tag.upper() +"\n")



def main():

    # main parser
    parser = argparse.ArgumentParser(description="get protein locus_tags from CopraRNA top targets")
    parser.add_argument("input_CopraRNA_result", help="input CopraRNA prediction result in tsv format")
    parser.add_argument("-p", "--pvalue_cutoff", type=float, default=1.0, help="p-value cutoff, default=1.0")
    parser.add_argument("-f", "--fdr_cutoff", type=float, default=1.0, help="fdr cutoff, default=1.0")
    parser.add_argument("-o", "--out_folder", help="output directory, default=./", default="./")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    if len(sys.argv) < 2:
        sys.stderr.write("\nERROR: Not enough parameters were provided, please refer to the usage.\n\n")
        sys.stderr.write(parser.format_help())
        sys.exit(1)

    args = parser.parse_args()

    # input and output handeling
    basename = os.path.basename(args.input_CopraRNA_result)
    prefix = basename.split(".")[0]
    out_file = os.path.join(args.out_folder, prefix+"_locustags.tsv")

    if os.path.exists(out_file):
        sys.stdout.write("Warning: output file exists, will be overwriten!")

    # convert
    get_locustags(args.input_CopraRNA_result, out_file, args.pvalue_cutoff, args.fdr_cutoff)


if __name__ == "__main__":
    main()
