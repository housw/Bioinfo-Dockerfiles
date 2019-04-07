#!/usr/bin/env python


# Copyright (C) 2018  Shengwei Hou, housw2010@gmail.com
#
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


import sys
import os
import argparse
from ete3 import Tree
import numpy as np
from collections import Counter
from collections import OrderedDict


def taxpath2ete(taxpath_list, sep=';', root='root'):
    """ parse taxa from taxonomy path seperated by 'sep', return a ete tree object 
    """

    # initialize an ete tree, set the root name
    t = Tree()
    t.name = root
    tree = t
    for taxpath in taxpath_list:
        path = taxpath.strip().split(sep)
        for tax in path:
            #tax = tax.strip().split("__")[-1]
            #print(tax)
            try:
                # from root to leaf to query node,
                # if it works, then this node already in t
                tree = t&tax
                continue
            except Exception as e:
                # if cannot query, then add this node
                tree.add_child(name=tax)
                tree = t&tax
        # reset tree to root for each path 
        tree = t

    # get the root of entire tree
    tree_root = t.get_tree_root()
 
    return tree_root



def get_newick_tree_data(GLASSgo_tsv, output_nwk, output_count, ranks=['class', 'order', 'family', 'genus', 'species']):

    # {tax:{count=100, identity = [], length = [], gc_content=[], rank={'phylum':'Firmicutes', 'class':'Bacilli', ...} }}
    taxpath2features = {}
    
    with open(GLASSgo_tsv, "r") as ih:
        ih.readline()
        for line in ih:
            line = line.strip().split("\t")
            if line[8] == "NA":  # no identity
                continue
            ranktax_list = line[3].strip().split(";")  # rank__tax list
            if ranktax_list[0] == "NA":
                continue
            rank_list = []
            tax_list = []
            for rank_tax in ranktax_list:
                _rank, _tax = rank_tax.split("__")
                _tax = _tax.replace(" ", "_").replace("[", "").replace("]", "")
                rank_list.append(_rank)
                tax_list.append(_tax)
            rank_dict = {}
            for _rank, _tax in zip(rank_list, tax_list):
                #print(_rank +":"+ _tax)
                if _rank != "no rank":
                    rank_dict[_rank] = _tax
            taxpath = ";".join(tax_list) 
            identity = float(line[8])
            length = int(line[9])
            gc_content = float(line[10])
            if taxpath not in taxpath2features:
                taxpath2features[taxpath] = {'count':1, 
                                             'identity':[identity], 'length':[length], 
                                             'gc_content':[gc_content], 'rank':rank_dict}
            else:
                taxpath2features[taxpath]['count'] += 1
                taxpath2features[taxpath]['identity'].append(identity)
                taxpath2features[taxpath]['length'].append(length)
                taxpath2features[taxpath]['gc_content'].append(gc_content)

    # write to newick 
    ete_tree = taxpath2ete(taxpath2features.keys(), sep=';', root='root')
    ete_tree.write(format=1, outfile=output_nwk)

    # write values 
    with open(output_count, "w") as oh:
        oh.write("Name\tCount\tIdentity\tLength\tGC_Content\t" +"\t".join(ranks) +"\n")
        for taxpath, features in taxpath2features.items():
            name = taxpath.strip().split(';')[-1]
            line = name + "\t{ct}\t{id:.2f}\t{le:.0f}\t{gc:.2f}".format(
                   ct=features['count'], 
                   id=np.mean(features['identity']), 
                   le=np.mean(features['length']),
                   gc=np.mean(features['gc_content']))
            rank_dict = taxpath2features[taxpath]['rank']
            #print(rank_dict)
            for rank_name in ranks:
                rank_tax = rank_dict.get(rank_name, 'NA')
                line += '\t' + rank_tax
            oh.write(line +"\n")


def main():

    # main parser
    parser = argparse.ArgumentParser(description="convert GLASSgo tsv to newick tree with associated feature values")
    parser.add_argument("input_GLASSgo_tsv", help="input GLASSgo tsv file")
    parser.add_argument("-p", "--prefix", help="output prefix")
    parser.add_argument("-o", "--out_folder", help="output directory [./]", default="./")
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
    if not args.prefix:
        basename = os.path.basename(args.input_GLASSgo_tsv)
        args.prefix = os.path.splitext(basename)[0]

    output_nwk = os.path.join(args.out_folder, args.prefix+"_tree.nwk")
    output_count =  os.path.join(args.out_folder, args.prefix+"_tree.tsv")

    if os.path.exists(output_count):
        if args.force:
            print("Warning: output file exists, will be overwritten!")
        else:
            sys.stderr.write("Error: output file detected, please remove it or use --force option to overwrite it")

    # convert
    get_newick_tree_data(args.input_GLASSgo_tsv, output_nwk, output_count)

if __name__ == "__main__":
    main()
