#!/usr/bin/env python
#
# Copyright (C) 2019  Shengwei Hou, housw2010@gmail.com
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
from collections import Counter
from collections import OrderedDict
from ete3 import NCBITaxa
ncbi = NCBITaxa()


class Fasta(object):

    def __init__(self, header, seq):
        self.header = header
        self.seq = seq.upper()

    def __str__(self):
        return ">"+self.header+"\n"+self.seq+"\n"


class GLASSgo(Fasta):

    def __init__(self, header, seq):
        super().__init__(header, seq)

    @property
    def name(self):
        return self.header.split()[0]

    @property
    def gc_content(self):
        ret = 0
        if len(self.seq) == 0:
            return ret
        for c in self.seq:
            if c == "G" or c == "C":
                ret += 1
        return float(ret)/len(self.seq)

    @property
    def length(self):
        return len(self.seq)

    @property
    def taxid(self):
        taxid = ""
        if "taxID:" in self.header:
            taxid = self.header.strip().split("taxID:")[-1]
            # for taxid like this one: taxID:1639;1230340
            if ";" in taxid:
                taxid = taxid.split(";")[-1]
        if taxid:
            return str(taxid)
        else:
            return "NA"

    @property
    def lineage(self):
        if self.taxid != "NA":
            try:
                lineage = ncbi.get_lineage(self.taxid)
            except ValueError as e:
                lineage = ["NA"]
                print("WARNING: {taxid} was not found in current database, 'NA' returned for its lineage".format(taxid=self.taxid))
        else:
            lineage = ["NA"]
        return ";".join([str(item) for item in lineage])

    @property
    def taxonomy(self):
        if self.taxid != "NA":
            try:
                lineage = ncbi.get_lineage(self.taxid)
                taxid2name = ncbi.get_taxid_translator(lineage)
                taxonomy = [taxid2name.get(taxid, 'NA') for taxid in lineage]
            except ValueError as e:
                taxonomy = ["NA"]
                print("WARNING: {taxid} was not found in current database, 'NA' returned for its taxonomy".format(taxid=self.taxid))
        else:
            taxonomy = ["NA"]
        return ";".join([str(item) for item in taxonomy])

    @property
    def identity(self):
        identity = None
        if "VAL:" in self.header:
            identity = self.header.split("VAL:")[-1].split("%-taxID")[0]
            identity = float(identity)
        if identity:
            return identity
        else:
            return "NA"

    @property
    def accession(self):
        accession = ""
        if ":" in self.header:
            accession = self.header.split(":")[0]
        if accession:
            return accession
        else:
            return "NA"

    @property
    def position(self):
        start = end = strand = ''
        accession_position = self.header.split()[0]
        if ":" in accession_position:
            position = accession_position.split(':')[-1]
            if position.startswith('c'):
                position = position[1:]
                strand = '-'
                end, start = position.split('-')
            else:
                strand = '+'
                start, end = position.split('-')
        return (start, end, strand)
    
    @property    
    def start(self):
        if self.position[0]:
            return self.position[0]
        else:
            return "NA"
    
    @property    
    def end(self):
        if self.position[1]:
            return self.position[1]
        else:
            return "NA"

    @property    
    def strand(self):
        if self.position[2]:
            return self.position[2]
        else:
            return "NA"

    def to_tsv(self):
        attr = [self.name, self.taxid, self.lineage, self.taxonomy, self.accession, self.start, self.end, self.strand, self.identity, "{:d}".format(self.length), "{:.2f}".format(self.gc_content), self.seq]
        return "\t".join([ str(it) for it in attr])
        

def parse_GLASSgo_fasta(fasta_file):
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
                    yield GLASSgo(header, "".join(seq))
                header = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        yield GLASSgo(header, "".join(seq))
        


def main():

    # main parser
    parser = argparse.ArgumentParser(description="convert GLASSgo fasta into tabular output with additional taxa information")
    parser.add_argument("input_GLASSgo_file", help="input GLASSgo output file in fasta format")
    parser.add_argument("-p", "--prefix", help="output prefix")
    parser.add_argument("-o", "--out_folder", help="output directory [./]", default="./")
    parser.add_argument("-f", "--force", action="store_true", help="force to overwrite the output")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.0.1")

    # show usage
    if len(sys.argv) < 2:
        sys.stderr.write("\nError: Not enough parameters were provided, please refer to the usage.\n")
        sys.stderr.write(parser.format_help())
        sys.exit(1)

    # parse args
    args = parser.parse_args()    
        
    # input and output handling
    if not args.prefix:
        basename = os.path.basename(args.input_GLASSgo_file)
        args.prefix = os.path.splitext(basename)[0]
    output_file = os.path.join(args.out_folder, args.prefix+".tsv")
    if os.path.exists(output_file):
        if args.force:
            print("Warning: output file exists, will be overwritten!")
        else:
            sys.stderr.write("Error: output file detected, please remove it or use --force option to overwrite it")
    
    # convert to tsv
    with open(output_file, "w") as oh:
        oh.write("Name\tTaxID\tLineage\tTaxonomy\tAccession\tStart\tEnd\tStrand\tIdentity\tLength\tGC_Content\tSeq\n")
        for sRNA in parse_GLASSgo_fasta(args.input_GLASSgo_file):
            line = sRNA.to_tsv()
            oh.write(line + "\n")


if __name__ == "__main__":
    main()
