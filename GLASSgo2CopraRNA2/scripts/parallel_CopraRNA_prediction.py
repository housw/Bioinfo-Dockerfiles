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
import subprocess
import argparse
import time
import logging
from multiprocessing import Pool


def prepare_input_sRNAs(input_parent_folder, out_folder="CopraRNA_Results", suffix="fasta", prefix=None):
    """
    """
    subprocess.call(['mkdir', '-p', out_folder])
    for sRNA_folder in os.listdir(input_parent_folder):
        path_to_sRNA_folder = os.path.join(input_parent_folder,sRNA_folder)
        if os.path.isdir(path_to_sRNA_folder):
            path_to_sRNA_output_folder = os.path.join(out_folder, sRNA_folder)
            subprocess.call(['mkdir', '-p', path_to_sRNA_output_folder])

            # only one sRNA in each sRNA_folder
            for sRNA in os.listdir(path_to_sRNA_folder):
                path_to_sRNA = os.path.join(path_to_sRNA_folder, sRNA)
                if not prefix:
                    if sRNA.endswith(suffix):
                        subprocess.call(['cp', path_to_sRNA, os.path.join(path_to_sRNA_output_folder, "sRNA.fa")])
                        break
                else:
                    if sRNA.startswith(prefix) and sRNA.endswith(suffix):
                        subprocess.call(['cp', path_to_sRNA, os.path.join(path_to_sRNA_output_folder, "sRNA.fa")])
                        break


# work horse for each CopraRNA run
def run_CopraRNA(CopraRNA_cmd, path_to_sRNA_folder):
    print("Running CopraRNA for ", path_to_sRNA_folder)
    try:
        os.chdir(path_to_sRNA_folder)
        os.listdir(path_to_sRNA_folder)
        ret_code = subprocess.call(CopraRNA_cmd, shell=False)
    except FileNotFoundError as e:
        print("Folder was not found", path_to_sRNA_folder)
        print(e)
        ret_code = 1
    return ret_code

# wrapp into one argument
def run_CopraRNA_Wrapper(args):
    return run_CopraRNA(*args)

# multiprocessing
def run_CopraRNA_in_parallel(batch=20, cores=2, out_folder="CopraRNA_Results", 
                             ntup=200, ntdown=100, region="5utr", enrich=200, topcount=200, websrv=False, noclean=False):
    """ parallel wrapper of CopraRNA
    """
    # command to run CopraRNA
    CopraRNA_cmd = ['CopraRNA2.pl', '-srnaseq', 'sRNA.fa',
                    '-ntup', str(ntup), '-ntdown', str(ntdown), 
                    '-region', region, '-enrich', str(enrich), 
                    '-topcount', str(topcount), '-cores', str(cores)]
    if websrv:
        CopraRNA_cmd.append("-websrv")
    if noclean:
        CopraRNA_cmd.append("-noclean")
    #CopraRNA_cmd = ['sleep', '20'] 

    #logging.basicConfig(format="%(asctime)-15s %(message)s", datefmt="%F %T",
    #                    level=logging.INFO)

    # create a batch pool
    pool = Pool(processes=batch)
    
    # list of (cmd, path_to_sRNA_folder)
    jobs = []
    for sRNA_folder in os.listdir(out_folder):
        path_to_sRNA_folder = os.path.join(os.path.abspath(sys.path[0]), out_folder, sRNA_folder)
        jobs.append((CopraRNA_cmd, path_to_sRNA_folder))
    
    for ret_code in pool.map(run_CopraRNA_Wrapper, jobs):
        print ("ret code is ", ret_code)


def main():

    # main parser
    parser = argparse.ArgumentParser(description="Wrapper of running CopraRNA_v2 in parallel")
    parser.add_argument("input_parent_folder", type=str, help="input parent folder contains sRNA subfolders,\
                        within each subfolder, there should be one input sRNA file in fasta format")
    parser.add_argument("-s", "--suffix", type=str, help="suffix for sRNA file, default is fasta", 
                        default="fasta")
    parser.add_argument("-p", "--prefix", type=str, help="prefix for sRNA file")
    parser.add_argument("-c", "--cores", type=int, help="number of cores will be used for each CopraRNA run, default=2", 
                        default=2)
    parser.add_argument("-b", "--batch", type=int, help="number of parallel CopraRNA runs each batch, default=20", 
                        default=20)
    parser.add_argument("-o", "--out_folder", type=str, help="output directory, default=CopraRNA_Results", default="CopraRNA_Results")
    parser.add_argument("--ntup", type=int, help="-ntup parameter, upstream length, default=200", default=200)
    parser.add_argument("--ntdown", type=int, help="-ntdown parameter, downstream length, default=200", default=100)
    parser.add_argument("--region", type=str, help="-region parameter, default=5utr", default="5utr")
    parser.add_argument("--enrich", type=int, help="-enrich parameter, default=200", default=200)
    parser.add_argument("--websrv", action='store_true', help=" switch to provide webserver output files")
    parser.add_argument("--noclean", action='store_true', help=" switch to prevent removal of temporary files")
    parser.add_argument("--topcount", type=int, help="-topcount parameter, default=200", default=200) 
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    if len(sys.argv) < 2:
        print("\nERROR: Not enough parameters were provided, please refer to the usage.\n", file=sys.stderr)
        print(parser.format_help(), file=sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # total number of sRNAs
    prepare_input_sRNAs(args.input_parent_folder, args.out_folder, args.suffix, args.prefix)

    # run CopraRNA
    run_CopraRNA_in_parallel(batch=args.batch, cores=args.cores, 
                             out_folder=args.out_folder, ntup=args.ntup, ntdown=args.ntdown, 
                             region=args.region, enrich=args.enrich, topcount=args.topcount, websrv=args.websrv, noclean=args.noclean)



if __name__ == "__main__":
    main()
