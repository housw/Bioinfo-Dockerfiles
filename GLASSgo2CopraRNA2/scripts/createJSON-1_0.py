#!/usr/bin/env python3

__author__ = 'steffen'
import argparse
import math
from ete3 import NCBITaxa
import sys
import os
import warnings
warnings.filterwarnings('ignore', '.*was translated into.*',)


# absolute path of this script
scriptPath = os.path.dirname(os.path.realpath(__file__))

# absolute path to local NCBI taxa DB file
NCBITaxaDbFile = scriptPath + "/taxa.sqlite"

def anlyse_input_file(in_file_new_line, in_file_glassgo):
    unique_ids = dict()
    counter = 0
    if args.in_taxid == "" and args.in_GLASSgo_result == "":
        print("Please specify your input!")
        exit()
    if args.in_taxid != "" and args.in_GLASSgo_result == "":
        handle = open(in_file_new_line, "r")
        for tax_id in handle:
            tax_id = tax_id.rstrip()
            if tax_id in unique_ids:
                unique_ids[tax_id] += 1
            else:
                unique_ids[tax_id] = 1
            counter += 1
        handle.close()
    if args.in_taxid == "" and args.in_GLASSgo_result != "":
        handle = open(in_file_glassgo, "r")
        for line in handle:
            if line.startswith(">"):
                tax_id_all = line.split("-taxID:")
                if len(tax_id_all) == 2:
                    tax_id = tax_id_all[1].split(";")[0]
                    if tax_id in unique_ids:
                        unique_ids[tax_id] += 1
                    else:
                        unique_ids[tax_id] = 1
                    counter += 1
        handle.close()
    if counter == 0:
        sys.stdout.write("*** ERROR: No tax-ids were found! Please check your input file! ***\n")
        exit()
    return unique_ids


def compute_taxid_paths(unique_tax_id_hash, ):
    #ncbi = NCBITaxa()
    path_output = ""
    ncbi = NCBITaxa(NCBITaxaDbFile)
    pathways = list()
    tax_name_ctr = dict()
    max_scalable_hits = 1000
    max_value = 40
    for tax_id in unique_tax_id_hash:
        # save mode; because the tax id can also be a not parsable string
        try:
            # get pathway (ete3 package) => "['root', 'bacteria', 'bac1']"
            global_scaling_val = unique_tax_id_hash[tax_id]
            lineage = ncbi.get_lineage(int(tax_id))

            # prepare output for CopraRNA
            path_output += str(ncbi.get_rank(lineage)) + "\n"
            path_output += str(lineage) + "\n\n"

            names = ncbi.get_taxid_translator(lineage)
            tmp_path = list()
            for tax_id2 in lineage:
                tax_name = str(tax_id2) + ":" + str(names[tax_id2])
                if tax_name in tax_name_ctr:
                    tax_name_ctr[tax_name][0] += global_scaling_val
                else:
                    tax_name_ctr[tax_name] = list()
                    tax_name_ctr[tax_name].append(global_scaling_val)
                    #tax_name_ctr[tax_name][0] += unique_tax_id_hash[tax_id]
                    tax_name_ctr[tax_name].append(0)
                    tax_name_ctr[tax_name].append(0)
                tmp_path.append(tax_name)
            # normalize node values
            for tax_name in tax_name_ctr:
                if (tax_name_ctr[tax_name][0]) <= max_scalable_hits:
                    tax_name_ctr[tax_name][1] = math.sqrt(float(tax_name_ctr[tax_name][0])) * 1.26
                    tax_name_ctr[tax_name][2] = "passed"
                else:
                    tax_name_ctr[tax_name][1] = max_value
                    tax_name_ctr[tax_name][2] = "failed"
            # append sub-pathway to pathways
            pathways.append(tmp_path)
        except ValueError:
            pass

    return pathways, tax_name_ctr, path_output


def build_tree_from_root_to_leaf(master_paths):
    masterTable = dict()
    max_depth = 0
    for tmp_list in master_paths:
        for i in range(0, len(tmp_list)):
            if i in masterTable:
                if not tmp_list[i] in masterTable[i]:
                    masterTable[i][tmp_list[i]] = list()
                if i+1 in masterTable:
                    if i+1 < len(tmp_list) and not tmp_list[i+1] in masterTable[i+1]:
                        if i+1 < len(tmp_list):
                            if len(masterTable[i][tmp_list[i]]) > 0:
                                masterTable[i][tmp_list[i]].append(tmp_list[i+1])
                            else:
                                 masterTable[i][tmp_list[i]].append(tmp_list[i+1])
                        else:
                            if len(masterTable[i][tmp_list[i]]) > 0:
                                masterTable[i][tmp_list[i]].append("leaf")
                            else:
                                masterTable[i][tmp_list[i]].append("leaf")
                    else:
                        if len(masterTable[i][tmp_list[i]]) <= 0:
                            masterTable[i][tmp_list[i]] = list()
                            masterTable[i][tmp_list[i]].append("leaf")
                else:
                    if len(masterTable[i][tmp_list[i]]) <= 0:
                        masterTable[i][tmp_list[i]] = list()
                        if i+1 < len(tmp_list):
                            masterTable[i][tmp_list[i]].append(tmp_list[i+1])
                        else:
                            masterTable[i][tmp_list[i]].append("leaf")
            else:
                masterTable[i] = dict()
                masterTable[i][tmp_list[i]] = list()
                if i+1 < len(tmp_list):
                    masterTable[i][tmp_list[i]].append(tmp_list[i+1])
                else:
                    masterTable[i][tmp_list[i]].append("leaf") #todo leaf problem
            if i+1 > max_depth:
                max_depth = i+1
    return masterTable, max_depth


def build_tree_from_leaf_to_root(master_paths):
    masterTable = dict()
    max_depth = 0
    for tmp_list in master_paths:
        for i in range(0, len(tmp_list)):
            if i in masterTable:
                if not tmp_list[i] in masterTable[i]:
                    masterTable[i][tmp_list[i]] = list()
                if i-1 in masterTable:
                    if i+1 < len(tmp_list) and not tmp_list[i-1] in masterTable[i-1]:
                        if i+1 < len(tmp_list):
                            if len(masterTable[i][tmp_list[i]]) > 0:
                                masterTable[i][tmp_list[i]].append(tmp_list[i-1])
                            else:
                                 masterTable[i][tmp_list[i]].append(tmp_list[i-1])
                        else:
                            if len(masterTable[i][tmp_list[i]]) > 0:
                                masterTable[i][tmp_list[i]].append(tmp_list[i-1])
                            else:
                                masterTable[i][tmp_list[i]].append(tmp_list[i-1])
                    else:
                        if len(masterTable[i][tmp_list[i]]) <= 0:
                            masterTable[i][tmp_list[i]] = list()
                            masterTable[i][tmp_list[i]].append(tmp_list[i-1])
                else:
                    if len(masterTable[i][tmp_list[i]]) <= 0:
                        masterTable[i][tmp_list[i]] = list()
                        if i+1 < len(tmp_list):
                            masterTable[i][tmp_list[i]].append(tmp_list[i-1])
                        else:
                            masterTable[i][tmp_list[i]].append(tmp_list[i-1])
            else:
                masterTable[i] = dict()
                masterTable[i][tmp_list[i]] = list()
                if i+1 < len(tmp_list):
                    masterTable[i][tmp_list[i]].append(tmp_list[i-1])
                else:
                    masterTable[i][tmp_list[i]].append(tmp_list[i-1])
            if i+1 > max_depth:
                max_depth = i+1
    return masterTable, max_depth



def setup_block(name, parent, value, leaf_or_inner, available_children, max_children, storage):
    block = ""
    type_c = "lightgrey"
    level_c = "lightgrey"
    type_c_error = "pink"

    if leaf_or_inner == "root":
        block += "[{"
        block += "\"name\": \"" + str(name) + "{" + str(value[0]) + "}\","
        block += "\"parent\": \"null\","
        block += "\"value\": \"" + str(value[1]) + "\","
        if value[2] == "passed":
            block += "\"type\": \"" + str(type_c) + "\","
        else:
            block += "\"type\": \"" + str(type_c_error) + "\","
        block += "\"level\": \"" + str(level_c) + "\","
        block += "\"children\": "
        return block
    else:
        # leaf_level == 0 means, that
        if available_children+1 == max_children:
            block += "["
        block += "{"
        block += "\"name\": \"" + str(name) + "{" + str(value[0]) + "}\","
        block += "\"parent\": \"" + str(parent) + "\","
        block += "\"value\": \"" + str(value[1]) + "\","
        if value[2] == "passed":
            block += "\"type\": \"" + str(type_c) + "\","
        else:
            block += "\"type\": \"" + str(type_c_error) + "\","

        if leaf_or_inner == "leaf_node":
            block += "\"level\": \"" + str(level_c) + "\""
            #block += "}]"
        else:
            block += "\"level\": \"" + str(level_c) + "\","
            block += "\"children\": "
            if name in storage:
                block += str(storage[name])
        if available_children > 0:
            block += "}"
            block += ","
        else:
            block += "}]"
        return block


def analyse_tree_from_leaf(masterTable, max_depth, masterTable_f, tax_name_ctr):
    json_output = dict()
    final_json_string = ""

    for i in range(0, max_depth):
        collection_hash = dict()
        depth = (max_depth-1) - i
        sub_tree = masterTable[depth]
        for node_name in sub_tree:
            parent = sub_tree[node_name][0]
            if node_name == "1:root":
                collection_hash["final"] = list()
                break

            if parent in collection_hash:
                collection_hash[parent].append(node_name)
            else:
                collection_hash[parent] = list()
                collection_hash[parent].append(node_name)

        # analyse substructure
        for parent in collection_hash:
            if parent != "final":
                child_list = sorted(collection_hash[parent])
                av_child = len(child_list)
                for child in child_list:
                    av_child -= 1
                    if masterTable_f[depth][child][0] == "leaf":
                        sub_json = setup_block(child, parent, tax_name_ctr[child], "leaf_node", av_child, len(child_list), json_output) # todo value should be changed
                    else:
                        sub_json = setup_block(child, parent, tax_name_ctr[child], "inner_node", av_child, len(child_list), json_output)
                    # remove ' from names
                    sub_json = sub_json.replace("'","")
                    if parent in json_output:
                        json_output[parent] += sub_json
                    else:
                        json_output[parent] = sub_json
            elif parent == "final":
                # finalize json output -> add root node
                sub_json = setup_block("1:root", "null", tax_name_ctr["1:root"], "root", -1, len(child_list), json_output)
                # remove ' from names
                sub_json = sub_json.replace("'","")
                final_json_string = str(sub_json)
                final_json_string += json_output["1:root"]
                final_json_string += "}]"
                break
    return final_json_string


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_taxid", help="Taxonomic input - one taxonomic id per line.", type=str, default="")
    parser.add_argument("-g", "--in_GLASSgo_result", help="Read GLASSgo result (FASTA).", type=str, default="")
    parser.add_argument("-o", "--out_file", help="Write JSON to user specified file.", type=str, default="")
    parser.add_argument("-p", "--out_paths", help="Write PATHS used for CopraRNA.", type=str, default="")
    parser.add_argument("-u", "--update_db", help="Install/Update taxonomic database. (Default=false [true,false])", type=str, default="false")
    args = parser.parse_args()

    # INIT DB or update DB
    if args.update_db == "true":
        #ncbi = NCBITaxa()
        ncbi = NCBITaxa(NCBITaxaDbFile)
        ncbi.update_taxonomy_database()
        exit()

    # (1) SETUP PATHS
    # read taxonomic ids and compute weights
    unique_ids = anlyse_input_file(args.in_taxid, args.in_GLASSgo_result)

    # read unique taxonomic ids and compute paths => "root.Bacteria.BacA"
    master_paths, tax_name_ctr, path_output = compute_taxid_paths(unique_ids)

    # (2) SETUP TREES
    f_tree, f_depth = build_tree_from_root_to_leaf(master_paths)
    r_tree, r_depth = build_tree_from_leaf_to_root(master_paths)

    # (3) START TREE ANALYSIS AND SETUP JSON FILE
    json_output = analyse_tree_from_leaf(r_tree, r_depth, f_tree, tax_name_ctr)

    # (4) PRINT OUTPUT TO STDOUT OR WRITE TO FILE - JSON FILE
    if args.out_file == "":
        print(json_output)
    else:
        handle = open(args.out_file, "w")
        handle.write(json_output)
        handle.close()

    # (5) PRINT OUTPUT TO STDOUT OR WRITE TO FILE - PATH FILE
    if args.out_paths == "":
        print(path_output)
    else:
        handle = open(args.out_paths, "w")
        handle.write(path_output)
        handle.close()
