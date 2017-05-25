#!/usr/bin/env python3
# Author: Junjie Zhu

from __future__ import print_function
from os import path
import numpy as np
import time
import sys

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def load_tissue_map(fname):
    id_map = {}
    with open(fname, 'r') as f:
        for line in f:
            fields = line.rstrip().split("\t")
            id_map[fields[1]] = fields[0]
    eprint("Loaded {} tissues".format(len(id_map)))
    return(id_map)
    
def load_sel_map(fname_pfx):
    id_map = {}
    for i in range(1,23):
        fname = fname_pfx + "_part{}.txt".format(i) 
        with open(fname, 'r') as f:
            for line in f:
                id_map[line.rstrip()] = i
        eprint(fname)
    return(id_map)

def load_conversion_map(fname, sel_map):
    id_map = {}
    header=True
    with open(fname, 'r') as f:
        for line in f:
            if header:
                header=False
                continue
            fields = line.rstrip().split("\t")
            if fields[6] in sel_map:
                id_map[fields[2]] = fields[6]
    eprint("Loaded ID selected conversion map with {} keys".format(len(id_map)))
    return(id_map)

def load_ntissues(fname):
    return(sum(1 for line in open(fname)))


def write_sel_index_set(chunk_path, sel_map):
    # load gene list
    with open(path.join(chunk_path, "Genes.txt"), "r") as g:
        g_list = [line.rstrip() for line in g]

    # iterate through gene list 
    for gene in g_list:
        v_fname = path.join(chunk_path, "Vars_{}.txt".format(gene))
        s_fname = path.join(chunk_path, "prunedset_{}.npy".format(gene))
        prune_list = []
        counter = 1 # index from 1 for R processing
        with open(v_fname, "r") as v: 
            for line in v:
                var = line.rstrip()
                if var in sel_map:
                    prune_list.append(counter)
                counter += 1
        print("{}\t{}\t{}".format(gene, len(prune_list), counter))
        np.save(s_fname, np.array(prune_list))
    



##########################
if __name__ == "__main__":

    main_path = "/share/PI/sabatti/controlled_access_gtex_data/our_analysis"
    tmp_path  = "/scratch/PI/sabatti/controlled_access_data/fastqtl_tmp" 

    tissue = "Muscle_Skeletal"
    tissue_map = load_tissue_map(path.join(main_path, "variant_pruning", "tissue_alt_names.txt"))

    tissue_fname_pfx = path.join(main_path,"variant_pruning","ld-subset", tissue_map[tissue])
    sel_set = load_sel_map(tissue_fname_pfx)

    lookup_fname = "GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_chr1-22+X_genot_imput_info04_maf01_HWEp1E6_variant_id_lookup.txt"
    sel_map = load_conversion_map(path.join(main_path, "genotype_data", lookup_fname), sel_set)

    # do multicore processing here for speed if necessary
    for chunk in range(1, 101):
    # for chunk in range(1, 2):
        start = time.time()
        ntissues = load_ntissues(path.join(main_path, "sample_list", "intersect_ids_{}.txt".format(tissue)))
        chunk_path = path.join(tmp_path, tissue, "{0}_{1}_chunk{2:03d}_mtx".format(tissue, ntissues, chunk))
        if (path.exists(chunk_path)): 
            eprint(chunk_path)
            write_sel_index_set(chunk_path, sel_map)
        else:
            eprint("Error: {} does not exists!".format(chunk_path))
            exit(1)
        eprint(time.time() - start)


