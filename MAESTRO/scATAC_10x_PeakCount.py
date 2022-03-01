#! /usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Gali Bai
# @E-mail: gali.bai@hotmail.com
# @Date:   2020-08-30 22:26:54
# @Last Modified by:   Gali Bai
# @Last Modified time: 2021-08-30 20:46:06

from scipy import sparse
import os,sys
import time
import shutil
import argparse
import itertools
import multiprocessing as mp
import numpy as np
from collections import defaultdict
from functools import partial

from MAESTRO.scATAC_utility import *
from MAESTRO.scATAC_H5Process import *

def peakcount_parser(subparsers):
    """
    Add main function peakcount argument parsers.
    """

    workflow = subparsers.add_parser("scatac-peakcount",
        help = "Generate peak-cell count matrix. ")
    group_input = workflow.add_argument_group("Input arguments")
    group_input.add_argument("--filtered_count", dest = "filtered_count", default = " ", type = str, required = True,
        help = "Filtered Fragement file with peak count ")
    group_input.add_argument("--binary", dest = "binary", action = "store_true",
        help = "Whether or not to generate binary peak count matrix. If set, "
        "MAESTRO will binarize the peak-cell count matrix. "
        "If not (by default), MAESTRO will generate the original peak-count matrix. ")
    group_input.add_argument("--species", dest = "species", default = "GRCh38", required = False,
        choices = ["GRCh38", "GRCm38"], type = str,
        help = "Genome assembly of either 'GRCh38' for human or 'GRCm38' for mouse. DEFAULT: GRCh38.")

    group_output = workflow.add_argument_group("Output and running arguments")
    group_output.add_argument("--directory", dest = "directory", default = "MAESTRO",
        help = "Path to the directory where the result file shall be stored. DEFAULT: MAESTRO.")
    group_output.add_argument("--outprefix", dest = "outprefix", default = "10x-genomics",
        help = "Prefix of output files. DEFAULT: MAESTRO.")


def make_sparse_matrix(frag_file, binary, barcode_pos = 3, peak_pos = 0, count_pos = 4):
    barcodes = {}
    peaks = {}
    def get_index(value, encountered):
        if value in encountered:
            return encountered[value]
        else:
            i = len(encountered)
            encountered[value] = i
            return i
    data, rows, cols = [],[],[]
    for line in frag_file:
        line = line.strip().split('\t')
        barcode = line[barcode_pos]
        peak = "_".join(line[peak_pos: peak_pos + 3])
        i = get_index(barcode, barcodes)
        j = get_index(peak, peaks)
        count = line[count_pos]
        if binary:
            #print("binary")
            data.append(1)
        else:
            data.append(count)
        cols.append(i)
        rows.append(j)
    #print(data)
    matrix = sparse.coo_matrix((data, (rows,cols)), shape = (len(peaks),len(barcodes)))
    return matrix, list(barcodes.keys()), list(peaks.keys())


def peak_count_matrix(filtered_count, directory, outprefix, binary, species = 'GRCh38'):
    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass
    outfile = os.path.join(directory, outprefix + "_peak_count.h5")
    with open(filtered_count, 'r') as f:
        matrix, barcode_list, peak_list = make_sparse_matrix(f, binary)
        write_10X_h5(outfile, matrix, peak_list, barcode_list, genome = species, datatype = 'Peaks')
        return(matrix)
#peak_count_matrix(filtered_count, directory, outprefix, species, binary)
