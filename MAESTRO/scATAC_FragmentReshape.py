# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-07-17 07:37:07
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-08-05 20:12:22


import argparse
import time, os
from collections import defaultdict

def CommandLineParser():
    parser=argparse.ArgumentParser(description = "Reshape the 10X fragment file to get tid of ")
    parser.add_argument("-F","--frag", dest = "fragfile", default = "",help = "The fragment file generated by 10X CellRanger ATAC, "
        "with barcodes and counts.")
    parser.add_argument("-O", "--output", dest = "outfile", default = "", help = "The output fragments file.")

    return parser.parse_args()


parser = CommandLineParser()
fragfile = parser.fragfile
outfile = parser.outfile


chr_list = ['chr'+str(i) for i in list(range(1,23))]
chr_list = chr_list + ['chrX', 'chrY']

start_time = time.time()
print("Start to reshape fragment file.",time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
with open(fragfile, "r") as frag_in, open(outfile, "w") as frag_out:
    for line in frag_in.readlines():
        items = line.strip().split("\t")
        if str(items[0]) in chr_list :
            frag_out.write("\t".join(items[0:]) + "\n")
end_time = time.time()
print("End:", end_time-start_time)
