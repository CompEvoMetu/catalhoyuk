#!/bin/bash -l
for i in {1..40..1}; do
        touch ./haproh_027x_$i.py
        cat > ./haproh_027x_$i.py << EOF
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os as os
import sys as sys
import multiprocessing as mp

from hapsburg.PackagesSupport.hapsburg_run import hapsb_ind


iids=["cch289_merged","cch328_merged","cch251_merged","cch170_merged","cch294_merged","cch376_merged",
"cch144_merged","cch186_merged","cch312_merged","cch311_merged","cch151_merged","cch153_merged",
"cch245_merged","cch174_merged","cch510_merged","cch285_merged","cch416_merged","cch3418_merge","cch3208","cch205_merged"]

for iid in iids:
    hapsb_ind(iid=iid, chs=range(1,23), processes=15,
              path_targets='/mnt/NEOGENE1/projects/catalhoyuk_2020/roh/reich_roh/Final_Version_ROH/multitrial_ROH/eigen_file/all_bamlist.$i',
              h5_path1000g='/mnt/NEOGENE1/projects/catalhoyuk_2020/1240K/1000g1240khdf5/all1240/chr',
              meta_path_ref='/mnt/NEOGENE1/projects/catalhoyuk_2020/roh/elifnaz_roh/meta_df_all.csv',
              folder_out='/mnt/NEOGENE1/projects/catalhoyuk_2020/roh/reich_roh/Final_Version_ROH/multitrial_ROH/ROH_files/ROH_Result_Trial_$i',
              e_model="haploid", p_model="EigenstratUnpacked", n_ref=2504,
              random_allele=True, readcounts=False,
              delete=False, logfile=True, combine=True)
# post-processing
from hapsburg.PackagesSupport.pp_individual_roh_csvs import pp_individual_roh
# hapROHanalysis.py file end

EOF
done
