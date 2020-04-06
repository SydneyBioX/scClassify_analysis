#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:47:47 2020

@author: yue
"""

                                
                                 
import os

trainlist =  [ "train_100", "train_200", "train_500", "train_10000", "train_2000", 
                 "train_5000", "train_10000", "train_20000", "train_30000",
                "test_100", "test_200", "test_500", "test_10000", "test_2000", 
                 "test_5000", "test_10000", "test_20000", "test_30000", 
                 "train_2class", "train_4class", "train_6class", 
                 "train_8class", "train_10class", "train_12class" ]

# change the current working directory 
os.chdir("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark")


#convert to hdf5 file 
#needed by ACTINN

dir = "/dora/nobackup/yuec/scclassify/benchmark/ACTINN/actinn_format.py"

dir_count = "/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/"

for i in range(0,len(trainlist)):
    testname = trainlist[i]
    os.system("python " + dir + " -i " + dir_count + testname  + ".csv -o "+ testname+" -f csv")
