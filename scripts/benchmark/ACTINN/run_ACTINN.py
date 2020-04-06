#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 14:59:37 2020

@author: yue
"""


import os 
import numpy as np
import pandas as pd
import time as tm

import tables
import tracemalloc

import time

def run_ACTINN(traindata, testdata, trainlabels, testlabels, n):

   

    dir = "/Users/yue/Dropbox (Sydney Uni)/PhD/scclassify/countmatrix/logcount/"
    
    truelab = []
    pred = []
    y_train=trainlabels
    y_test=testlabels
    y_train.to_csv("train_lab.csv", header = False, index = True, sep = '\t')
    y_test.to_csv("test_lab.csv", header = False, index = True, sep = '\t')
        
      
      

    
    dir = '"/Users/yue/Dropbox (Sydney Uni)/PhD/scclassify/ACTINN/actinn_predict.py"'
    
                               
    dir_data = '"/Users/yue/Dropbox (Sydney Uni)/PhD/scclassify/ACTINN/"'    
    

    #tracemalloc.start() 
    now = time.time()      
    os.system("python -m memory_profiler " + dir + " -trs "+ dir_data + traindata + ".h5 -trl " + dir_data + "train_lab.csv -ts " + dir_data + testdata+".h5")
    later = time.time()
    time_usage = int(later - now)
    
    
    data = pd.read_csv('/Users/yue/Dropbox (Sydney Uni)/PhD/scclassify/ACTINN/memory_profiler.log', 
                       error_bad_lines=False, skiprows=5, delim_whitespace=True)
    mem = data.iloc[:, 1]
    mem = np.array(mem)
    mem_index = np.argmax(mem)
    mem_usage = max(mem)
    mem_unit = data.iloc[mem_index, 2]
    
    
    #snapshot = tracemalloc.take_snapshot()
   

    ###############################################3
    truelab.extend(y_test.values)
    predlabels = pd.read_csv('/Users/yue/Dropbox (Sydney Uni)/PhD/scclassify/ACTINN/predicted_label.txt',header=0,index_col=None, sep='\t', usecols = [1])            
    pred.extend(predlabels.values)
    
            
    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)
   

    truelab.to_csv(traindata + "_" + testdata + "_ACTINN_True.csv", index = False)
    pred.to_csv(traindata + "_"+testdata + "_ACTINN_Pred.csv", index = False)
      
        
    return mem_usage, mem_unit , time_usage