#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 18:36:23 2020

@author: yuec
"""

                                
                                 
import os



import os 
import numpy as np
import pandas as pd


import time

def run_ACTINN(traindata, testdata, n):

   

    dir = "/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/"
    trainlabels =  pd.read_csv(dir + traindata + "_label.csv", header=0,index_col=0, sep=',')
    testlabels =  pd.read_csv(dir + testdata + "_label.csv", header=0,index_col=0, sep=',')
    truelab = []
    pred = []
    y_train=trainlabels
    y_test=testlabels
    
    
    os.chdir("/dora/nobackup/yuec/scclassify/benchmark/ACTINN/vary_celltype")
    y_train.to_csv("train_lab.csv", header = False, index = True, sep = '\t')
    y_test.to_csv("test_lab.csv", header = False, index = True, sep = '\t')
        
      
    
    dir = '"/dora/nobackup/yuec/scclassify/benchmark/ACTINN/actinn_predict.py"'
    
                               
    dir_data = '"/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/"'    
    dir_label = '"/dora/nobackup/yuec/scclassify/benchmark/ACTINN/vary_celltype/"'    

    #tracemalloc.start() 
    now = time.time()      
    os.system("python -m memory_profiler " + dir + " -trs "+ dir_data + traindata + ".h5 -trl " + dir_label + "train_lab.csv -ts " + dir_data + testdata+".h5")
    later = time.time()
    time_usage = int(later - now)
    
    
    data = pd.read_csv('/dora/nobackup/yuec/scclassify/benchmark/ACTINN/vary_celltype/memory_profiler.log', 
                       error_bad_lines=False, skiprows=5, delim_whitespace=True)
    mem = data.iloc[:, 1]
    mem = np.array(mem)
    mem_index = np.argmax(mem)
    mem_usage = max(mem)
    mem_unit = data.iloc[mem_index, 2]
    
    
    #snapshot = tracemalloc.take_snapshot()
   

    ###############################################3
    truelab.extend(y_test.values)
    predlabels = pd.read_csv('/dora/nobackup/yuec/scclassify/benchmark/ACTINN/vary_celltype/predicted_label.txt',header=0,index_col=None, sep='\t', usecols = [1])            
    pred.extend(predlabels.values)
    
            
    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)
   

    truelab.to_csv(n + "_ACTINN_True.csv", index = False)
    pred.to_csv(n + "_ACTINN_Pred.csv", index = False)
      
        
    return mem_usage, mem_unit , time_usage









trainlist =  [ "train_2class", "train_4class", "train_6class", "train_8class", "train_10class", "train_12class" ]


df = pd.DataFrame()

for i in range(0, len(trainlist)):
  train = trainlist[i]
  print(train)
  mem_usage, mem_unit , time_usage = run_ACTINN( train , "test_5000", train  )
  
  if df.empty :
    df = df.append({'n': train, 'mem': mem_usage, 'unit': mem_unit, 'time': time_usage},
                   ignore_index=True)
   
  else:
    df = df.append ({'n': train, 'mem': mem_usage, 'unit': mem_unit, 'time': time_usage},
                    ignore_index=True)
            
os.chdir("/dora/nobackup/yuec/scclassify/benchmark/ACTINN/vary_celltype")
            
df.to_csv("cpu_mem.csv")
            
            
            
            
         
