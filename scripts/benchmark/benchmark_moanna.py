#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 19:05:49 2020

@author: yuec
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 18:08:49 2019

@author: yue
"""

import os
import pandas as pd
import numpy as np

import moana
from moana.core import ExpMatrix
from moana.core import CellAnnVector
from moana.classify import CellTypeClassifier
import time as tm
import time

from get_mem_python import display_top
import tracemalloc

def run_moana(trainname, testname, n):

    
    
    DataPath  = '/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/'+trainname+'.csv'  
    matrix = ExpMatrix.read_tsv(DataPath, sep = ',')    
    LabelsPath = '/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/'+trainname+'_label.csv'
    truelab = pd.read_csv(LabelsPath, header=0,index_col=0, sep=',')
    data = ExpMatrix(X = matrix.X, genes = matrix.genes, cells = matrix.cells)
    
    data.genes.name = 'Genes'
    data.cells.name = 'Cells'
    data.index.name = 'Genes'
    data.columns.name = 'Cells'
    
    l = CellAnnVector(cells=data.cells, data=truelab['x'].values)
    
    
    
    
    now = time.time()
    
    tracemalloc.start() 
    clf = CellTypeClassifier()
    clf.fit(matrix = data, cell_labels = l)
    snapshot = tracemalloc.take_snapshot()
    mem_train = display_top(snapshot)
    
    later = time.time()
    time_train = int(later - now)
  

    
    DataPath  = '/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/'+testname+'.csv'  
    matrix = ExpMatrix.read_tsv(DataPath, sep = ',') 
    data = ExpMatrix(X = matrix.X, genes = matrix.genes, cells = matrix.cells)
    data.genes.name = 'Genes'
    data.cells.name = 'Cells'
    data.index.name = 'Genes'
    data.columns.name = 'Cells'
        
    
    
    now = time.time()
    
    tracemalloc.start() 
    predictions = clf.predict(data) 
    snapshot = tracemalloc.take_snapshot()
    mem_test = display_top(snapshot)
    
    later = time.time()
    time_test = int(later - now)
  

    
    
    
    predictions = np.asarray(predictions)
    pred = pd.DataFrame(predictions)
        
        
    LabelsPath = '/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/'+testname+'_label.csv'
    truelab =  pd.read_csv(LabelsPath, header=0,index_col=None, sep=',')
    
    os.chdir("/dora/nobackup/yuec/scclassify/benchmark/moanna/vary_train")
       
        
    truelab.to_csv(n + "_moana_True.csv", index = False)
    pred.to_csv( n  + "_moana_Pred.csv", index = False)
    
      
    return mem_train, time_train, mem_test, time_test






os.chdir("/dora/nobackup/yuec/scclassify/benchmark/moanna/vary_train")
  
trainlist =  [ "train_100", "train_200", "train_500", "train_1000", "train_2000", 
                 "train_5000", "train_10000", "train_20000", "train_30000" ]

df = pd.DataFrame()

for i in range(0, len(trainlist)):
  train = trainlist[i]
  print(train)
 
  mem_train, time_train, mem_test, time_test = run_moana(train, "test_2000" , train)
 
  if df.empty :
    df = df.append({'n': train, 'mem_train': mem_train, 'mem_test': mem_test, 
                    'unit': "MB", 'time_train': time_train, 'time_test': time_test },
                   ignore_index=True)
   
  else:
    df = df.append ({'n': train, 'mem_train': mem_train, 'mem_test': mem_test, 
                    'unit': "MB", 'time_train': time_train, 'time_test': time_test},
                    ignore_index=True)
            
os.chdir("/dora/nobackup/yuec/scclassify/benchmark/moanna/vary_train")
           
df.to_csv("cpu_mem.csv")
            


