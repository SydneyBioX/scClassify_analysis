#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 21:25:25 2020

@author: yuec
"""



import numpy as np
import pandas as pd

from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV

from get_mem_python import display_top
import tracemalloc

import time
import os

def run_SVMreject( trainname, testname,   n):

    trainDataPath = '/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/' +  trainname + '.csv'
    trainLabelsPath = '/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/' + trainname +  '_label.csv'
    
    testDataPath = '/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/' +  testname + '.csv'
    testLabelsPath = '/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/' + testname +  '_label.csv'
   
    
    
    # read the data
    train = pd.read_csv(trainDataPath,  index_col = 0,  sep=',')
    test =  pd.read_csv(testDataPath, index_col = 0, sep = ',' )
    
    y_train = pd.read_csv(trainLabelsPath, header=0, index_col= 0 ,  sep=',')
    y_train = y_train['x'].ravel()
    y_test = pd.read_csv(testLabelsPath, header=0, index_col= 0 ,  sep=',')
    y_test = y_test['x'].ravel()


    truelab = []
    pred = []

    train = train.transpose()
    test = test.transpose()






    now = time.time() 
    tracemalloc.start() 
    
    Classifier = LinearSVC()
    clf = CalibratedClassifierCV(Classifier)
    

    clf.fit(train, y_train)

    snapshot = tracemalloc.take_snapshot()
    mem_train = display_top(snapshot)   
    
    later = time.time()
    time_train = int(later - now)
    
    
    
    
    
    now = time.time() 
    tracemalloc.start() 

    predicted = clf.predict(test)
    
    snapshot = tracemalloc.take_snapshot()
    mem_test = display_top(snapshot)   
    
    later = time.time()
    time_test = int(later - now)
    
    
    
    
    prob = np.max(clf.predict_proba(test), axis = 1)
 
    unlabeled = np.where(prob < 0.7)
    predicted[unlabeled] = 'Unassigned'


    truelab = y_test
    pred = predicted

    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)


    os.chdir("/dora/nobackup/yuec/scclassify/benchmark/SVMreject/vary_train")
 
    truelab.to_csv( n + "_SVMreject_true.csv", index = False)
    pred.to_csv( n +"_SVMreject_pred.csv", index = False)
  
    return mem_train, time_train, mem_test, time_test
  
  
  
  
  



os.chdir("/dora/nobackup/yuec/scclassify/benchmark/SVMreject/vary_train")
  
trainlist =  [ "train_100", "train_200", "train_500", "train_1000", "train_2000", 
                 "train_5000", "train_10000", "train_20000", "train_30000" ]

df = pd.DataFrame()

for i in range(0, len(trainlist)):
  train = trainlist[i]
  print(train)
 
  mem_train, time_train, mem_test, time_test = run_SVMreject(train, "test_2000" , train)
 
  if df.empty :
    df = df.append({'n': train, 'mem_train': mem_train, 'mem_test': mem_test, 
                    'unit': "MB", 'time_train': time_train, 'time_test': time_test },
                   ignore_index=True)
   
  else:
    df = df.append ({'n': train, 'mem_train': mem_train, 'mem_test': mem_test, 
                    'unit': "MB", 'time_train': time_train, 'time_test': time_test},
                    ignore_index=True)
            
os.chdir("/dora/nobackup/yuec/scclassify/benchmark/SVMreject/vary_train")
           
df.to_csv("cpu_mem.csv")
            
  