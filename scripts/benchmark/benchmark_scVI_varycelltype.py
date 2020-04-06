#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 21:13:21 2020

@author: yuec
"""



from scvi.dataset import CsvDataset
import os
import pandas as pd
from scvi.models import SCANVI
from scvi.inference import SemiSupervisedTrainer
import time

from get_mem_python import display_top
import tracemalloc


def run_scVI(trainname, testname , n ):

    #trainDataPath = "/Users/yue/Dropbox (Sydney Uni)/scclassify/scRNAseq_Benchmark_datasets/Pancreatic_data/Segerstolpe/Filtered_Segerstolpe_HumanPancreas_data.csv"
    #train = pd.read_csv(trainDataPath,index_col=0,sep=',')
    #trainLabelsPath =  "/Users/yue/Dropbox (Sydney Uni)/scclassify/scRNAseq_Benchmark_datasets/Pancreatic_data/Segerstolpe/Labels.csv"
    #trainlabels = pd.read_csv(trainLabelsPath, header=0,index_col=None, sep=',')
    
    
    #testDataPath = "/Users/yue/Dropbox (Sydney Uni)/scclassify/scRNAseq_Benchmark_datasets/Pancreatic_data/Xin/Filtered_Xin_HumanPancreas_data.csv"
    #test = pd.read_csv(testDataPath,index_col=0,sep=',')
    #testLabelsPath =  "/Users/yue/Dropbox (Sydney Uni)/scclassify/scRNAseq_Benchmark_datasets/Pancreatic_data/Xin/Labels.csv"
    #testlabels = pd.read_csv(testLabelsPath, header=0,index_col=None, sep=',')
    
    

    train = pd.read_csv('/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/'+trainname+'.csv',index_col=0,sep=',')
    test = pd.read_csv('/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/'+testname+'.csv',index_col=0,sep=',')
    trainlabel = pd.read_csv('/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/'+trainname+'_label.csv', header=0,index_col=0, sep=',')
    testlabel = pd.read_csv('/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/'+testname+'_label.csv',header=0,index_col=0, sep=',')
    
    
    commongene =   train.index & test.index
    train = train.loc[commongene ]
    test = test.loc[commongene ]
    
    newdata = pd.concat([train, test], axis =1 )
    newlabel = pd.concat([trainlabel, testlabel], axis=0)
    
   
    #train = '/Users/yue/Dropbox (Sydney Uni)/scclassify/countmatrix/logcount/xin.csv'
   
 
    #save labels as csv file with header and index column
    #trainlabels.to_csv('trainLabels_scvi.csv')
    #train.to_csv('trainData_scvi.csv')    
    
    #testlabels.to_csv('testLabels_scvi.csv')
    #test.to_csv('testData_scvi.csv')
    
    
        
    os.chdir("/dora/nobackup/yuec/scclassify/benchmark/scVI/vary_celltype")
    
    newdata.to_csv('data_scvi.csv')
    newlabel.to_csv('labels_scvi.csv')
    data = CsvDataset('data_scvi.csv', save_path = "", sep = ",", labels_file = "labels_scvi.csv", gene_by_cell = True)
       
 
    n_epochs = 100
    
    truelab = []
    pred = []
    
    
    
    
    
    ## this semisupervised trainer automatically uses a part of the input data for training and a part for testing
    
    now = time.time() 
    tracemalloc.start() 
    
    scanvi = SCANVI(data.nb_genes,data.n_batches, data.n_labels)
    trainer_scanvi = SemiSupervisedTrainer(scanvi, data, frequency=5)
  
    trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(list(range(0,trainlabel.shape[0]))), shuffle = False)
    trainer_scanvi.labelled_set.to_monitor = ['ll','accuracy']
    trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=(list(range(trainlabel.shape[0],trainlabel.shape[0] + testlabel.shape[0]))), shuffle = False)
    trainer_scanvi.unlabelled_set.to_monitor = ['ll','accuracy']
      
    trainer_scanvi.train(n_epochs)
    
    snapshot = tracemalloc.take_snapshot()
    mem_train = display_top(snapshot)   
    
    later = time.time()
    time_train = int(later - now)
    
    
    
        ## labels of test set are in y_pred
        ## labels are returned in numbers, should be mapped back to the real labels
        ## indices are permutated
    
    
    now = time.time()
    tracemalloc.start() 
    
    y_true, y_pred = trainer_scanvi.unlabelled_set.compute_predictions()
     
    snapshot = tracemalloc.take_snapshot()
    mem_test = display_top(snapshot)   
    
    later = time.time()
    time_test = int(later - now)
    
    
   
    truelab.extend(y_true)
    pred.extend(y_pred)
    
    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)
    
    os.chdir("/dora/nobackup/yuec/scclassify/benchmark/scVI/vary_celltype")
    
 
    truelab.to_csv(n + "_scVI_True.csv", index = False)
    pred.to_csv(n + "_scVI_Pred.csv", index = False)
    
    
    return mem_train, time_train, mem_test, time_test
    
    






os.chdir("/dora/nobackup/yuec/scclassify/benchmark/scVI/vary_celltype")
  
trainlist =  [ "train_2class", "train_4class", "train_6class", "train_8class", "train_10class", "train_12class" ]



df = pd.DataFrame()

for i in range(0, len(trainlist)):
  train = trainlist[i]
  print(train)
 
  mem_train, time_train, mem_test, time_test = run_scVI(train , "test_5000", train)
 
  if df.empty :
    df = df.append({'n': train, 'mem_train': mem_train, 'mem_test': mem_test, 
                    'unit': "MB", 'time_train': time_train, 'time_test': time_test },
                   ignore_index=True)
   
  else:
    df = df.append ({'n': train, 'mem_train': mem_train, 'mem_test': mem_test, 
                    'unit': "MB", 'time_train': time_train, 'time_test': time_test},
                    ignore_index=True)
            
os.chdir("/dora/nobackup/yuec/scclassify/benchmark/scVI/vary_celltype")
           
df.to_csv("cpu_mem.csv")
            


       