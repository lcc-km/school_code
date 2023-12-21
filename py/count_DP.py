import os,re
import pandas as pd
import numpy as np
import gc
os.chdir("D:\\backup\\our_data\\DNA\\normaliz")
file= ['11']
for index in range(len(file)):
  data=pd.read_table("111",sep='\t',header=None)
  #data[4]=data[4].astype(int)
  data1=data[data[4]>=10]
  #data1=data1[data1[3]=="<*>"]
  print(data1.dtypes)
  name1= "/mnt/nas1/lucc/work/ourdata/zebarfish_gDNA/rwa_data/clean/sam/bam/full_data/SNPfisher_get_bam/try/add_QC/"+file[index] + "_QC"
  name2 ="/mnt/nas1/lucc/work/ourdata/zebarfish_gDNA/rwa_data/clean/sam/bam/full_data/SNPfisher_get_bam/try/add_QC/"+file[index] + "_add"
  data1.to_csv(name1, sep='\t', index=False,header=None)
  data2.to_csv(name2, sep='\t', index=False,header=None)
  del data1
  del data2
  del data
  del name1
  del name2
  gc.collect()
  index = index + 1