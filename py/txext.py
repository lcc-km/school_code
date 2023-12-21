import os,re
import pandas as pd
import numpy as np
os.chdir("D:\\backup\\our_data\\DNA\\get_SNP\\text")
file= ['1','2']
for index in range(len(file)):
  data=pd.read_table(file[index],sep='\t',header=None)
  data1=data[data[4]>=10]
  data2=data1[data1[3]=="<*>"]
  name1= "D:\\backup\\our_data\\DNA\\get_SNP\\text\\"+file[index] + "_QC"
  name2 ="D:\\backup\\our_data\\DNA\\get_SNP\\text\\"+file[index] + "_add"
  data1.to_csv(name1, sep='\t', index=False,header=None)
  data2.to_csv(name2, sep='\t', index=False,header=None)
  del data1
  del data2
  del data
  del name1
  del name2
  index = index + 1