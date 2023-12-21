import os,re
import pandas as pd
import numpy as np
import gc
from collections import OrderedDict
os.chdir("D:\\backup\\our_data\\DNA\\QC")
df1=pd.read_table("gD1-6.SNP_homozygosis",sep='\t',header=None)
df1['pos']=df1[0].map(str)+"_"+df1[1].map(str)
file= ['111']
for index in range(len(file)):
  df2=pd.read_table(file[index],sep='\t',header=None)
  df2['pos']=df2[0].map(str)+"_"+df2[1].map(str)
  x= pd.merge(df1,df2,on=['pos'],how = "inner")
  y=x[x['3_x']!=x['2_y']]
  name="1_in_6_"+file[index]
  y.to_csv(name, sep='\t', index=False,header=None)
  del y,x,df2,name
  gc.collect()
  index = index + 1


