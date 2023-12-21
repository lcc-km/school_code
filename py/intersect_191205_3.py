import os,re
import pandas as pd
import numpy as np
import gc
os.chdir("D:\\backup\\our_data\\RNA\\VQSR\\10\\gene_area\\intersect_191205")
df1=pd.read_table("139ID_RNA_link_DNA_0_1",sep='\t',header=None)
df1['pos']=df1[8].map(str)+"_"+df1[9].map(str)
df2=pd.read_table("integrate_pos_5.work",sep='\t')
x= pd.merge(df1,df2,on=['pos'],how = "inner")
name = "139ID_RNA_link_DNA_0_1.work"
x.to_csv(name, sep='\t', index=False)