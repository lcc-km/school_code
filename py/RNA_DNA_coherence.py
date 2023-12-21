import os,re
import pandas as pd
import numpy as np
import gc
os.chdir("D:\\backup\\our_data\\RNA\\VQSR\\10\\gene_area\\expression_quantity\\integrate_locate\\all\\logMm\\RNA_DNA_coherence")

df1=pd.read_table("RMA-1_FRRB190272334.pos_basic",sep='\t',header=None)
df1[2].replace({'A':0 ,'T':1,'C':7,'G':10},inplace= True)
df1[3].replace({'A':0 ,'T':1,'C':7,'G':10},inplace= True)
df1['pos']=df1[0].map(str)+"_"+df1[1].map(str)
df1['GT']=df1.sum(axis=1)
#df1.to_csv("df1", sep='\t', index=False,header=True)

df2=pd.read_table("A.work_clean",sep='\t',header=None)
df2[2].replace({' A':0 ,' T':1,' C':7,' G':10},inplace= True)
df2[5].replace({' A':0 ,' T':1,' C':7,' G':10},inplace= True)
df2['pos']=df2[0].map(str)+"_"+df2[1].map(str)

df2.to_csv("df2", sep='\t', index=False,header=True)

x = pd.merge(df1, df2, on=['pos'], how="inner")



x["2_y"].replace({'A':0 ,'T':1,'C':7,'G':10},inplace= True)
x.rename(columns={5:'5_y'},inplace=True)
x["5_y"].replace({'A':0 ,'T':1,'C':7,'G':10},inplace= True)

x.to_csv("x1", sep='\t', index=False,header=True)
