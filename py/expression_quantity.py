import os,re
import pandas as pd
import numpy as np
import gc
os.chdir("D:\\backup\\our_data\\RNA\\VQSR\\10\\gene_area\\expression_quantity")
df1=pd.read_table("integrate_locate_196",sep='\t',header=None)
file= ['RME-1-2_FRRB190276052-1a','RME-2-2_FRRB190276053-1a','RME-3_FRRB190272360','RME-4_FRRB190272361','RME-5_FRRB190272362','RME-6_FRRB190272363','RMF-1_FRRB190272364','RMF-2_FRRB190272365','RMF-3_FRRB190272366','RMF-4_FRRB190272367','RMF-5_FRRB190272368','RMF-6_FRRB190272369']
for index in range(len(file)):
  df2=pd.read_table(file[index],sep='\t',header=None)
  pos1=['pos']
  pos2=['id','chr','locate','ale1','ale2','x','y','pos']
  df1.columns=pos1
  df2.columns=pos2
  x= pd.merge(df1,df2,on=['pos'],how = "inner")
  name=file[index]+"_locate_196"
  x.to_csv(name, sep='\t', index=False,header=None)
  del x,df2,name
  gc.collect()
  index = index + 1