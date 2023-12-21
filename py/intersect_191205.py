import os,re
import pandas as pd
import numpy as np
import gc
os.chdir("D:\\backup\\our_data\\RNA\\VQSR\\10\\gene_area\\intersect_191205")
df1=pd.read_table("139_ID.txt.work",sep='\t')
df2=pd.read_table("ALL_RNA-seq_heterozygosisgene_result.txt3",sep='\t')
x = pd.merge(df1, df2, on=['ID'], how="inner")
name = "139ID_RNA_link_DNA"
x.to_csv(name, sep='\t', index=False)