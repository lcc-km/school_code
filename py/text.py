import os,re
import pandas as pd
import numpy as np
import gc
pd.set_option("display.max_columns",15)
os.chdir("D:\\backup\\our_data\\RNA\\VQSR\\10\\gene_area\\intersect_202016")
#file=["RMA-1_FRRB190272334.heterozygosisgene_result.txt_0_1","RMD-1_FRRB190272352.heterozygosisgene_result.txt_0_1","RMA-2_FRRB190272335.heterozygosisgene_result.txt_0_1","RMD-2_FRRB190272353.heterozygosisgene_result.txt_0_1","RMA-3_FRRB190272336.heterozygosisgene_result.txt_0_1","RMD-3_FRRB190272354.heterozygosisgene_result.txt_0_1","RMA-4_FRRB190272337.heterozygosisgene_result.txt_0_1","RMD-4_FRRB190272355.heterozygosisgene_result.txt_0_1","RMA-5_FRRB190272338.heterozygosisgene_result.txt_0_1","RMD-5_FRRB190272356.heterozygosisgene_result.txt_0_1","RMA-6_FRRB190272339.heterozygosisgene_result.txt_0_1","RMD-6_FRRB190272357.heterozygosisgene_result.txt_0_1","RMB-1_FRRB190272340.heterozygosisgene_result.txt_0_1","RME-1-2_FRRB190276052-1a.heterozygosisgene_result.txt_0_1","RMB-2_FRRB190272341.heterozygosisgene_result.txt_0_1","RME-2-2_FRRB190276053-1a.heterozygosisgene_result.txt_0_1","RMB-3_FRRB190272342.heterozygosisgene_result.txt_0_1","RME-3_FRRB190272360.heterozygosisgene_result.txt_0_1","RME-4_FRRB190272361.heterozygosisgene_result.txt_0_1","RMB-5_FRRB190272344.heterozygosisgene_result.txt_0_1","RME-5_FRRB190272362.heterozygosisgene_result.txt_0_1","RMB-6_FRRB190272345.heterozygosisgene_result.txt_0_1","RME-6_FRRB190272363.heterozygosisgene_result.txt_0_1","RMC-1_FRRB190272346.heterozygosisgene_result.txt_0_1","RMF-1_FRRB190272364.heterozygosisgene_result.txt_0_1","RMC-2_FRRB190272347.heterozygosisgene_result.txt_0_1","RMF-2_FRRB190272365.heterozygosisgene_result.txt_0_1","RMC-3_FRRB190272348.heterozygosisgene_result.txt_0_1","RMF-3_FRRB190272366.heterozygosisgene_result.txt_0_1","RMC-4_FRRB190272349.heterozygosisgene_result.txt_0_1","RMF-4_FRRB190272367.heterozygosisgene_result.txt_0_1","RMC-5_FRRB190272350.heterozygosisgene_result.txt_0_1","RMF-5_FRRB190272368.heterozygosisgene_result.txt_0_1","RMC-6_FRRB190272351.heterozygosisgene_result.txt_0_1","RMF-6_FRRB190272369.heterozygosisgene_result.txt_0_1"]
#for index in range(len(file)):
df1=pd.read_table("RMA-1_FRRB190272334.heterozygosisgene_result.txt_0_1",sep='\t',header=None)
df1[5].replace({'A':0 ,'T':1,'C':7,'G':10},inplace= True)
df1[6].replace({'A':0 ,'T':1,'C':7,'G':10},inplace= True)
name = "RMA-1_FRRB190272334.heterozygosisgene_result.txt_0_1" + "work"
df1.to_csv(name, sep='\t', index=False,header=None)
del df1,name
gc.collect()
#index = index + 1