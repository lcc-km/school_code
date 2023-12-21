import os,re
import pandas as pd
import numpy as np
os.chdir("D:\\backup\\our_data\\DNA\\get_SNP")
data=pd.read_table("gD1-6.snp.clean.homozygosis",sep='\t',header=None)
data1=data[data[3]>=10]
#data2=data1[data1[3]=="<*>"]
data1.to_csv('gD1-6.snp.clean.homozygosis.QC', sep='\t', index=False,header=None)