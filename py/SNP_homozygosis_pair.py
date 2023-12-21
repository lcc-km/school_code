import os,re
import pandas as pd
import numpy as np
import gc
os.chdir("/mnt/nas2/lcc/SNP2/done/logMm")
df1=pd.read_table("all_homo_add_locate",sep='\t',header=None)
df1['pos']=df1[0].map(str)+"_"+df1[1].map(str)
file=['gD6-7-1.vcf_clean_done_rmIndelab_QC.cutac','gD1-6.vcf_clean_done_rmIndelac_QC.cutac','gD2-6.vcf_clean_done_rmIndelac_QC.cutac','gD3-6.vcf_clean_done_rmIndelac_QC.cutac','gD4-4.vcf_clean_done_rmIndelac_QC.cutac','gD5-6.vcf_clean_done_rmIndelac_QC.cutac','gD6-7-1.vcf_clean_done_rmIndelac_QC.cutac','gD1-6.vcf_clean_done_rmIndelad_QC.cutac','gD2-6.vcf_clean_done_rmIndelad_QC.cutac','gD3-6.vcf_clean_done_rmIndelad_QC.cutac','gD4-4.vcf_clean_done_rmIndelad_QC.cutac','gD5-6.vcf_clean_done_rmIndelad_QC.cutac','gD6-7-1.vcf_clean_done_rmIndelad_QC.cutac','gD1-6.vcf_clean_done_rmIndelae_QC.cutac','gD2-6.vcf_clean_done_rmIndelae_QC.cutac','gD3-6.vcf_clean_done_rmIndelae_QC.cutac','gD4-4.vcf_clean_done_rmIndelae_QC.cutac','gD5-6.vcf_clean_done_rmIndelae_QC.cutac','gD6-7-1.vcf_clean_done_rmIndelae_QC.cutac']
for index in range(len(file)):
  df2=pd.read_table(file[index],sep='\t',header=None)
  df2['pos']=df2[0].map(str)+"_"+df2[1].map(str)
  x= pd.merge(df1,df2,on=['pos'],how = "inner")
  del df1,df2
  gc.collect()
  name = file[index]+".cut"
  x.to_csv(name,sep='\t', index=False,header=None)
  del x,name
  gc.collect()