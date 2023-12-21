import os,re
os.chdir("/mnt/vol2_ws/test/wdd/projects/HRR/F300007230/snp_indel")
x={}

file= ['/mnt/vol2_ws/test/wdd/projects/HRR/F300007230/snp_indel/IVD-CA2021062809-phoenix-one-4-A1-210628L2/IVD-CA2021062809-phoenix-one-4-A1-210628L2.hg19_multianno_muts_somatic.xls',
'/mnt/vol2_ws/test/wdd/projects/HRR/F300007230/snp_indel/IVD-CA2021062809-phoenix-one-4-A2-210628L3/IVD-CA2021062809-phoenix-one-4-A2-210628L3.hg19_multianno_muts_somatic.xls',
'/mnt/vol2_ws/test/wdd/projects/HRR/F300007230/snp_indel/IVD-CA2021062809-phoenix-one-4-A3-210628L4/IVD-CA2021062809-phoenix-one-4-A3-210628L4.hg19_multianno_muts_somatic.xls',
'/mnt/vol2_ws/test/wdd/projects/HRR/F300007230/snp_indel/IVD-CA2021062809-phoenix-one-4-D3-210628L1/IVD-CA2021062809-phoenix-one-4-D3-210628L1.hg19_multianno_muts_somatic.xls',
'/mnt/vol2_ws/test/wdd/projects/HRR/F300007230/snp_indel/IVD-CA2021062816-phoenix-one-4-210628L17/IVD-CA2021062816-phoenix-one-4-210628L17.hg19_multianno_muts_somatic.xls',
'/mnt/vol2_ws/test/wdd/projects/HRR/F300007230/snp_indel/IVD-CA2021062816-phoenix-one-4-210628L26/IVD-CA2021062816-phoenix-one-4-210628L26.hg19_multianno_muts_somatic.xls']
for index in range(len(file)):
  sample={}
  l1 = index.rstrip("\n").split("\t")
  chr=l1[0]
  start=l1[1]
  end=l1[2]
  ref=l1[3]
  alt=l1[4]
  AF_Tumor=l1[8]
  Gene=l1[10]
  cDNA_change=l1[14]
  protein_change=l1[15]



  index = index + 1