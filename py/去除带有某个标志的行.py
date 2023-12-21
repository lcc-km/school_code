import os,re
os.chdir("/mnt/nas1/lucc/work/ourdata/zebarfish_gDNA/rwa_data/clean/sam/bam/full_data/SNPfisher_get_bam/try")
with open ('gD1-6.vcf.work') as f:
    p = re.compile("#")
    lines = [line for line in f.readlines() if p.search(line) is None]
f1=open('gD1-6.SNP.add','a')
f1.writelines(lines)
f1.close()










