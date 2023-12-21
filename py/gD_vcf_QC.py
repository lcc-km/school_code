import os,re
os.chdir("/mnt/nas1/lucc/work/ourdata/zebarfish_gDNA/rwa_data/clean/sam/bam/full_data/SNPfisher_get_bam/try")
vcf={}
txt=open("gD1-6.vcf,QC","w")
with open("gD1-6.vcf.work") as f:
    p = re.compile("#")
    lines = [line for line in f.readlines() if p.search(line) is None]
    for l1 in lines:
        l11=l1.rstrip("\n").split("\t")
        chr=l11[0]
        pos=l11[1]
        ref=l11[2]
        ALT=l11[3]
        DP=l11[4]
        dp=DP.replace("DP=","")
        dp=int(dp)
        I16=l11[5]
        vcf[chr,pos]=ref,ALT,dp,I16
        print("read finish")
for a,b in vcf.items():
    if b[2] >= 10:
        print(a,b)
        print(a,b,file=txt)
txt.close()