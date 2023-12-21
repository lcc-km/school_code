import os,re
os.chdir("D:\\backup\\our_data\\DNA")
vcfqc={}
file= ['gD1-6_BDSW190046442.vcf','gD2-6_BDSW190046443.vcf','gD3-6_BDSW190046444.raw.vcf','gD6-7-1_BDSW190048305.raw.vcf']
for index in range(len(file)):
 with open ("file[index]") as f:
    p = re.compile("#")
    lines = [line for line in f.readlines() if p.search(line) is None]
 for l1 in lines:
    l11=l1.rstrip("\n").split("\t")
    #chrom=l11[0]
    #pos=l11[1]
    x=l11[7]
    xx=x.split(";")
    MQ=""
    name = file[index] + "SOR"
    txt = open(name, "a")
    for a in xx:
        if a.startswith("SOR="):
            MQ=print(a,file=txt)
            MQ = print(a)
        txt.close()

