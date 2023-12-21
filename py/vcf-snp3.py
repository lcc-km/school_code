import os,re
os.chdir("C:/Users/lcc/Desktop/snp")
vcf={}
with open ("SRR1524241_01.raw_hisat2") as f:
    for l1 in f:
        l11 = l1.rstrip("\n").split("\t")
        chr=l11[0]
        site=l11[1]
        ref=l11[3]
        alt=l11[4]
        addp=l11[9]
        addp=addp.split(':')
        AD=addp[1]
        AD=AD.split(',')
        REF=AD[0]
        ALT=AD[1]
        v=ref,alt,REF,ALT
        vcf[chr,site]=v
gene={}
with open("C:/Users/lcc/Desktop/myself/Danio_rerio.GRCz11.94.gene.gtf") as f2:
    for l2 in f2:
        l22 = l2.rstrip("\n").split("\t")
        chrg=l22[0]
        start=l22[1]
        end=l22[2]
        id=l22[4]
        name=l22[5]
        gene[chrg,start]=start,end,id,name
txt=open("SRR1524241_01.raw_hisat2_result.txt","w")
for a,b in vcf.items():
    for c,d in gene.items():
        if a[0]==c[0]:
            if int(d[0]) < int(a[1]) <int(d[1]):
                print(str(d[2:3])+"\t"+str(a[::])+"\t"+str(b[::]),file=txt)
txt.close