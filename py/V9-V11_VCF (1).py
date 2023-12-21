import os,re
os.chdir("D:\\backup\\our_data\\RNA\\paper_data")
with open("SNPfisher_FLI.vcf.bed") as f1:
    for l1 in f1:
        l11=l1.rstrip("\n").split("\t")
        chr_vcf = l11[0]