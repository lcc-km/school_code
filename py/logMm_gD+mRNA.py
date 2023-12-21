import os,re
rna={}
dna={}
os.chdir("/mnt/nas1/lucc/work/ourdata/zebarfish_mRNA/raw_data/raw_clean/sam/bam2/bam/add_head/vcf/snp/VQSR_10/vcf/heterozygosis/work/logMm")
with open("RNA_pos.txt") as f1:
    for l1 in f1:
        l11 = l1.rstrip("\n").split("\t")
        id1=l11[1]
        chr1=l11[2]
        POS=l11[3]
        ALE1=l11[4]
        ALE2=l11[5]
        READ1=l11[6]
        READ2=l11[7]
        ALL=l11[8]
        LOCATE=l11[9]
        INDEX=l11[10]
        start=l11[11]
        end=l11[12]
        rna[id1,INDEX,LOCATE]=chr1,POS,ALE1,ALE2,READ1,READ2,ALL,start,end
with open("D.work_clean") as f2:
    for l2 in f2:
        l22 = l2.rstrip("\n").split("\t")
        chr2=l22[0]
        pos=l22[1]
        ale1=l22[2]
        read1=l22[3]
        index1=l22[4]
        ale2=l22[5]
        read2=l22[6]
        index2=l22[7]
        dna[chr2,pos]=ale1,read1,index1,ale2,read2,index2
for a,b in rna.items():
    for c,d in dna.items():
        if b[0]==c[0] and b[7] <= c[1] <= b[8]:
            txt = open("RNA_pos.D", "a")
            print(a,b,c,d,"D",file=txt)
txt.close()