import os,re
rna={}
dna={}
os.chdir("/mnt/nas1/lucc/work/ourdata/zebarfish_mRNA/raw_data/raw_clean/sam/bam2/bam/add_head/vcf/snp/VQSR_10/vcf/heterozygosis/work/logMm/new/intersect_191205")
with open("183_gene_RNA_link_DNA_0_1") as f1:
    for l1 in f1:
        l11 = l1.rstrip("\n").split("\t")
        id1=l11[0]
        chr_ENSDARG=l11[1]
        s_10w=l11[2]
        s=l11[3]
        e=l11[4]
        e_10w=l11[5]
        name=l11[7]
        chr_RNA=l11[8]
        POS=l11[9]
        ALE1=l11[10]
        ALE2=l11[11]
        READ1=l11[13]
        READ2=l11[14]
        #all_read=l11[15]
        class_=l11[16]
        rna[id1,name,chr_RNA,POS,class_]=chr_ENSDARG,s_10w,s,e,e_10w,ALE1,ALE2,READ1,READ2
print("finish RNA")
with open("A.work_clean") as f2:
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
print("finish DNA")
for a,b in rna.items():
    for c,d in dna.items():
        if b[0]==c[0] and int(b[1]) <= int(c[1])  <= int(b[4]) :
            txt = open("RNA_pos.A", "a")
            print(a,b,c,d,file=txt)
txt.close()