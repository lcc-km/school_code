import os,re
os.chdir("/mnt/nas1/lucc/work/ourdata/zebarfish_gDNA/rwa_data/clean/sam/bam/full_data/SNPfisher_get_bam/try/result/done2/work/work2")
with open("F.work") as f1:
    for l1 in f1:
        l11 = l1.rstrip("\n").split("\t")
        x=len(l11)
        txt = open("F.work_clean","a")
        if x == 7:
            print(l11,file=txt)
txt.close()