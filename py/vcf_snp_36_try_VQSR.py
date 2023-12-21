import os,re
os.chdir("/mnt/nas1/lucc/work/ourdata/zebarfish_mRNA/raw_data/raw_clean/sam/bam2/bam/add_head/vcf/snp/VQSR_10/vcf/heterozygosis/work")
gene={}
with open("Danio_rerio.GRCz11.94_gene.gtf") as f2:
    for l2 in f2:
        l22 = l2.rstrip("\n").split("\t")
        chr=l22[0]
        start=l22[1]
        end=l22[2]
        strand=l22[3]
        gene_id=l22[4]
        gene_name=l22[5]
        gene[chr,start,end]=strand,gene_id,gene_name
file= ['RMB-2_FRRB190272341.heterozygosis','RMC-4_FRRB190272349.heterozygosis','RMD-6_FRRB190272357.heterozygosis','RMF-2_FRRB190272365.heterozygosis','RMA-1_FRRB190272334.heterozygosis','RMB-3_FRRB190272342.heterozygosis','RMC-5_FRRB190272350.heterozygosis','RME-1-2_FRRB190276052-1a.heterozygosis','RMF-3_FRRB190272366.heterozygosis','RMA-2_FRRB190272335.heterozygosis','RMB-4_FRRB190272343.heterozygosis','RMC-6_FRRB190272351.heterozygosis','RME-2-2_FRRB190276053-1a.heterozygosis','RMF-4_FRRB190272367.heterozygosis','RMA-3_FRRB190272336.heterozygosis','RMB-5_FRRB190272344.heterozygosis','RMD-1_FRRB190272352.heterozygosis','RME-3_FRRB190272360.heterozygosis','RMF-5_FRRB190272368.heterozygosis','RMA-4_FRRB190272337.heterozygosis','RMB-6_FRRB190272345.heterozygosis','RMD-2_FRRB190272353.heterozygosis','RME-4_FRRB190272361.heterozygosis','RMF-6_FRRB190272369.heterozygosis','RMA-5_FRRB190272338.heterozygosis','RMC-1_FRRB190272346.heterozygosis','RMD-3_FRRB190272354.heterozygosis','RME-5_FRRB190272362.heterozygosis','RMA-6_FRRB190272339.heterozygosis','RMC-2_FRRB190272347.heterozygosis','RMD-4_FRRB190272355.heterozygosis','RME-6_FRRB190272363.heterozygosis','RMB-1_FRRB190272340.heterozygosis','RMC-3_FRRB190272348.heterozygosis','RMD-5_FRRB190272356.heterozygosis','RMF-1_FRRB190272364.heterozygosis']      
for index in range(len(file)):
   vcf={} 
   with open (file[index]) as f:
     for l1 in f:
       l11 = l1.rstrip("\n").split("\t")
       chr = l11[0]
       site = l11[1]
       ref = l11[3]
       alt = l11[4]
       addp = l11[9]
       v = ref,alt,addp
       vcf[chr, site] = v
       #print('读取完成 line'+" "+file[index]+" "+ str(lines.index(l1)))
     name_file = file[index] + "gene_result.txt"
     txt=open(name_file,"w")
     for a,b in vcf.items():
        for c,d in gene.items():
               if a[0]==c[0]:
                  if int(c[1]) < int(a[1]) <int(c[2]):
                    print(str(d)+"\t"+str(a[::])+"\t"+str(b[::]),file=txt)
     txt.close()
   index=index+1
