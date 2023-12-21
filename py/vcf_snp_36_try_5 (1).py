import os,re
os.chdir("/mnt/nas1/lucc/work/ourdata/zebarfish_mRNA/raw_data/clean/sam/bam2/bam/add_head/vcf")
gene={}
with open("Danio_rerio.GRCz11.94_gene.gtf") as f2:
    for l2 in f2:
        l22 = l2.rstrip("\n").split("\t")
        chr=l22[0]
        start=l22[1]
        end=l22[2]
        gene_id=l22[3]
        gene_name=l22[4]
        gene[chr,start,end]=gene_id,gene_name

file= ['RMA-1_FRRB190272334.vcf','RMB-6_FRRB190272345.vcf','RMC-4_FRRB190272349.vcf','RME-1-2_FRRB190276052-1a.vcf','RMB-2_FRRB190272341.vcf','RMC-3_FRRB190272348.vcf','RMD-2_FRRB190272353.vcf','RME-2-2_FRRB190276053-1a.vcf']      
for index in range(len(file)):
   vcf={} 
   with open (file[index]) as f:
     p = re.compile("#")
     lines = [line for line in f.readlines() if p.search(line) is None]
     for l1 in lines:
       l11 = l1.rstrip("\n").split("\t")
       chr = l11[0]
       site = l11[1]
       ref = l11[3]
       alt = l11[4]
       addp = l11[9]
       v = ref,alt,addp
       vcf[chr, site] = v
       #print('读取完成 line'+" "+file[index]+" "+ str(lines.index(l1)))
     name_file = file[index] + "_result.txt"
     txt=open(name_file,"w")
     for a,b in vcf.items():
        for c,d in gene.items():
               if a[0]==c[0]:
                  if int(c[1]) < int(a[1]) <int(c[2]):
                    print(str(d[0])+"\t"+str(d[1])+"\t"+str(a[::])+"\t"+str(b[::]),file=txt)
     txt.close()
   index=index+1
