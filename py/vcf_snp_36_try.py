import os,re
os.chdir("C:/Users/lcc/Desktop/our_data/RNA/clean/vcf")
file= ['RMA-5_FRRB190272338.raw.vcf','RMD-1_FRRB190272352.raw.vcf','RME-5_FRRB190272362.raw.vcf','RME-6_FRRB190272363.raw.vcf','RMD-6_FRRB190272357.raw.vcf','RMC-1_FRRB190272346.raw.vcf','RMD-4_FRRB190272355.raw.vcf','RMF-4_FRRB190272367.raw.vcf','RMD-3_FRRB190272354.raw.vcf','RMF-3_FRRB190272366.raw.vcf','RMB-1_FRRB190272340.raw.vcf','RMA-6_FRRB190272339.raw.vcf','RMB-4_FRRB190272343.raw.vcf','RMB-5_FRRB190272344.raw.vcf','RMD-5_FRRB190272356.raw.vcf','RME-3_FRRB190272360.raw.vcf','RMA-3_FRRB190272336.raw.vcf','RME-4_FRRB190272361.raw.vcf','RMF-1_FRRB190272364.raw.vcf','RMA-2_FRRB190272335.raw.vcf','RMF-2_FRRB190272365.raw.vcf','RMC-6_FRRB190272351.raw.vcf','RMC-5_FRRB190272350.raw.vcf','RMF-6_FRRB190272369.raw.vcf','RMC-2_FRRB190272347.raw.vcf','RMA-4_FRRB190272337.raw.vcf','RMB-3_FRRB190272342.raw.vcf','RMF-5_FRRB190272368.raw.vcf']
vcf={}
for index in range(len(file)):
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
     v = ref, alt, addp
     vcf[chr, site] = v
     name = file[index] + "_result.txt"
     for a,b in vcf.items():
         txt = open(name, "a")
         print(a,b,file=txt)
         txt.close()