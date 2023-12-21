import os,re
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
os.chdir("C:/Users/49248/Desktop/meme")
fa=open(r'Danio_rerio.GRCz11.dna.chromosome.1.fa','r')
#fasta1 = SeqIO.to_dict(SeqIO.parse(fa,"fasta"))
#print(fasta1.id)
#print(repr(line.seq))
#print(len(line))
fasta={}
for line in SeqIO.parse(fa,"fasta"):
   print (line.id)
   print(repr(line.seq))
   print(len(line))
   seq=Seq(line.seq._data,IUPAC.unambiguous_dna)
   long = len(seq)
   ID=line.id
   fasta[ID]=seq
gtf183={}
with open('183.gtf') as result183:
      for line183 in result183:
         end183 = []
         start183 = []
         chr183 = []
         ID183 = []
         line183=line183.rstrip("\n").split("\t")
         chr183=line183[1]
         start183=line183[3]
         start183=int(start183)
         end183=line183[4]
         end183=int(end183)
         ID183=line183[6]
         gtf183[ID183]=([ID183],[chr183],[start183],[end183])
for k,v in gtf183:

       #break

















