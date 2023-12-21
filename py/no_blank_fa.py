#!/usr/bin/env python
import os
import sys
import re
import linecache
from collections import OrderedDict

"python get_fasta_seq XXX.fa XXX.blank.fa"
#python get_fasta_seq.py Danio_rerio.GRCz10.dna.toplevel.fa Danio_rerio.GRCz10.dna.toplevel.blank.fa
input_fa=sys.argv[1]
output_blank_fa=sys.argv[2]
fl=open(input_fa,"r")

fl_tem=open(output_blank_fa,"w")

for ln in fl:
	if ln.startswith(">"):
		fl_tem.write("\n"+ln)
	else:
		fl_tem.write(ln.rstrip("\n"))
fl.close()
fl_tem.close()