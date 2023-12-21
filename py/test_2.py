#!/usr/bin/env python
import re,sys,os
from collections import OrderedDict

output_file="result.txt"
output=open(output_file,"w")

similar_dict=OrderedDict()
with open("similar.exam100","r") as similar:
	for line in similar:
		ke=line.split("\t")[0]+"_"+line.split("\t")[1]
		similar_dict[ke]=line

with open("pairs.exam") as pairs:
	for query_line in pairs:
		query_ele=query_line.split("\t")[0]+"_"+query_line.split("\t")[1]
		if query_ele in similar_dict:
			output.write(similar_dict[query_ele])
			#print(query_ele)
