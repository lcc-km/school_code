#!/usr/bin/python
import sys,os,re
os.chdir("D:\\backup\\our_data\\DNA\\get.bam.fa")
id={}
with open ("A") as f2:
    for l2 in f2:
        id[l2]=l2
def process_file(reader):
    '''Open, read,and print a file'''
    names=[]
    index=0
    dict={}
    for line in reader:
      for a, b in id.items():
        if line.startswith('@'):
           if index >=1:
               names.replace("@HD\tSO:coordinate\t\t\t\n",a)
           index =index+1
           name=line[:-1]
           seq = ''
        else:
           seq +=line[:-1]
           dict[name]=seq
    return dict
if __name__ == "__main__":
    input_file=open(sys.argv[1],"r")
    reader=input_file.readlines()
    items=process_file(reader)
    for key in items:
        length=int(len(items[key]))
        print(key,length)
    input_file.close()