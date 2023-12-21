#!/usr/bin/python
import sys,os,re
os.chdir("C:/Users/49248/Desktop")
def process_file(reader):
    '''Open, read,and print a file'''
    names=[]
    index=0
    dict={}
    for line in reader:
        if line.startswith('>'):
           if index >=1:
               names.append(line)
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
        print"%s\t%d" %(key,length)
    input_file.close()