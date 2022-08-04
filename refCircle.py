import pandas as pd
import os
import argparse
import re
from Bio import SeqIO
import pysam
import os.path
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
import time
from cigar import Cigar
from collections import Counter
import numpy as np
import warnings


warnings.filterwarnings('ignore')
pd.set_option("display.max_column",40)

parser=argparse.ArgumentParser()
parser.add_argument("-bam","--inFile")
parser.add_argument("-js","--JunSize",default=100)
parser.add_argument("-Prefix","--Prefix",default=None)

args=parser.parse_args()


pd.set_option("display.max_columns",40)

bamFile=args.bamFile
pre=args.Prefix

from bamConverter import bamConverter
bamConverter().ConverAlignment(bamFile)

records=list(SeqIO.parse(ref,"fasta"))
d=dict(zip([rec.id for rec in records],[len(str(rec.seq)) for rec in records]))

def map_ratio(sub_f):
  r_len=int(list(sub_f["ReadLen"])[0])
  l=zip(sub_f["ReadStart"],sub_f["ReadEnd"])
  a=np.array([0]*r_len)
  for i in l:
    a[i[0]:i[1]+1]=1
  return list(a).count(1)

def getChimeric_reads(infile):
  f=pd.read_table(infile,header=0,sep="\t")
  f=f.sort_values(["Readname","ReadStart","ReadEnd"])
  s1=f.drop_duplicates(["Readname"],keep="first").shape[0]
  f["m"]=0
  f["infor"]=f["Readname"]+"_"+f["Refname"]
  for r in set(f["infor"]):
    sub_f=f.loc[f["infor"]==r]
    p=map_ratio(sub_f)
    f.loc[f["infor"]==r,"m"]=p
  f["p"]=f["m"]/f["ReadLen"]
  re=f.loc[f["p"]<0.9]
  f=f.loc[f["p"]>=0.9]
  f.to_csv(pre+"_rc90.tsv",index=None,sep="\t")
  s2=f.drop_duplicates(["Readname"],keep="first").shape[0]
  print("Total aligned reads %s; Fully aligned reads %s (%s)"%(s1,s2,round(s2/s1*100,2)))



def FilterReads(File):
    f=pd.read_table(File,header=0,sep="\t")
    f["Refname"]=f["Refname"].apply(str)
    f["Reflen"]=f["Refname"].apply(lambda x: d[x])
    LinearAlignment=f.loc[(f["ReadStart"]<=100) & (f["ReadEnd"]>=f["ReadLen"]-100)]
    LinearAlignment.to_csv(pre+"_LiAg.tsv",index=None,sep="\t")
    candidateReads=f.loc[f["Readname"].apply(lambda x: x not in list(LinearAlignment["Readname"]))]
    candidateReads.to_csv(pre+".candi.tsv",index=None,sep="\t")


def Junction(list1,list2):
    n1,n2,n3,n4=list1[0],list1[1],list2[0],list2[1]
    if n3>n2 and n3-n2>100:
        return False
    elif abs(n3-n1)<100 or abs(n2-n4)<100:
        return False
    else:
        return True


def isCircle(list1,list2):

    n1,n2,n3,n4=list1[0],list1[1],list2[0],list2[1]
    t1,t2,t3,t4=list1[2],list1[3],list2[2],list2[3]
    n_max=max(n1,n2,n3,n4)
    n_min=min(n1,n2,n3,n4)
    d1,d2=list1[4],list2[4]
    t_min=min(t1,t2,t3,t4)
    t_max=max(t1,t2,t3,t4)
    j=Junction(list1,list2)
    readlength=list1[-1]
    if d1!=d2 or j==False:
        return False
    else:
        if d1=="+" and t_max==t2:
            return True
        if d1=="-" and t_max==t4:
            return True
        else:
            return False

def CircleType(list1,list2):

    n1,n2,n3,n4=list1[0],list1[1],list2[0],list2[1]
    t1,t2,t3,t4=list1[2],list1[3],list2[2],list2[3]
    t_max=max(t1,t2,t3,t4)
    t_min=min(t1,t2,t3,t4)
    n_max=max(n1,n2,n3,n4)
    n_min=min(n1,n2,n3,n4)
    reflength=list1[-2]
    readlength=list1[-1]
    if isCircle(list1,list2) ==False:
        return "NC"
    else:
		return str(t_min)+"-"+str(t_max)


def getCircle(f):
    d={}
    infors=list(set(list(f["infor"])))
	for infor in infors:
		read=infor.spli("_")[0]
        d[infor]="NC"
        sub=f.loc[f["infor"]==infor]
        sub=sub.sort_values(["Refname","ReadStart","ReadEnd"])
        l=list(zip(sub["ReadStart"],sub["ReadEnd"],sub["RefStart"],sub["RefEnd"],sub["Strand"],sub["Reflen"],sub["ReadLen"]))
        i=0
        while i<len(l)-1:
            list1=l[i]
			list2=l[i+1]
            cirType=CircleType(list1,list2)
            if cirType!="NC":
               d[infor]=cirTyped
               break
            else:
               i+=1
    f["Circle"]=f["Readname"].apply(lambda x: d[x])
    return f


def GetReads(canfile):
    f_c=pd.read_table(canfile,header=0,sep="\t")
    f_c["infor"]=f_c["Readname"]+"_"+f["Refname"]
	f_c=f_c.sort_values(["infor","ReadStart","ReadEnd"])
	f_circle=getCircle(f_c)
    f_circle.to_csv(canfile+"_circleAnalyze.txt",index=None,sep="\t")
    fc=f_circle.loc[f_circle["Circle"]!="NC"]
    fc.to_csv(pre+"_circles.txt",index=None,sep="\t")
    nc=f_circle.loc[f_circle["Circle"]=="NC"]
    circle=f_circle.drop_duplicates(["Readname"],keep="first").groupby(["Circle"],as_index=False).count()[["Circle","Readname"]].sort_values(["Circle"])



getChimeric_reads(bamFile+"_AligTable.tsv")
FilterReads(pre+"_rc90.tsv")
Getreads(pre+".candi.tsv")
