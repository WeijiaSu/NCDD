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
parser.add_argument("-js","--JunSize",default=700)
parser.add_argument("-Prefix","--Prefix",default=None)

args=parser.parse_args()


pd.set_option("display.max_columns",40)

parser=argparse.ArgumentParser()
parser.add_argument("-bam","--bamFile")

args=parser.parse_args()

bamFile=args.bamFile


from bamConverter import bamConverter
bamConverter().ConverAlignment(bamFile)

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
  f.to_csv(infile+"_rc90.tsv",index=None,sep="\t")
  s2=f.drop_duplicates(["Readname"],keep="first").shape[0]
  print("Total aligned reads %s; Fully aligned reads %s (%s)"%(s1,s2,round(s2/s1*100,2)))


records=list(SeqIO.parse(ref,"fasta"))
d=dict(zip([rec.id for rec in records],[len(str(rec.seq)) for rec in records]))

def GetTEreads(File):
    f=pd.read_table(File,header=0,sep="\t")
    f["Refname"]=f["Refname"].apply(str)
    f["Reflen"]=f["Refname"].apply(lambda x: d[x])
    LinearAlignment=f.loc[(f["RefStart"]<=100) & (f["RefEnd"]>=f["ReadLen"]-100)]
    full_reads=LinearAlignment
    LinearAlignment.to_csv(pre+"_LiAg.tsv",index=None,sep="\t")
    Linear_ref_reads=full_reads.loc[(full_reads["RefStart"]<=100) & (full_reads["RefEnd"]>=full_reads["Reflen"]-100)]
    Linear_ref_reads.to_csv(pre+"_LiTEreads.tsv",index=None,sep="\t")
    short=full_reads.loc[~full_reads["Readname"].isin(list(Linear_ref_reads["Readname"]))]
    short.to_csv(pre+"_shortReads.tsv",index=None,sep="\t")
    #print("Total:", len(set(f["Readname"])))
    #print("TEreads:", len(set(Linear_ref_reads["Readname"])))
    #print("short reads:", len(set(short["Readname"])))
    candidateReads=f.loc[f["Readname"].apply(lambda x: x not in list(LinearAlignment["Readname"]))]
    #print("reads_left:", len(set(candidateReads["Readname"])))
    candidateReads.to_csv(pre+".candi.tsv",index=None,sep="\t")


def Junction(list1,list2):
    n1,n2,n3,n4=list1[0],list1[1],list2[0],list2[1]
    if n3>n2 and n3-n2>100:
        return False
    elif n3<n2 and n3-n2<-defJ:
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
        if d1=="+" and t_max==t2 and t_min==t3:
            return True
        if d1=="-" and t_max==t4 and t_min==t1:
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
        if t_max>=reflength-100 and t_min<=100 and n2-n3>=100:
            return "FC_1LTR"
        elif t_max>=reflength-100 and t_min<=100 and n2-n3<100:
            return "FC_2LTR"
        elif t_max<reflength-100 and t_min>100 and (t1!=t3 or t2!=t4):
            return "PC_nonLTR"
        elif (t_max>=reflength-100 or t_min<=100) and (t1!=t3 or t2!=t4):
            return "PC_1LTR"



def getCircle(f):
    d={}
    reads=list(set(list(f["Readname"])))
    for read in reads:
        d[read]="NC"
        sub=f.loc[f["Readname"]==read]
        sub=sub.sort_values(["Refname","ReadStart","ReadEnd"])
        l=list(zip(sub["ReadStart"],sub["ReadEnd"],sub["RefStart"],sub["RefEnd"],sub["Strand"],sub["Reflen"],sub["ReadLen"]))
        i=0
        while i<len(l)-1:
            j=i+1
            while j<len(l):
                list1=l[i]
                list2=l[j]
                cirType=CircleType(list1,list2)
                if cirType!="NC":
                    d[read]=cirType
                    break
                    i=len(l)
                else:
                    j=j+1
            i+=1
    f["Circle"]=f["Readname"].apply(lambda x: d[x])
    return f


def GetTEreads(canfile):
    f_c=pd.read_table(canfile,header=0,sep="\t")
    f_c=f_c.sort_values(["Readname","ReadStart","ReadEnd"])
    f_circle=getCircle(f_c)
    f_circle.to_csv(canfile+"_circleAnalyze.txt",index=None,sep="\t")
    fc=f_circle.loc[f_circle["Circle"]!="NC"]
    fc.to_csv(pre+"_circles.txt",index=None,sep="\t")
    nc=f_circle.loc[f_circle["Circle"]=="NC"]
    circle=f_circle.drop_duplicates(["Readname"],keep="first").groupby(["Circle"],as_index=False).count()[["Circle","Readname"]].sort_values(["Circle"])





getChimeric_reads(InputFile)
GetTEreads(inputFile)
#GetTEreads(inFile)
