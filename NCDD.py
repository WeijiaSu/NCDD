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
import pandas as pd

warnings.filterwarnings('ignore')
pd.set_option("display.max_column",40)

parser=argparse.ArgumentParser()
parser.add_argument("-bam","--bamFile")
parser.add_argument("-buffer","--Buffer",default=100)
parser.add_argument("-Prefix","--Prefix",default=None)
parser.add_argument("-ref","--reference")
args=parser.parse_args()


pd.set_option("display.max_columns",40)

bamFile=args.bamFile
pName=args.Prefix
ref=args.reference
fl=args.Buffer

from bamToPaf import bamConverter

#bamConverter().ConverAlignment(bamFile,pre)

#
def filterFullreads(mapped_paf):
	f=pd.read_table(mapped_paf)
	s1=f.drop_duplicates(["QName"],keep="first").shape[0]

	bedA=f[["QName","QLen"]].drop_duplicates(["QName"],keep="first")
	bedA["s"]=0
	bedA[["QName","s","QLen"]].to_csv(pName+".bedA.bed",header=None,index=None,sep="\t")
	f[["QName","QStart","QEnd"]].to_csv(pName+".bedB.bed",header=None,index=None,sep="\t")
	bedtools="bedtools coverage -a %s -b %s > %s"%(pName+".bedA.bed",pName+".bedB.bed",pName+".bed")
	os.system(bedtools)
	
	f_bed=pd.read_table(pName+".bed",header=None)
	f_bed["cutoff"]=f_bed[2]-f_bed[2]*f_bed[6]
	f_bed=f_bed.loc[(f_bed["cutoff"]<=fl) | (f_bed[6]>=0.9)]
	
	f=f.loc[f["QName"].isin(f_bed[0])]
	s2=f.drop_duplicates(["QName"],keep="first").shape[0]
	f.to_csv(pName+".filter.full.paf",index=None,sep="\t")
	os.remove(pName+".bedA.bed")
	os.remove(pName+".bedB.bed")
	os.remove(pName+".bed")
	print("Total aligned reads %s; Fully aligned reads %s (%s)"%(s1,s2,round(s2/s1*100,2)))	
#filterFullreads(pName+".paf")


def LinearReads(File):
    f=pd.read_table(File,header=0,sep="\t")
    f["RName"]=f["RName"].apply(str)
    LinearAlignment=f.loc[(f["QStart"]<=fl) & (f["QEnd"]>=f["QLen"]-fl)]
    #LinearAlignment=LinearAlignment.groupby("QName").filter(lambda x: len(x)==1)
    LinearAlignment.to_csv(pName+"_LiAg.tsv",index=None,sep="\t")
    candidateReads=f.loc[f["QName"].apply(lambda x: x not in list(LinearAlignment["QName"]))]
    candidateReads.to_csv(pName+".candi.tsv",index=None,sep="\t")
    s1=f.drop_duplicates(["QName"],keep="first").shape[0]
    s2=LinearAlignment.drop_duplicates(["QName"],keep="first").shape[0]
    s3=candidateReads.drop_duplicates(["QName"],keep="first").shape[0]
    
    print("Fully aligned reads %s: Linear Alignment %s (%s); Chimeric Alignment %s (%s)"%(s1,s2,round(s2/s1,2),s3,round(s3/s1,2)))
#LinearReads(pName+".filter.full.paf")

def ReadConfigure(list1,list2):
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
    j=ReadConfigure(list1,list2)
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

def CompleteCopy(sub,coordinate):
    sub["s"]=sub["RStart"].apply(lambda x: x>=int(coordinate[0])-100 and x<=int(coordinate[0])+100)
    sub["e"]=sub["REnd"].apply(lambda x: x>=int(coordinate[1])-100 and x<=int(coordinate[1])+100)
    sub=sub.loc[(sub["s"]==True) & (sub["e"]==True)]
    if sub.shape[0]>0:
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
    Overlap=n3-n2
    reflength=list1[-2]
    readlength=list1[-1]
    if isCircle(list1,list2) ==False:
        return "NC"
    else:
        return str(t_min)+"_"+str(t_max)+"_"+str(Overlap)


def getCircle(f):
    d={}
    infors=list(set(list(f["infor"])))
    for infor in infors:
        d[infor]="NC"
        sub=f.loc[f["infor"]==infor]
        sub=sub.sort_values(["QName","QStart","QEnd"])
        l=list(zip(sub["QStart"],sub["QEnd"],sub["RStart"],sub["REnd"],sub["Strand"],sub["RLen"],sub["QLen"]))
        i=0
        while i<len(l)-1:
            list1=l[i]
            j=i+1
            list2=l[i+1]
            cirType=CircleType(list1,list2)
            if cirType!="NC":
                if sub.shape[0]==2:
                    d[infor]=cirType
                    break
                else:
                    coordinate=cirType.split("_")
                    comp=CompleteCopy(sub,coordinate)
                    if comp==True:
                        d[infor]=cirType
                        break
                    else:
                        i+=1
            else:
                i+=1
    f["Circle"]=f["infor"].apply(lambda x: d[x])
    return f


def GetReads(canfile):
    f_c=pd.read_table(canfile,header=0,sep="\t")
    f_c["infor"]=f_c["QName"]+"_"+f_c["RName"]
    f_c=f_c.sort_values(["infor","QStart","QEnd","RStart","REnd"])
    f_circle=getCircle(f_c)
    f_circle.to_csv(canfile+"_circleAnalyze.txt",index=None,sep="\t")
    fc=f_circle.loc[f_circle["Circle"]!="NC"]
    fc.to_csv(pre+"_circles.txt",index=None,sep="\t")
    nc=f_circle.loc[f_circle["Circle"]=="NC"]
    circle=f_circle.drop_duplicates(["Readname"],keep="first").groupby(["Circle"],as_index=False).count()[["Circle","Readname"]].sort_values(["Circle"])
#
files=os.listdir("./")
if pName+".filter.full.paf" not in files:
    print("Geting full read")
    filterFullreads(pName+".paf")
if pName+".candi.tsv" not in files:
    print("Filtering reads")
    LinearReads(pName+".filter.full.paf")
GetReads(pName+".filter.full.paf")
