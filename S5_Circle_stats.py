import pandas as pd
import os
import argparse
import re
from Bio import SeqIO
import pysam

pd.set_option("display.max_columns",40)

parser=argparse.ArgumentParser()
parser.add_argument("-file","--inFile")
parser.add_argument("-Prefix","--Prefix",default=None)
args=parser.parse_args()

pre=args.Prefix
#if pre ==None:
#  pre=inputFile
inputFile=args.inFile

def CountCircles(File):
	f=pd.read_table(File,sep="\t")
	f=f.drop_duplicates(["Readname"],keep="first")
	l=[]
	for r in list(set(f["Refname"])):
		l_=[r]
		for i in ["1LTR_FL","2LTR_FL","1LTR_Frg","nonLTR_Frg"]:
			sub=f.loc[(f["Refname"]==r) & (f["Type"]==i)]
			sub=sub.drop_duplicates(["Readname"],keep="first")
			l_.append(sub.shape[0])
		l.append(l_)
	typeCount=pd.DataFrame(l,columns=["Refname","1LTR_FL","2LTR_FL","1LTR_Frg","nonLTR_Frg"])
	typeCount.to_csv(pre+"_circleType.tsv",index=None,sep="\t")

CountCircles(inputFile)
