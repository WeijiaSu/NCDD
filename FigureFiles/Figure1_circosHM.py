import pandas as pd
import os
import os.path
import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
from cigar import Cigar
from collections import Counter
import pysam
import math

import warnings
warnings.filterwarnings('ignore')


pd.set_option("display.max_column",40)

def getCount(file,TE,Type):
	f=pd.read_table(file,header=0,sep="\t")
	f=f.loc[f["Circle"]!="NC"]
	f=f.sort_values(["Refname","Readname","RefStart"])
	f=f.groupby(["Refname","Readname"]).filter(lambda x: len(x)>1)
	
	t=f.loc[f["Refname"]==TE]
	tp=t.loc[t["Type"]==Type]
	count=tp.drop_duplicates(["Readname"],keep="first")
	print(count.shape[0])
	return tp	

def LTR1_frag(file,TE,side):
	f=pd.read_table(file)
	f=f.sort_values(["Refname","Readname","RefStart"])
	f=f.groupby(["Refname","Readname"]).filter(lambda x: len(x)>1)
	t=f.loc[f["Refname"]==TE]
	tp=t.loc[t["Type"]=="1End_Frg"]
	read_list=[]
	refLen=int(set(f["Reflen"])[0])
	for r in set(tp["Readname"]):
		sub=tp.loc[tp["Readname"]==r]
		if side==1:
			start=list(sub["RefStart"])
			if min(start)<=200:
				read_list.append(r)
		else:
			end=list(sub["RefEnd"])
			if max(end)>=refLen-200:
				read_list.append(r)
	side_df=tp.loc[tp["Readname"].isin(read_list)]
	read_number=side_df.drop_duplicates(["Readname"],keep="first")
	print(read_number.shape[0])
	return side_df


def getBed(f,refLen,TEname,sizeFactor,sampleName,Type):
	tp=f
	if tp.shape[0]==0:
		l=[[m,m+1] for m in range(0,reLen)]
		new_f=pd.DataFrame(l)
		new_f["v"]=0
		new_f["Refname"]=TEname
		new_f[["Refname",0,1,"v"]].to_csv(sampleName+"_"+TEname+"_"+Type+".bedgraph.expend.tsv",header=None,index=None,sep="\t")
	else:
		tp_bed=tp.drop_duplicates(["Refname","RefStart","RefEnd","Readname"],keep="first")
		tp_bed=tp[["Refname","RefStart","RefEnd"]]
		tp_bed=tp_bed.sort_values(["Refname","RefStart","RefEnd"])
		tp_bed.to_csv(sampleName+"_"+TEname+"_"+Type+".bed.tsv",header=None,index=None,sep="\t")
		g=tp[["Refname","Reflen"]].drop_duplicates(["Refname","Reflen"],keep="first")
		g.to_csv(TEname+"_TEsize",header=None,index=None,sep="\t")
		bedtools="bedtools genomecov -bga -i %s -g %s -scale %s> %s"% (sampleName+"_"+TEname+"_"+Type+".bed.tsv",TEname+"_TEsize",sizeFactor,sampleName+"_"+TEname+"_"+Type+".bedgraph.tsv")
		os.system(bedtools)
		bg=pd.read_table(sampleName+"_"+TEname+"_"+Type+".bedgraph.tsv",header=None)
		bg["c"]=bg[1].apply(str)+"_"+bg[2].apply(str)
		d=dict(zip(bg["c"],bg[3]))
		l=[]
		for k in d:
			v=d[k]
			for j in range(int(k.split("_")[0]),int(k.split("_")[1])):
				l.append([j,v])
		new_f=pd.DataFrame(l)
		new_f["Refname"]=TEname
		new_f["RefStart"]=new_f[0]
		new_f["RefEnd"]=new_f[0]+1
		new_f["Value"]=new_f[1]
		new_f=new_f[["Refname","RefStart","RefEnd","Value"]]
		new_f.to_csv(sampleName+"_"+TEname+"_"+Type+".bedgraph.expend.tsv",header=None,index=None,sep="\t")
		print(max(list(new_f["Value"])))	
		rm="rm %s %s %s "%(sampleName+"_"+TEname+"_"+Type+".bed.tsv",TEname+"_TEsize",sampleName+"_"+TEname+"_"+Type+".bedgraph.tsv")
		os.system(rm)

d2={}
d2["NegativeControlpreprocess"]=0.97
d2["PositiveControlpreprocess"]=1.06
d2["shcorollapreprocess"]=0.82
d2["shmei218preprocess"]=1.34
d2["shXpdpreprocess"]=1.06

dir="/data/zhanglab/Alexis_Bailey/"

for i in ["NegativeControlpreprocess","PositiveControlpreprocess","shcorollapreprocess","shmei218preprocess","shXpdpreprocess"]:
	HMS=getCount(dir+i+".fastq.minimap2.TE_full.sorted.bam_TE_full.fa.Type.txt_data.txt","HMS-Beagle","1End_FL")
	getBed(HMS,7062,"HMS-beagle",d2[i],i,"1End_FL")


#for i in range(1,6):
#	HMS=LTR1_frag("/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig4/rep1-210806_vas_Alt-EJ_circle-seq/barcode0%s.fastq.chop.fastq_HMS-Beagle_circleAnalyze.txt"%(i),"HMS-Beagle",1)
#	getBed(HMS,"HMS-beagle",d1[i],"barcode0%s_1_1LTR2"%(i))

#for i in range(1,6):
#	HMS=LTR1_frag("/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig4/rep2-211102-fly-vas_aEJ_KD-ovary-eccDNA/barcode0%s.fastq.chop.fastq_HMS-Beagle_circleAnalyze.txt"%(i),"HMS-Beagle",1)
#	getBed(HMS,"HMS-beagle",d2[i],"barcode0%s_2_1LTR2"%(i))


#for i in  ["FC_1LTR","FC_2LTR","PC_nonLTR"]:
#for i in ["1LTR2"]:
#	for j in range(1,6):
#		f1=pd.read_table("barcode0%s_1_%s_HMS-beagle_PC_1LTR.090922.20X.bedgraph.expend.tsv"%(j,i),header=None)
#		f2=pd.read_table("barcode0%s_2_%s_HMS-beagle_PC_1LTR.090922.20X.bedgraph.expend.tsv"%(j,i),header=None)
##		f1=pd.read_table("barcode0%s_1__HMS-beagle_%s.090922.20X.bedgraph.expend.tsv"%(j,i),header=None)
##		f2=pd.read_table("barcode0%s_2__HMS-beagle_%s.090922.20X.bedgraph.expend.tsv"%(j,i),header=None)
#		#f=f1
#		#f[3]=(f1[3]/2+f2[3]/2)
#		#f[3]=f[3].apply(lambda x: round(x,2))
#		f=f1.merge(f2,on=[0,1,2],how="inner")
#		f[5]=(f["3_x"]+f["3_y"])/2
#		f[5]=f[5].apply(lambda x: round(x,2))
#		f=f[[0,1,2,5]]
#		print(f1[0:10])
#		print(f2[0:10])
#		print(f[0:10])
#		f.to_csv("barcode0%s_HMS-beagle_%s.092122.bedgraph.expend.tsv"%(j,i),header=None,index=None,sep="\t")
#
