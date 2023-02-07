import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import pandas as pd
pd.set_option("display.max_columns",40)


res1=sys.argv[1]
name=sys.argv[2]


def getFile(r,Out_name):
	f=pd.read_table(r)
	f=f.loc[f["Circle"]!="NC"]
	f["D"]=f["Circle"].apply(lambda x: int(x.split("_")[-1]))
	f.loc[f["D"]>1,"JunFeature"]="Ins"
	f.loc[f["D"]<1."JunFeature"]="Del"
	f.loc[f["D"]==1,"JunFeature"]="Adj"
	f["CirS"]=f["Circle"].apply(lambda x: int(x.split("_")[0]))
	f["CirE"]=f["Circle"].apply(lambda x: int(x.split("_")[1]))
	f["CirSize"]=f["CirE"]-f["CirS"]+1
	f=f.drop(["m","infor","p","D","CirS","CirE"],axis=1)
	f.to_csv(Out_name,index=None,sep="\t",header=None)
	
	g=g.drop_duplicates(["Readname"],keep="first")
	print("Number of circle reads: %s"%(g.shape[0]))
	print("Largest circle size: %s"%(g["CirSize"].max()))
	print("Smallest circle size: %s"%(g["CirSize"].min()))
	print("Average circle size: %s"%(g["CirSize"].mean()))

getFile(res1,name)



