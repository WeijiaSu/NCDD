import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import pandas as pd
pd.set_option("display.max_columns",40)


res1=sys.argv[1]
name=sys.argv[2]

def circleType(x):
	cirCle_S=int(x.split("_")[0])
	cirCle_E=int(x.split("_")[1])
	cirCle_D=int(x.split("_")[2])
	cirCle_L=int(x.split("_")[-1])
	if cirCle_S <100 and cirCle_E>cirCle_L-100 and cirCle_D<=-100:
		return "1End_FL"
	elif cirCle_S <100 and cirCle_E>cirCle_L-100 and cirCle_D>-100:
		return "2Ends_FL"
	elif cirCle_S <100 or cirCle_E>cirCle_L-100:
		return "1End_Frg"
	else:
		return "noEnd_Frg"
		


def getFile(r,Out_name):
	f=pd.read_table(r)
	f=f.loc[f["Circle"]!="NC"]
	f=f.drop(["m","infor","p"],axis=1)
	f["circileInfor"]=f["Circle"]+"_"+f["Reflen"].apply(str)
	f=f.drop_duplicates(["Readname"],keep="first")
	
	f["Type"]=f["Circle"].apply(lambda x : circleType(x))
	
	d=dict(zip(list(f["Readname"]),list(f["Type"])))

	f_=pd.read_table(r)
	f_=f_.loc[f_["Circle"]!="NC"]
	f_["Type"]=f_["Readname"].apply(lambda x: d[x])
	f_.to_csv(Out_name+"_data.txt",index=None,sep="\t")
	g=f.groupby(["Type"],as_index=False).count()
	g=g[["Type","Readname"]]
	print(g)
	g.to_csv(Out_name,index=None,sep="\t",header=None)

getFile(res1,name)



