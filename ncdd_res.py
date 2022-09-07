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
	cirCle_D=int(x.split("_")[-1])
	if cirCle_S <300 and cirCle_E>8557 and cirCle_D<-100:
		return "1LTR_FL"
	elif cirCle_S <300 and cirCle_E>8557 and cirCle_D>-100:
		return "2LTR_FL"
	elif cirCle_S <300 or cirCle_E>8557:
		return "1LTR_Frg"
	else:
		return "nonLTR_Frg"




def getFile(r,Out_name):
	f=pd.read_table(r)
	f=f.loc[f["Circle"]!="NC"]
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



