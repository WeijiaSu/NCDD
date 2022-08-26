import pandas as pd
import sys


pd.set_option("display.max_column",40)

def GetMappingBed(Mapping):
	f=pd.read_table(Mapping)
	f=f.drop_duplicates(["Readname","ReadStart","ReadEnd"],keep="first")
	f=f.drop_duplicates(["Readname","Refname","RefStart","RefEnd"],keep="first")
	f[["Refname","RefStart","RefEnd"]].to_csv(Mapping+".bed",header=None,index=None,sep="\t")


def GetCircleBed(Circle):
	f=pd.read_table(Circle)
	f=f.loc[f["Circle"]!="NC"]
	f=f.drop_duplicates(["Refname","Readname"],keep="first")
	f["c_s"]=f["Circle"].apply(lambda x: int(x.split("-")[0]))
	f["c_e"]=f["Circle"].apply(lambda x: int(x.split("-")[1]))
	f[["Refname","c_s","c_e"]].to_csv(Circle+".bed",header=None,index=None,sep="\t")



Mapping=sys.argv[1]
Circle=sys.argv[2]

GetMappingBed(Mapping)
GetCircleBed(Circle)
