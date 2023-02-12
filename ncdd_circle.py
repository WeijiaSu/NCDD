import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import pandas as pd
pd.set_option("display.max_columns",40)

# Get the input file name from command line arguments
res1=sys.argv[1]

def getFile(r):
	# Read the input file and store it in a pandas dataframe
	f=pd.read_table(r)
	f=f.drop_duplicates(["Readname"],keep="first")
	# Create columns "CirS" and "CirE" which contain the start and end values of the "Circle" column after splitting it
	f["CirS"]=f["Circle"].apply(lambda x: int(x.split("_")[0]))
	f["CirE"]=f["Circle"].apply(lambda x: int(x.split("_")[1]))
	circle=f[["Refname","CirS","CirE"]]
	circle=circle.sort_values(["Refname","CirS","CirE"])
	circle=circle.to_csv(res1+"circleCoor.tsv",header=None,index=None,sep="\t")

	bedcluster="bedtools cluster -i %s >%s"%(res1+"circleCoor.tsv",res1+"circleCoorCluster.tsv")
	os.system(bedcluster)
	
	cluster=pd.read_table(res1+"circleCoorCluster.tsv",header=None)
	cluster=cluster.groupby([3],as_index=False).count().sort_values([0],ascending=[False])
	print(cluster[0:10])
	print("Number of circles: %s"%(cluster.shape[0]))
	print("Most frequent circle: %s"%(cluster[0].max()))
getFile(res1)

