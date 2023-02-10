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
	circle=circle.to_csv(res1+"circleCoor.tsv",index=None,sep="\t")

getFile(res1)

