import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import pandas as pd
pd.set_option("display.max_columns",40)

res1=sys.argv[1]

def getFile(r):
	f=pd.read_table(r)
	f=f.loc[f["Circle"]!="NC"]
	f["D"]=f["Circle"].apply(lambda x: int(x.split("_")[-1]))
	f.loc[f["D"]>1,"JunFeature"]="Ins"
	f.loc[f["D"]<1,"JunFeature"]="Del"
	f.loc[f["D"]==1,"JunFeature"]="Adj"
	f["CirS"]=f["Circle"].apply(lambda x: int(x.split("_")[0]))
	f["CirE"]=f["Circle"].apply(lambda x: int(x.split("_")[1]))
	f["CirSize"]=f["CirE"]-f["CirS"]+1
	f=f.drop(["m","infor","p","D","CirS","CirE"],axis=1)
	f.to_csv(r+".datatype.tsv",index=None,sep="\t")
	m=f.groupby(["Readname"],as_index=False).filter(lambda x:len(x)>2)	
	m=m.drop_duplicates(["Readname"],keep="first")
	g=f.drop_duplicates(["Readname"],keep="first")
	print("Number of circle reads: %s"%(g.shape[0]))
	print("Largest circle size: %s"%(g["CirSize"].max()))
	print("Smallest circle size: %s"%(g["CirSize"].min()))
	print("Average circle size: %s"%(g["CirSize"].mean()))
	print("Tandem repeat>1 reads: %s" %(m.shape[0]))

getFile(res1)


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

	# Filter out all the rows where the "Circle" column has the value "NC"
	f=f.loc[f["Circle"]!="NC"]

	# Create a new column "D" which contains the value from the "Circle" column after splitting it by "_"
	# and taking the last element
	f["D"]=f["Circle"].apply(lambda x: int(x.split("_")[-1]))

	# Create a new column "JunFeature" and fill it with "Ins" if the value of "D" is greater than 1,
	# "Del" if the value of "D" is less than 1, and "Adj" if the value of "D" is equal to 1
	f.loc[f["D"]>1,"JunFeature"]="Ins"
	f.loc[f["D"]<1,"JunFeature"]="Del"
	f.loc[f["D"]==1,"JunFeature"]="Adj"

	# Create columns "CirS" and "CirE" which contain the start and end values of the "Circle" column after splitting it
	f["CirS"]=f["Circle"].apply(lambda x: int(x.split("_")[0]))
	f["CirE"]=f["Circle"].apply(lambda x: int(x.split("_")[1]))

	# Create a column "CirSize" which contains the size of the circle
	f["CirSize"]=f["CirE"]-f["CirS"]+1

	# Drop the columns "m", "infor", "p", "D", "CirS", and "CirE"
	f=f.drop(["m","infor","p","D","CirS","CirE"],axis=1)

	# Write the resulting dataframe to a new file with .datatype.tsv extension
	f.to_csv(r+".datatype.tsv",index=None,sep="\t")

	# Filter out all the rows where the same "Readname" has more than 2 occurrences
	m=f.groupby(["Readname"],as_index=False).filter(lambda x:len(x)>2)	

	# Keep only the first unique "Readname"
	m=m.drop_duplicates(["Readname"],keep="first")

	# Keep only the first unique "Readname" in the original dataframe
	g=f.drop_duplicates(["Readname"],keep="first")

	# Print the results
	print("Number of circle reads: %s"%(g.shape[0]))
	print("Largest circle size: %s

