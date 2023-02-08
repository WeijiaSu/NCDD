import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import pandas as pd
pd.set_option("display.max_columns",40)

# Get file name as an argument
res1=sys.argv[1]

def getFile(r):
    # Read the input file
    f=pd.read_table(r)
    # Filter the rows where the "Circle" column has a value of "NC"
    f=f.loc[f["Circle"]!="NC"]
    # Add a new column "D" to the data frame which is the last part of the "Circle" column value split by "_"
    f["D"]=f["Circle"].apply(lambda x: int(x.split("_")[-1]))
    # Add a new column "JunFeature" to the data frame based on the value of the "D" column
    f.loc[f["D"]>1,"JunFeature"]="Ins"
    f.loc[f["D"]<1,"JunFeature"]="Del"
    f.loc[f["D"]==1,"JunFeature"]="Adj"
    # Add new columns "CirS" and "CirE" to the data frame which are the first and second part of the "Circle" column value split by "_"
    f["CirS"]=f["Circle"].apply(lambda x: int(x.split("_")[0]))
    f["CirE"]=f["Circle"].apply(lambda x: int(x.split("_")[1]))
    # Add a new column "CirSize" to the data frame which is the difference between "CirE" and "CirS" plus 1
    f["CirSize"]=f["CirE"]-f["CirS"]+1
    # Remove the columns "m", "infor", "p", "D", "CirS", "CirE" from the data frame
    f=f.drop(["m","infor","p","D","CirS","CirE"],axis=1)
    # Write the modified data frame to a new file with a specified format
    f.to_csv(r+".datatype.tsv",index=None,sep="\t")
    # Filter the data frame to only include rows where the "Readname" column has more than 2 unique values
    m=f.groupby(["Readname"],as_index=False).filter(lambda x:len(x)>2)    
    # Remove duplicate rows based on the "Readname" column and keep only the first one
    m=m.drop_duplicates(["Readname"],keep="first")
    # Remove duplicate rows based on the "Readname" column and keep only the first one
    g=f.drop_duplicates(["Readname"],keep="first")
    # Print out various statistics about the filtered data frame
    print("Number of circle reads: %s"%(g.shape[0]))
	print("Largest circle size: %s"%(g["CirSize"].max()))
    print("Small



