import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import pandas as pd
pd.set_option("display.max_columns",40)

# get file name and output name from the command line arguments
res1=sys.argv[1]
name=sys.argv[2]

# Function to determine the type of circle
def circleType(x):
	# Split the string by "_" and extract each component as integer
	cirCle_S=int(x.split("_")[0])
	cirCle_E=int(x.split("_")[1])
	cirCle_D=int(x.split("_")[2])
	cirCle_L=int(x.split("_")[-1])
	# Check the conditions to determine the type of circle
	if cirCle_S <100 and cirCle_E>cirCle_L-100 and cirCle_D<=-100:
		return "1End_FL"
	elif cirCle_S <100 and cirCle_E>cirCle_L-100 and cirCle_D>-100:
		return "2Ends_FL"
	elif cirCle_S <100 or cirCle_E>cirCle_L-100:
		return "1End_Frg"
	else:
		return "noEnd_Frg"

# Function to get the file, process and output the data
def getFile(r,Out_name):
	# Read the input file
	f=pd.read_table(r)
	# Filter the data for "Circle" column which are not "NC"
	f=f.loc[f["Circle"]!="NC"]
	# Drop the unnecessary columns "m", "infor", "p"
	f=f.drop(["m","infor","p"],axis=1)
	# Concatenate "Circle" and "Reflen" columns to get "circileInfor"
	f["circileInfor"]=f["Circle"]+"_"+f["Reflen"].apply(str)
	# Drop the duplicate rows based on "Readname" column, keeping the first one
	f=f.drop_duplicates(["Readname"],keep="first")
	# Apply the function "circleType" to the "Circle" column to get the type of circle
	f["Type"]=f["Circle"].apply(lambda x : circleType(x))
	# Create a dictionary with "Readname" as key and "Type" as value
	d=dict(zip(list(f["Readname"]),list(f["Type"])))

	# Read the input file again
	f_=pd.read_table(r)
	# Filter the data for "Circle" column which are not "NC"
	f_=f_.loc[f_["Circle"]!="NC"]
	# Add the "Type" column based on the dictionary created above
	f_["Type"]=f_["Readname"].apply(lambda x: d[x])
	# Write the output to the file "Out_name_data.txt"
	f_.to_csv(Out_name + "_data.txt", index=None, sep="\t")
	g = f.groupby(["Type"], as_index=False).count()
	g = g[["Type", "Readname"]]
	print(g)
	g.to_csv(Out_name, index=None, sep="\t", header=None)

getFile(res1, name)	
