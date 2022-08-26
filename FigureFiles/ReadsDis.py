import pandas as pd
import sys
import os

pd.set_option("display.max_column",40)


def ReadDis(file1,file2,file3):
	f=pd.read_table(file1)
	f=f.drop_duplicates(["Readname"],keep="first")
	print(f.shape[0])

	f=pd.read_table(file2)
	f=f.drop_duplicates(["Readname"],keep="first")
	print(f.shape[0])
	
	f=pd.read_table(file3)
	f=f.drop_duplicates(["Readname"],keep="first")
	print(f.shape[0])
	f=f.loc[f["Circle"]!="NC"]
	print(f.shape[0])





file1=sys.argv[1]
file2=sys.argv[2]
file3=sys.argv[3]

bam=sys.argv[4]

ReadDis(file1,file2,file3)

samtools="samtools stats %s | head -n 30 > %s"%(bam,bam+".stats")
os.system(samtools)

