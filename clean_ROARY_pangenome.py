### LIBS

import sys
import os
from Bio import SeqIO
import pandas
import re #this time only
from BCBio import GFF
import math

### DEFS

def make_multifastas(bits_list,ROARY_matrix,GFF_folder):
	tab_pd=pandas.pandas.read_csv(ROARY_matrix,index_col=0)
	tab_pd=tab_pd.drop(['Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate', 'Genome Fragment','Order within Fragment', 'Accessory Fragment','Accessory Order with Fragment', 'QC', 'Min group size nuc','Max group size nuc', 'Avg group size nuc'],axis=1)
	tab=tab_pd.to_dict(orient='index')
	for x in bits_list:
		outname="%s.fasta" %(x)
		outF=open(outname,"w")
		x=x.strip().split("__")[1] #this time only
		x=re.sub(r'(\d+)', '_\\1', x) #this time only
		#print tab[x]#print x
		for y in tab[x].keys():
			if str(tab[x][y])=="nan":
				continue
			fetch=tab[x][y].split()
			#print fetch
			inname="%s/%s.gff" %(GFF_folder,y)
			inF=open(inname,"r")
			for j in GFF.parse(inF):
				for z in j.features:
					if str(z.id) in fetch:
						S=str(z.extract(j).seq)
						outF.write(">%s\n%s\n" %(str(z.id),S))
			inF.close()
		outF.close()


def test_alignment(bit,IS_file,OUTFOLDER):
	cmd="blastn -task blastn -query %s.fasta -db %s -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send' | awk '$3>90' | awk '$4>=30'| cut -f1 | sort | uniq > IS_bits" %(bit,IS_file)
	os.system(cmd)
	lisF=open("IS_bits","r")
	genes_with_bits=[]
	for i in lisF.xreadlines():
		genes_with_bits.append(i.strip())
	lisF.close()
	os.system("rm IS_bits")
	inname="%s.fasta" %(bit)
	inF=open(inname,"r")
	DICT={}
	for w in SeqIO.parse(inF,"fasta"):
		DICT[str(w.id)]={}
		DICT[str(w.id)]["IS_bit"]=0
		if str(w.id) in genes_with_bits:
			DICT[str(w.id)]["IS_bit"]=1
		DICT[str(w.id)]["LENGTH"]=len(str(w.seq))
		DICT[str(w.id)]["START"]=0
		starts=["ATG","GTG","TTG"]
		if str(w.seq)[0:3] in starts:
			DICT[str(w.id)]["START"]=1
		DICT[str(w.id)]["STOP"]=0
		stops=["TAA","TAG","TGA"]
		if str(w.seq)[-3:] in stops:
			DICT[str(w.id)]["STOP"]=1
		DICT[str(w.id)]["SEQUENCE"]=str(w.seq)
	inF.close()
	DIR="%s/IS_bits_alignments" %(OUTFOLDER)
	if os.path.isdir(DIR)==False:
		os.mkdir(DIR)
	cmd="mv %s.fasta %s/IS_bits_alignments" %(bit,OUTFOLDER)
	os.system(cmd)
	return DICT

		

### INPUTS

OUTFOLDER=sys.argv[1].strip() # ISCut_GGMMYY_hhmm

ROARY_matrix=sys.argv[2].strip() # gene_presence_absence.csv

GFF_folder=sys.argv[3].strip() # used as input for ROARY

IS_file=sys.argv[4].strip() # fasta file with IS

### MAIN

bitsfilename="%s/WARNING_IS_bits_in_PANGENOME" %(OUTFOLDER)
lisF=open(bitsfilename,"r")
bits_list=[]
for i in lisF.xreadlines():
	bits_list.append(i.strip())
lisF.close()

make_multifastas(bits_list,ROARY_matrix,GFF_folder)

xlsname="%s/IS_bits_check.xlsx" %(OUTFOLDER)
writer = pandas.ExcelWriter(xlsname)
for bit in bits_list:
	DICT=test_alignment(bit,IS_file,OUTFOLDER)
	df=pandas.DataFrame.from_dict(DICT, orient="index")
	order=["LENGTH","START","STOP","IS_bit","SEQUENCE"]
	df=df[order]
	df.to_excel(writer,bit)
writer.save()
