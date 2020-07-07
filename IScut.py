#### LIBRARIES
import sys
import os
import glob
from Bio import SeqIO
import pandas
import gzip
import time
import multiprocessing
from functools import partial

#### DEFS

def parse_org_list(FILE):
	genome_hash={}

	for i in FILE.xreadlines():
		i=i.strip().split()
		name=i[0]
		genome_hash[name]={}
		reads="%s*" %(i[1])
		READS=glob.glob(reads)
		genome_hash[name]["R1"]=READS[0]
		genome_hash[name]["R2"]=READS[1]
		genome_hash[name]["number"]=int(sum(1 for line in gzip.open(READS[0]))/2)
	return genome_hash

def trim_reads(h,HH,read_len):
	Read1=HH[h]["R1"]
	Read2=HH[h]["R2"]

	inF=gzip.open(Read1,"r")
	outname="%s_R1_trimmed%i.fastq.gz" %(h,read_len)
	outF=gzip.open(outname, "w")
	for record in SeqIO.parse(inF, "fastq"):
		if len(str(record.seq))>read_len:
			trimmed_record=record[:read_len]
			SeqIO.write(trimmed_record, outF, "fastq")
		else:
			SeqIO.write(record,outF,"fastq")
	outF.close()
	inF.close()
	inF=gzip.open(Read2,"r")
	outname="%s_R2_trimmed%i.fastq.gz" %(h,read_len)
	outF=gzip.open(outname, "w")
	for record in SeqIO.parse(inF, "fastq"):
		if len(str(record.seq))>read_len:
			trimmed_record=record[:read_len]
			SeqIO.write(trimmed_record, outF, "fastq")
		else:
			SeqIO.write(record,outF,"fastq")
	outF.close()
	inF.close()


def parallel_trim(data_list,threads,HASH,LEN):
    pool = multiprocessing.Pool(processes=threads)
    TTT=partial(trim_reads, HH=HASH, read_len=LEN) #######################
    pool.map(TTT, data_list)



def bowtieDB(db):
	cmd="bowtie2-build %s %s" %(db,db)
	os.system(cmd)



def parse_ISAba(File):
	hash={}
	inF=open(File,"r")
	for i in SeqIO.parse(inF,"fasta"):
		fam=str(i.description).strip().split()
		fam=fam[1]+"_"+fam[2]
		fam=fam.replace("/","_").replace(":","")
		hash[str(i.id)]=fam

	inF.close()
	return hash

def remove_IS_from_PANGENOME(ISname,PANG,OUTFOLDER):
	PANG2="%s_no_IS.fasta" %(PANG)
	cmd="makeblastdb -in %s -dbtype nucl" %(ISname)
	os.system(cmd)
	cmd="blastn -task blastn -query %s -db %s -outfmt '6 qseqid sseqid pident length qlen slen' | awk '$3>60' | awk '$5/$4>0.6' | awk '$5/$4<1.4' | cut -f1 | sort | uniq > REMOVED_IS_list" %(PANG,ISname)
	os.system(cmd)
	lisF=open("REMOVED_IS_list","r")
	bad_list=[]
	for i in lisF.xreadlines():
		bad_list.append(i.strip())
	lisF.close()

	cmd="blastn -task blastn -query %s -db %s -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send' | awk '$3>90' | awk '$4>=30'| cut -f1 | sort | uniq > IS_bits" %(PANG,ISname)
	os.system(cmd)
	lisF=open("IS_bits","r")
	bits_list=[]
	for i in lisF.xreadlines():
		bits_list.append(i.strip())
	lisF.close()
	bits_list2=[]

	for ii in bits_list:
		if ii not in bad_list:
			bits_list2.append(ii)
	bits_list=[]
	bits_list=bits_list2
	
	error_message="WARNING: %s IS bits were found in your PANGENOME, please see the WARNING file in the results folder!\nThese may alter your results... or just be your result, you decide" %(len(bits_list))
	print error_message

	
	os.system("rm IS_bits")
	bitsfilename="%s/WARNING_IS_bits_in_PANGENOME" %(OUTFOLDER)
	bitsF=open(bitsfilename,"w")
	for ii in bits_list:
		bitsF.write("%s\n" %(ii))
	bitsF.close()
	
	outF=open(PANG2,"w")
	inF=open(PANG,"r")
	for j in SeqIO.parse(inF,"fasta"):
		if str(j.id) not in bad_list:
			SeqIO.write(j,outF,"fasta")
	outF.close()
	inF.close()
	cmd="mv REMOVED_IS_list %s" %(OUTFOLDER)
	os.system(cmd)	
	return PANG2
	



def run_bowtie2(name,db,threads,read_len): ################## WITH UNALS
	cmd="bowtie2 -x %s -1 %s_R1_trimmed%i.fastq.gz -2 %s_R2_trimmed%i.fastq.gz --very-sensitive-local -S %s.sam -p %i -N 2" %(db,name,read_len,name,read_len,name,threads)
	os.system(cmd)


def simplify_sam(h,ISAbaLIST):
	cmd="grep -v '^@' %s.sam > samlike.tmp" %(h)
	os.system(cmd)
	inF=open("samlike.tmp","r")
	outname="%s.simple_sam" %(h)
	outF=open(outname,"w")
	while True:
		R1=inF.readline()
		if not R1:
			break
		L1=R1.strip().split()
		R2=inF.readline()
		if not R2:
			break
		L2=R2.strip().split()
		if L1[0]!=L2[0]:
			print "Read pair ERROR!!!"
			break
		if L1[2]!=L2[2]:
			if ((L1[2] in ISAbaLIST) or (L2[2] in ISAbaLIST)) and ((L1[2] not in ISAbaLIST) or (L2[2] not in ISAbaLIST)):
				if (L1[2]!="*") and (L2[2]!="*"):
					outF.write(R1)
					outF.write(R2)
	inF.close()
	outF.close()
	os.system("rm samlike.tmp")


def Sam2Matches(h,ISAbaHASH):
	inname="%s.simple_sam" %(h)
	inF=open(inname,"r")
	outname="%s.matches" %(h)
	outF=open(outname,"w")
	outF.write("# ISAba\tISAba_pos\tISAba_read_FR\tISAba_RevComp\tGene\tGene_pos\tGene_RevComp\n")
	for i in inF.xreadlines():
		i=i.strip().split()
		if i[2] in ISAbaHASH.keys():
			IS=ISAbaHASH[i[2]]
			ISpos=i[3]
			gene=i[6]
			genepos=i[7]
			
			flag=bin(int(i[1].strip()))
			flag=str(flag).split("b")[1]

			FR="F"
			if len(flag)>=8:
				if flag[-8]=="1":
					FR="R"

			geneCOMPL=False
			if len(flag)>=6:
				if flag[-6]=="1":
					geneCOMPL=True

			ISABACOMPL=False
			if len(flag)>=5:
				if flag[-5]=="1":
					ISABACOMPL=True

			outF.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(IS,ISpos,FR,ISABACOMPL,gene,genepos,geneCOMPL,flag,i[1]))
	inF.close()
	outF.close()

def get_interrupted(h,OUTFOLDER):
	hash={}
	inname="%s.matches" %(h)
	inF=open(inname,"r")
	for i in inF.xreadlines():
		if i.strip()[0]=="#":
			pass
		else:
			i=i.strip().split()
			ISname=i[0]
			genename=i[4]
			if ISname not in hash.keys():
				hash[ISname]={}
			if genename not in hash[ISname].keys():
				hash[ISname][genename]=[]
			OR="+"
			if i[3]=="True":
				OR="-"
			SS="%s%s" %(i[2],OR)
			hash[ISname][genename].append(SS)
	inF.close()
	for x in hash.keys():
		for y in hash[x].keys():
			ll=hash[x][y]
			if (("F+" in ll) and ("R-" in ll)) or (("F-" in ll) and ("R+" in ll)):     ######## CHECK CHECK CHECK
				hash[x][y]="X"
			else:
				hash[x][y]="?"
	DF=pandas.DataFrame.from_dict(hash, orient='columns')
	DF = DF.fillna('')
	outname="%s/by_GENOME/%s.csv" %(OUTFOLDER,h)
	DF.to_csv(outname, sep='\t', encoding='utf-8')



#################
##### MAIN ######
#################

### read inputs ###

inF=open(sys.argv[1].strip(),"r") #org_list
HH=parse_org_list(inF)
inF.close()

ISname=sys.argv[2].strip() #ISABA_file "ISAba.ffn"
hash_ISAba=parse_ISAba(ISname)

PANG=sys.argv[3].strip() #pangenome file

bowtie_threads=int(sys.argv[4].strip()) # even number

read_len=int(sys.argv[5].strip())



### prepare output folder ###

OUTFOLDER="ISCut_%s" %(str(time.strftime("%d%m%y_%H%M")))
os.mkdir(OUTFOLDER)

hitsfolder="%s/matches" %(OUTFOLDER)
os.mkdir(hitsfolder)

ISfolder="%s/by_IS" %(OUTFOLDER)
os.mkdir(ISfolder)

GENOMEfolder="%s/by_GENOME" %(OUTFOLDER)
os.mkdir(GENOMEfolder)

### do stuff ###

PANG2=remove_IS_from_PANGENOME(ISname,PANG,OUTFOLDER)

cmd="cat %s %s > OneReference" %(ISname,PANG2)
os.system(cmd)

cmd="rm %s" %(PANG2)
os.system(cmd)

bowtieDB("OneReference")




parallel_trim(HH.keys(),bowtie_threads,HH,read_len)



for h in HH.keys():
	run_bowtie2(h,"OneReference",bowtie_threads,read_len)
	cmd="rm %s_R1_trimmed%i.fastq.gz %s_R2_trimmed%i.fastq.gz" %(h,read_len,h,read_len)
	os.system(cmd)
	simplify_sam(h,hash_ISAba.keys())

	cmd="rm %s.sam" %(h)
	os.system(cmd)



os.system("rm *.bt2")
cmd="mv OneReference %s" %(OUTFOLDER)
os.system(cmd)


for h in HH.keys():
	Sam2Matches(h,hash_ISAba)

os.system("rm *.simple_sam")


target_list=[]
for h in HH.keys():
	inname="%s.matches" %(h)
	inF=open(inname,"r")
	for i in inF.xreadlines():
		if i.strip()[0]=="#":
			pass
		else:
			t=i.strip().split()[4]
			target_list.append(t)
	inF.close()

target_list=list(set(target_list))



for h in HH.keys():
	get_interrupted(h,OUTFOLDER)

cmd="mv *.matches %s/matches/" %(OUTFOLDER)
os.system(cmd)


#OUTFOLDER="ISCut_260318_1157" #FOR_TEST

### results in the ISAba point of view ###


IS_hash={}
glob_hash={}
for h in HH.keys():
	inname="%s/by_GENOME/%s.csv" %(OUTFOLDER,h)
	tab_pd=pandas.pandas.read_csv(inname, sep='\t',index_col=0)
	tab=tab_pd.to_dict()
	for i in tab.keys():
		if i not in IS_hash.keys():
			IS_hash[i]={}
		IS_hash[i][h]=tab[i]


	glob_hash[h]={}
	t2=tab_pd.to_dict(orient='index')
	for j in t2.keys():
		LIST=[]
		for z in t2[j].keys():
			LIST.append(t2[j][z])
		if "X" in LIST:
			glob_hash[h][j]="X"
		elif "?" in LIST:
			glob_hash[h][j]="?"
		else:
			glob_hash[h][j]=""	

for x in IS_hash.keys():
	HHH=IS_hash[x]
	DF=pandas.DataFrame.from_dict(HHH, orient='columns')
	DF = DF.fillna('')
	outname="%s/%s.csv" %(ISfolder,x)
	DF.to_csv(outname, sep='\t', encoding='utf-8')

### global result ###

DF=pandas.DataFrame.from_dict(glob_hash, orient='columns')
DF = DF.fillna('')
outname="%s/SUMMARY.csv" %(OUTFOLDER)
DF.to_csv(outname, sep='\t', encoding='utf-8')

