#### LIBRARIES
import sys
import os
import glob
from Bio import SeqIO
import pandas
import gzip
import time

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
		#genome_hash[name]["number"]=int(sum(1 for line in gzip.open(READS[0]))/2)
	return genome_hash


def bowtieDB(db):
	cmd="bowtie2-build %s %s" %(db,db)
	os.system(cmd)


def run_bowtie2(name,readF,readR,db,threads):
	cmd="bowtie2 -x %s -1 %s -2 %s --very-sensitive-local --no-unal -S %s.sam -p %i" %(db,readF,readR,name,threads)
	os.system(cmd)

def parse_ISAba(File):
	hash={}
	inF=open(File,"r")
	for i in SeqIO.parse(inF,"fasta"):
		fam=str(i.description).strip().split()
		fam=fam[1]+"_"+fam[2]
		fam=fam.replace("/","_").replace(":","")
		hash[str(i.id)]=[fam,len(str(i.seq))]
	
	inF.close()
	return hash

def parse_sam(sfile,hISA,hRES,h,NORMlengths):
	samfile=open(sfile,"r")
	Lclasses=[]
	hreads={}
	hreads["total"]=0.0
	NORMhash={}
	for y in hISA.keys():
		#print hISA[y]
		Lclasses.append(hISA[y][0])
	Lclasses=list(set(Lclasses))
	for x in Lclasses:
		hreads[x]=0.0
	for i in samfile.xreadlines():
		if i.strip()[0]!="@":
			i=i.strip().split()
			hit=i[2].strip()
			#print hit
			if hit in hISA.keys():
				hreads["total"]+=1/float(hISA[hit][1])
				Fhit=hISA[hit][0]
				hreads[Fhit]+=1/float(hISA[hit][1])
			else:
				if hit not in NORMhash.keys():
					NORMhash[hit]=0
				NORMhash[hit]+=1
	#print NORMhash
	#print hISA.keys()
	#print hreads

	VALS=[]
	for i in NORMhash.keys():
		val=float(NORMhash[i])/float(NORMlengths[i])
		VALS.append(val)
		#print "\n\n\n"
		#print VALS
		#print "\n\n\n"
	#print VALS
	NORM=sum(VALS)/len(VALS)
	for j in hreads.keys():
		index="%s_norm" %(j)
		hreads[index]=float(hreads[j])/NORM
	
	#print hreads
	
	hRES[h]=hreads
	samfile.close()
	return hRES




def start_bowtie(HH,refname,bowtie_threads,hRES,hISA,OUTFOLDER,h,NORMlengths):
	run_bowtie2(h,HH[h]["R1"],HH[h]["R2"],refname,bowtie_threads)
	sfile="%s.sam" %(h)
	hRES=parse_sam(sfile,hISA,hRES,h,NORMlengths)
	cmd="rm %s.sam" %(h)
	os.system(cmd)
	return hRES




def build_ref_4_normalization(refname,coregen):
	NORMlengths={}
	badlist=[]
	cmd="makeblastdb -in %s -dbtype nucl"	%(coregen)
	os.system(cmd)
	cmd="blastn -task blastn -query %s -db %s -outfmt '6 sseqid' | sort | uniq > seqs_to_remove" %(refname,coregen)
	os.system(cmd)
	inF=open("seqs_to_remove","r")
	f=inF.read().strip().split("\n")
	for i in f:
		badlist.append(i.strip())
	inF.close()
	inF=open(coregen,"r")
	tot=0
	outF=open("QUANT_REFERENCE","w")
	for x in SeqIO.parse(inF,"fasta"):
		tot+=1
		if str(x.id) not in badlist:
			SeqIO.write(x,outF,"fasta")
			NORMlengths[str(x.id)]=len(str(x.seq))

	inF.close()
	outF.close()
	cmd="cat %s >> QUANT_REFERENCE" %(refname)
	os.system(cmd)
	bowtieDB("QUANT_REFERENCE")
	print "Total core genes used for normalization: %i" %(tot-len(badlist))
	return NORMlengths

##### MAIN
		
inF=open(sys.argv[1].strip(),"r") #org_list
HH=parse_org_list(inF)
inF.close()

refname=sys.argv[2].strip() #ref_file "ISAba.ffn"

hash_ISAba=parse_ISAba(refname)

num_threads=int(sys.argv[3].strip())

coregen=sys.argv[4].strip()


NORMlengths=build_ref_4_normalization(refname,coregen)




bowtie_threads=num_threads

hashRES={}

OUTFOLDER="MAP_IS_%s" %(str(time.strftime("%d%m%y_%H%M")))
os.mkdir(OUTFOLDER)

hitsfolder="%s/HITS" %(OUTFOLDER)
os.mkdir(hitsfolder)



for h in HH.keys():
	hashRES=start_bowtie(HH, "QUANT_REFERENCE", bowtie_threads,hashRES,hash_ISAba,OUTFOLDER,h,NORMlengths)



#print hashRES

##### MAIN OUTPUT #####

hashRES_pd=pandas.DataFrame.from_dict(hashRES, orient='columns')
outFILE="%s/Report_mapping.txt" %(OUTFOLDER)
hashRES_pd.to_csv(outFILE, sep='\t', encoding='utf-8')

"""
#### CLEAN UP
cmd="rm %s.*.bt2 %s.*.bt2" %(refname,PANG)
os.system(cmd)

"""

