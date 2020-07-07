import sys
import gzip

F1=gzip.open(sys.argv[1].strip(),"r")
F2=gzip.open(sys.argv[2].strip(),"r")

c=0
d=0

for l1 in F1.readlines():
	c+=1
	l2=F2.readline()
	if l1!=l2:
		d+=1

print float(d)/float(c)

F1.close()
F2.close()
	
