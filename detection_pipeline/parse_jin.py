## Usage: python parse_jin.py <FILTERED DENOVOS> ~/mosaic/jin_et_al/jin_all_samples.TRIOS_ONLY.txt > <FILTERED DENOVOS ONLY FOR SAMPLES IN JIN ET AL>
## Purpose: Parse only variants that are called in samples found in Jin et al. 2017 publication cohort
## Dependencies: list of IDs 
import sys

f = open(sys.argv[1],'r') # samtools calls
f2 = open(sys.argv[2],'r') # jin ids

j = []
for line in f2:
	tmp = line.strip()
	j.append(tmp)
f2.close()

tot = 0
ct = 0
for line in f:
	tmp = line.strip().split('\t')
	if tmp[0] == "id":
		print '\t'.join(tmp)
	else:
		tot += 1
		if tmp[0] in j:
			ct += 1
			print '\t'.join(tmp)
f.close()
#print "Total samtools calls: %d  Calls in Jin samples: %d"%(tot, ct)