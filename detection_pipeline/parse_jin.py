## Usage: python parse_jin.py <FILTERED DENOVOS> ~/mosaic/jin_et_al/jin_all_samples.TRIOS_ONLY.txt > <FILTERED DENOVOS ONLY FOR SAMPLES IN JIN ET AL>
## Purpose: Parse only variants that are called in samples found in Jin et al. 2017 publication cohort
## Dependencies: list of IDs 
import sys

f = open(sys.argv[1],'r') # samtools calls

j = [] # list of ids in Jin 2017
with open(sys.argv[2],'r') as f2: # jin_2530_ids.txt
	for line in f2:
		tmp = line.strip()
		j.append(tmp)


tot = 0
ct = 0
with open(sys.argv[1],'r') as f:
	for line in f:
		tmp = line.strip().split('\t')
		if tmp[0] == "id":
			print '\t'.join(tmp)
		else:
			tot += 1
			if tmp[0].split('_')[0] in j:
				ct += 1
				print '\t'.join(tmp)

#print "Total samtools calls: %d  Calls in Jin samples: %d"%(tot, ct)