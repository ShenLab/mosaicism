## Usage: python check_cohortAF.py <cohort count table> <denovos> <cohortsize> > <denovos>.mkcont.txt
import sys

cohortd = {} # {chr:pos:ref:alt : n_carriers}
with open(sys.argv[1],'r') as f1: # pcgc_vcf_cohort_ct.txt
	for line in f1:
		tmp = line.strip().split('\t')
		if tmp[0] == 'KEY':
			idx = {col:index for index, col in enumerate(tmp)}
		else:
			key, n = tmp[idx['KEY']], tmp[idx['n_carriers']]
			cohortd[key] = int(n)

cohortsize = float(sys.argv[3])

with open(sys.argv[2],'r') as f2: # denovo_filt.txt
	for line in f2:
		tmp = line.strip().split('\t')
		if tmp[0] == 'id':
			idx = {col:index for index, col in enumerate(tmp)}
			print '\t'.join(tmp) + '\t' + 'germ.carriers' + '\t' + 'cohort_AF'
		else:
			chr, pos, ref, alt = tmp[idx['chr']], tmp[idx['pos']], tmp[idx['ref']], tmp[idx['alt']]
			key = ':'.join([chr, pos, ref, alt])
			# get number of variant carriers nalt>=2 in cohort
			try:
				n = cohortd[key]
				af = n/cohortsize
				#print '\t'.join(tmp) + '\t' + str(n) + '\t' + str(af)
			except:
				n = '.'
				af = '0.0'


			print '\t'.join(tmp) + '\t' + str(n) + '\t' + str(af)