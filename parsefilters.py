## Usage: python parsefilters.py candidates.denovo.txt rfilt.denovo.txt candidates.denovo.dp4.fishers.txt > candidates.denovo.filtered.txt
## Purpose: combines repeat filters with fisher's filters and extreme depth filters
## Output: candidates file cols + ';' separated filter hits.  '.' means survived
import sys

## parse intersectBed output
hitd = dict() # dictionary of repeat filter hit {chr|pos : [filters hit]}
with open(sys.argv[2],'r') as f2: #intersectBed output rfilt.*.txt
	for line in f2:
		tmp = line.strip().split('\t')
		chr = tmp[0][3:]
		pos = tmp[1]
		key = "|".join([chr,pos])

		filt = tmp[3]
		hit = ''
		if "mappability" in filt: #################### get score
			score = tmp[-1]
			hit = "mappability_"+score
		elif "LCR" in filt:
			hit = "LCR"
		elif "segdup" in filt: #################### get score
			score = tmp[-1]
			hit = "segdup_"+score

		if not key in hitd:
			hitd[key] = [hit]
		else:
			hitd[key].append(hit)

if len(sys.argv) > 2:
	sbd = {} # strand bias dictionary
	with open(sys.argv[3],'r') as f3: # strand bias annotation *.sb.txt
		for line in f3:
			tmp = line.strip().split('\t')
			if tmp[0] == "id":
				idx2 = {col:index for index, col in enumerate(tmp)}
			else:
				id, chr, pos, ref, alt = tmp[idx2['id']], tmp[idx2['chr']], tmp[idx2['pos']], tmp[idx2['ref']], tmp[idx2['alt']]
				key = '|'.join([id, chr, pos])
				sbflag, sbor, sbp = tmp[idx2['strand_bias_flag']], tmp[idx2['strand_bias_or']], tmp[idx2['strand_bias_p']]
				sbd[key] = [sbflag, sbor, sbp]



## parse candidates file and generate output 
with open(sys.argv[1],'r') as f1: # candidates file
	for line in f1:
		tmp = line.strip().split('\t')
		if tmp[0] == "id":
			head = '\t'.join(tmp)+"\t"+ 'filter'
			print head
			idx = {col:index for index, col in enumerate(tmp)}
		else:
			chr = tmp[idx['chr']]
			fullchr = "chr"+tmp[idx['chr']]
			pos = tmp[idx['pos']]
			id = tmp[idx['id']]
			
			key1 = "|".join([chr,pos]) # for hitd; repeat filters and dbsnp
			key2 = "|".join([id,chr,pos]) # for hitd2; strand filters

			rep, strand = ['.'], ['.']
			if key1 in hitd:
				rep = hitd[key1]
			if key2 in sbd:
				sbflag = sbd[key2][0]
				if sbflag == "1":
					strand.append('strand_bias')

			out = rep + strand

			if 'cohort_AF' in idx.keys():
				if float(tmp[idx['cohort_AF']]) >= 0.01:
					out.append('cohortAF_0.01')

			outstring = '\t'.join(tmp)+"\t"+ ';'.join(list(set(out)))
			print outstring
			




