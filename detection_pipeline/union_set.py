import sys

f1 = open(sys.argv[1],'r') # ~/mosaic/pcgc_pipeline/pcgc0322/ADfile.jin0322.txt
f2 = open(sys.argv[2],'r') # ~/mosaic/pcgc_pipeline/pcgc0322/sam0322.denovo.txt



idx1 = {}
jind = {}
for line in f1:
	tmp = line.strip().split('\t')
	if tmp[0] == "id":
		head = '\t'.join(tmp) + '\t' + 'set'
		for i in xrange(len(tmp)):
			idx1[tmp[i]] = i
	else:
		id, chr, pos, ref, alt = tmp[idx1['id']], tmp[idx1['chr']], tmp[idx1['pos']], tmp[idx1['ref']], tmp[idx1['alt']]
		key = '|'.join([id, chr, pos, ref, alt])
		jind[key] = '\t'.join(tmp)
		#if tmp[idx1['filter']] == ".":
		#	jind[key] = '\t'.join(tmp)
f1.close()

idx = {}
samd = {}
for line in f2:
	tmp = line.strip().split('\t')
	if tmp[0] == "id":
		for i in xrange(len(tmp)):
			idx[tmp[i]] = i
		samhead = tmp
	else:
		id, chr, pos, ref, alt = tmp[idx['id']], tmp[idx['chr']], tmp[idx['pos']], tmp[idx['ref']], tmp[idx['alt']]
		key = '|'.join([id, chr, pos, ref, alt])
		#out = [ tmp[idx['id']], tmp[idx['gene']], tmp[idx['chr']], tmp[idx['pos']], tmp[idx['ref']], tmp[idx['alt']], tmp[idx['refdp']],tmp[idx['altdp']],tmp[idx['aachange']],tmp[idx['maxfreq']],tmp[idx['vfunc']], tmp[idx['vtype']],tmp[idx['cadd']],tmp[idx['cadd_flag']],tmp[idx['meta']],tmp[idx['meta_flag']],]
		out = tmp[0:samhead.index('cohort_AF')+1]
		samd[key] = '\t'.join(out)
f2.close()


print "Jin: %d total dnSNVs"%(len(jind.keys()))
print "samtools: %d total dnSNVs"%(len(samd.keys()))
union = list(set(jind.keys()) | set(samd.keys()))
print "union: %d"%(len(union))
overlap = list(set(jind.keys()) & set(samd.keys()))
print "intersection: %d"%(len(overlap))
uniq_jin = list(set(jind.keys()) - set(overlap))
uniq_sam = list(set(samd.keys()) - set(overlap))
print "Jin-only: %d    samtools-only: %d"%(len(uniq_jin), len(uniq_sam))

outf = open(sys.argv[3], 'w')
outf.write(head+'\n')
for k in union:
	if k in overlap:
		src = 'overlap'
		out = samd[k] + '\t' + src 
	elif k in uniq_jin:
		src = 'jin_only'
		out = jind[k] + '\t' + src 
	elif k in uniq_sam:
		src = 'sam_only'
		out = samd[k] + '\t' + src
	outf.write(out + '\n')
outf.close()
print 'Done -- see %s'%(sys.argv[3])


