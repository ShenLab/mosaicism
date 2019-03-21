import sys

idx = {}
filtd = {'MD':[], 'LCR':[], 'mappability':[], 'segdup':[], 'readposbias':[], 'strandbias':[], 'cohortAF':[]}

tot = 0
survive = 0
with open(sys.argv[1],'r') as f: # *.filtered.txt
	for line in f:
		tmp = line.strip().split('\t')
		if tmp[0] == "id":
			for i in xrange(len(tmp)):
				idx[tmp[i]] = i
		else:
			tot += 1
			filt = tmp[idx['filter']]
			if 'exdp' in filt:
				filtd['MD'].append('\t'.join(tmp))
			if 'LCR' in filt:
				filtd['LCR'].append('\t'.join(tmp))
			if 'mappability' in filt:
				filtd['mappability'].append('\t'.join(tmp))
			if 'segdup' in filt:
				filtd['segdup'].append('\t'.join(tmp))
			if 'read' in filt:
				filtd['readposbias'].append('\t'.join(tmp))
			if 'strand' in filt:
				filtd['strandbias'].append('\t'.join(tmp))
			if 'cohort' in filt:
				filtd['cohortAF'].append('\t'.join(tmp))
			if filt == ".":
				survive += 1


print '\t' + "Total: " + '\t' + str(tot)
for k in filtd.keys():
	print '\t' + k + '\t' + str(len(filtd[k]))
print '\t' + "PASS: " + '\t' + str(survive)

