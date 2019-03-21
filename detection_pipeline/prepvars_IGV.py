## Purpose: given a list of mosaics or denovos
## 		1. parse sites with VAF <= 0.3
##		2. prepare files for automatic IGV screenshot generation
## Usage: python prepvars_IGV.py <VARIANTS FILE> <LIST OF BAMs/CRAMs> > <IGV INPUT FILENAME>
import sys

f1 = open(sys.argv[1],'r') # candidate mosaics or denovos file
f2 = open(sys.argv[2],'r') # crams list


def parsebams(blist1):
	out = {} # dict in form id:[bam file name, bam location]
	for line in blist1:
		fullpath = line.strip()
		tmp = line.strip().split('/')
		fname = tmp[-1]
		id = fname.split('.')[0]
		out[id] = [fname, fullpath]
	blist1.close()
	return out

bams = parsebams(f2)

bamout = open('bams.txt', 'w') # bams.txt containing full bam paths
idx = {}
for line in f1:
	tmp = line.strip().split('\t')
	if tmp[0] == 'id':
		for i in xrange(len(tmp)):
			idx[tmp[i]] = i
	else:
		chr = tmp[idx['chr']]
		pos = tmp[idx['pos']]
		id = tmp[idx['id']]
		
		fid = id.split('_')[0]
		faid = fid+'-01'
		moid = fid + '-02'

		refdp, altdp = tmp[idx['refdp']], tmp[idx['altdp']]
		N = int(refdp) + int(altdp)
		vaf = float(altdp)/N

		if vaf<=0.3: ## only process sites with VAF<=0.3 to reduce manual review burden

			pbloc = ''

			if id in bams:
				pbloc = bams[id][0]
				try:
					floc = bams[faid][0]
					bamout.write(bams[faid][1]+'\n')
				except:
					floc = ''
				try:
					mloc = bams[moid][0]
					bamout.write(bams[moid][1]+'\n')
				except:
					mloc = ''
				if floc != '' or mloc !='':
					bout = ','.join([pbloc, floc, mloc])
				else:
					bout = pbloc

				print '\t'.join([chr, pos, bout])
				bamout.write(bams[id][1]+'\n')

f1.close()
bamout.close()


