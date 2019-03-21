## Purpose: given a list of trio VCF files, call denovos
## Usage: python tissue.trio_to_denovo.py <VCF LIST> <OUTPUT FILENAME>
import sys

vcflist = open(sys.argv[1], 'r') # list of vcf files with full path
outname = sys.argv[2]

incompletetriof = open('tissue0130.incomplete_trios.txt', 'w')
zeroctf = open('tissue0130.zerocount_trios.txt', 'w')
## Function to parse samtools vcf files
## (1) get sample col idx (2) get sample AD (3) output sites where pb AD>=6 and both parents 0
def parsesam(vcffile):
	idx = {}
	sampid = vcffile.split('/')[-1].split('.')[0]
	faid = sampid.split('_')[0]+'-01'
	moid = sampid.split('_')[0]+'-02'
	v = open(vcffile, 'r')
	out = [] # holds strings containing each variant site for a given vcf
	for line in v:
		tmp = line.strip().split('\t')
		# parse header line + get col idx for each sample
		if not line.startswith("##"):
			if tmp[0] == "#CHROM":
				for i in xrange(len(tmp)): # iterate over header line
					vcfhead = tmp
					idx[tmp[i]] = i

			# parse each site
			else:
				pbid = sampid
				pbidx = idx[sampid]
				faidx = idx[faid]
				moidx = idx[moid]
				## handle incomplete trios
				try:
					if not (tmp[pbidx]=="." or tmp[faidx]=="." or tmp[moidx]=="."):
						testidxs = [pbidx, faidx, moidx]
				except:
					print "##### ERROR Incomplete trio: "+vcffile
					incompletetriof.write(vcffile+'\n')
					break

				chr, pos, ref, alt = tmp[idx['#CHROM']], tmp[idx['POS']], tmp[idx['REF']], tmp[idx['ALT']]
				tmp[idx['ID']] = pbid # for keeping track of samples when running on multiple vcfs
				
				stopidx = idx['FORMAT']
				## ignore sites for which we are missing information from proband or either parent
				if not (tmp[pbidx]=="." or tmp[faidx]=="." or tmp[moidx]=="."):
					## parse FORMAT section and create dict in form {GT:val, DP:val, AD:val, RO:val, QR:val, AO:val, QA:val}
					pbgeno = dict(zip(tmp[idx['FORMAT']].split(':'), tmp[pbidx].split(':')))
					fageno = dict(zip(tmp[idx['FORMAT']].split(':'), tmp[faidx].split(':')))
					mogeno = dict(zip(tmp[idx['FORMAT']].split(':'), tmp[moidx].split(':')))

					pb_ad = pbgeno['AD']
					fa_ad = fageno['AD']
					mo_ad = mogeno['AD']
					fa_dp = sum(map(int, fa_ad.split(',')))
					mo_dp = sum(map(int, mo_ad.split(',')))
					pb_altidx, fa_altidx, mo_altidx = 1, 1, 1
					## handle multiallelic sites
					if (len(pb_ad.split(','))>2 or len(fa_ad.split(','))>2 or len(mo_ad.split(','))>2):
						#print "pb: %s ; fa: %s ; mo: %s"%(pb_ad, fa_ad, mo_ad)
						# newalt_idx =  map(int, ad.split(',')[1:]).index(max(map(int, ad.split(',')[1:])))
						pb_altidx = map(int, pb_ad.split(',')[1:]).index(max(map(int, pb_ad.split(',')[1:]))) + 1
						fa_altidx = map(int, fa_ad.split(',')[1:]).index(max(map(int, fa_ad.split(',')[1:]))) + 1
						mo_altidx = map(int, mo_ad.split(',')[1:]).index(max(map(int, mo_ad.split(',')[1:]))) + 1
					## denovo = pb Nalt >= 6, parents Nalt==0, parents DP>=10
					if (int(pb_ad.split(',')[pb_altidx]) >= 6 and int(fa_ad.split(',')[fa_altidx]) == 0 and int(mo_ad.split(',')[mo_altidx]) == 0):
						if(fa_dp >= 10 and mo_dp >= 10):
							outstring = '\t'.join(tmp[:stopidx+1]) + '\t' + tmp[pbidx] + '\t' + tmp[faidx] + '\t' + tmp[moidx]
							out.append(outstring)
	v.close()
	return(out)
outf = open(outname, 'w')
total = 60 ## total number of trio vcfs to process
k = 0
head = ['#CHROM', 'POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','PROBAND','FATHER','MOTHER']
outf.write('\t'.join(head) + '\n')
progf = open('tissue.progress_log.txt', 'w')
for line in vcflist:
	k += 1
	vcf = line.strip()
	fname = vcf.split('/')[-1]
	tmp_prog = '( %d / %d ) %s processing ... '%(k, total, fname)
        progf.write(tmp_prog + '\n')

	variants = parsesam(vcf)
	
	outf.write('\n'.join(variants)+'\n')
	print '( %d / %d ) %s done ... %d denovos'%(k, total, fname, len(variants))
	progf.write('( %d / %d ) %s done ... %d denovos'%(k, total, fname, len(variants))+'\n')
        if len(variants) == 0:
		zeroctf.write(vcf+'\n')
outf.close()
progf.close()
