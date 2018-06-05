## Usage: python check_contamination.py <FILTERED DENOVOS FILE> > <FILTERED DENOVOS WITH COHORT_AF annotation>
## Purpose: For each variant, get the number of samples in this cohort that carry this variant and calculate cohort allele fraction (for downstream filtering)
## Dependencies: VCFs containing all variants (inh, denovo, etc) from a cohort
import sys
import gzip
import subprocess

f = open(sys.argv[1],'r') # filtered denovos  file

t2path = '/home/local/ARCS/hq2130/WES/CHD/CHD_NimbleGenV2/vcf0919/CHD_NimbleGenV2_VQSR.vcf.gz'
mepath = '/home/local/ARCS/hq2130/WES/CHD/CHD_MedExomeKit/vcf1118/CHD_MedExome_VQSR.vcf.gz'

cohortsize = 2645

## Function to get list of carriers for a given variant coordinate
## Dependencies: tabix ready cohort VCFs
def get_carriers(chr, pos, ref, alt, idx, path):
	carriers = []
	## tabix call 
	cmd = "tabix "+path+" "+chr+":"+pos+"-"+str(int(pos)+1)
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
	pout = p.communicate()[0].strip().split("\t")
	if len(pout) > 10 and len(pout)<6000:
		pout_ref, pout_alt = pout[3], pout[4]
		if ref == pout_ref and alt == pout_alt:
			colzip = dict(zip(pout, idx))
			for i in pout[9:]:
				ad = i.split(':')[1].split(',')
				refdp, altdp = ad[0], ad[1]
				if int(refdp) + int(altdp) > 0:
					vaf = float(altdp) / float(int(refdp)+int(altdp))
				else:
					vaf = 0.0
				if int(altdp) >= 2: # carriers = Nalt > 2
					carriers.append(colzip[i])
			return carriers
		else:
			return carriers
	else:
		return carriers

mef = gzip.open('/home/local/ARCS/hq2130/WES/CHD/CHD_MedExomeKit/vcf1118/CHD_MedExome_VQSR.vcf.gz', 'r')
meidx = {}
for line in mef:
	if line.startswith("#"):
		if line.startswith("#CHROM"):
			tmp = line.strip().split('\t')
			for i in xrange(len(tmp)):
				meidx[tmp[i]] = i
	else:
		break
mef.close()

t2f = gzip.open('/home/local/ARCS/hq2130/WES/CHD/CHD_NimbleGenV2/vcf0919/CHD_NimbleGenV2_VQSR.vcf.gz', 'r')
t2idx = {}
for line in t2f:
	if line.startswith("#"):
		if line.startswith("#CHROM"):
			tmp = line.strip().split('\t')
			for i in xrange(len(tmp)):
				t2idx[tmp[i]] = i
	else:
		break
t2f.close()


idx = {}
for line in f:
	tmp = line.strip().split('\t')
	if tmp[0] == "id":
		print '\t'.join(tmp) + '\t' + 'germ.carriers' + '\t' + 'cohort_AF'
		for i in xrange(len(tmp)):
			idx[tmp[i]] = i
	else:
		cflag = False
		id, chr, pos, ref, alt = tmp[idx['id']], tmp[idx['chr']], tmp[idx['pos']], tmp[idx['ref']], tmp[idx['alt']]
		me_carriers = get_carriers(chr, pos, ref, alt, meidx, mepath)
		t2_carriers = get_carriers(chr, pos, ref, alt, t2idx, t2path)
		# check if samples OTHER THAN ORIGINAL MUTATION CARRER carry as germline
		if id in me_carriers + t2_carriers:
			carriers = (me_carriers + t2_carriers).remove(id)
			if not carriers:
				carriers = []
		else:
			carriers = me_carriers + t2_carriers
			if not carriers:
				carriers = []
		if len(carriers) > 1:
			cflag = True
			print '\t'.join(tmp) + '\t' + str(len(carriers)) + '\t' + str(float(len(carriers))/cohortsize)
		else:
			print '\t'.join(tmp) + '\t' + '.' + '\t' + '0.0'
f.close()


