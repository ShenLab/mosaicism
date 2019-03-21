## Usage: python check_cohortAF.py <FILTERED DENOVOS FILE> > <FILTERED DENOVOS WITH COHORT_AF annotation>
## Purpose: For each variant, get the number of samples in this cohort that carry this variant 
## 			and calculate cohort allele fraction (for downstream filtering)
## Dependencies: VCFs containing all variants (inh, denovo, etc) from a cohort

import sys
import gzip

# '/home/local/ARCS/hq2130/WES/CHD/CHD_MedExomeKit/vcf1118/CHD_MedExome_VQSR.vcf.gz'
# '/home/local/ARCS/hq2130/WES/CHD/CHD_NimbleGenV2/vcf0919/CHD_NimbleGenV2_VQSR.vcf.gz'

## function to parse joint vcf and return dictionary object {chr:pos:ref:alt  :  n_carriers}
def parse_carriers(vcf):
	carriers = {}
	with gzip.open(vcf, 'r') as f:
		for line in f:
			tmp = line.strip().split('\t')
			if not line.startswith('##'):
				if tmp[0] == '#CHROM':
					idx = {col:index for index, col in enumerate(tmp)}
				else:
					chr, pos, ref, alt = tmp[idx['#CHROM']], tmp[idx['POS']], tmp[idx['REF']], tmp[idx['ALT']]
					key = ':'.join([chr, pos, ref, alt])
					
					format = tmp[idx['FORMAT']]
					
					n = 0 ## counter for number of individuals carrying this variant with Nalt>=2

					## only count SNVs to speed up
					if (len(ref)==1 and len(alt)==1) and (not '.' in ref+alt):
						'''
						shortlist = [i for i, s in enumerate(tmp) if (('0/1' in s) or ('1/1' in s))]
						for i in shortlist:
							gt = tmp[i]
							tmpgt = dict(zip(tmp[idx['FORMAT']].split(':'), gt.split(':')))
							if int(tmpgt['AD'].split(',')[1]) >= 2:
								n += 1
						'''
						## iterate over genotypes
						for gt in tmp[idx['FORMAT']+1:]:
							tmpgt = dict(zip(tmp[idx['FORMAT']].split(':'), gt.split(':')))
							if int(tmpgt['AD'].split(',')[1]) >= 2:
								n += 1

					if n >= 1:
						carriers[key] = n
						print key + '\t' + str(n)
	return(carriers)

print "KEY" + '\t' + "n_carriers"
me_carriers = parse_carriers('/home/local/ARCS/hq2130/WES/CHD/CHD_MedExomeKit/vcf1118/CHD_MedExome_VQSR.vcf.gz')
t2_carriers = parse_carriers('/home/local/ARCS/hq2130/WES/CHD/CHD_NimbleGenV2/vcf0919/CHD_NimbleGenV2_VQSR.vcf.gz')


