## Usage: python filter_denovo.py <INPUT VCF> <INPUT VCF ANNOVAR MULTIANNO> > <FILTERED DENOVO FILE>
## Purpose: Filter for exonic SNVs with population frequency < 1e-4 and passing samtools PV4 filter
## Dependencies: Requires ANNOVAR-annotated VCF file 
import sys

## function to parse ANNOVAR multianno file
def parse_annovar(annof):
	out = {} # output dictionary for parsed ANNOVAR information {variant key : parsed information}
	idx = {} # index dict {colname : index}
	#key_cols = ['Otherinfo', 'Chr', 'Start', 'Ref', 'Alt', 
	#			'Gene.refGene', 'AAChange.refGene', 'Func.refGene', 'ExonicFunc.refGene',
	#			'Gene.ensGene', 'AAChange.ensGene', 
	#			'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS',
	#			'gnomAD_exome_ALL', 'gnomAD_exome_AFR', 'gnomAD_exome_AMR', 'gnomAD_exome_ASJ', 'gnomAD_exome_EAS',
	#			'gnomAD_exome_FIN', 'gnomAD_exome_NFE', 'gnomAD_exome_OTH', 'gnomAD_exome_SAS',
	#			'CADD_phred', 'MetaSVM_score', 'MetaSVM_pred', 'M-CAP_score', 'REVEL', 'MVP_rank']
	key_cols = ['Otherinfo', 'Chr', 'Start', 'Ref', 'Alt', 
				'Gene.refGene', 'AAChange.refGene', 'Func.refGene', 'ExonicFunc.refGene',
				'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS',
				'gnomAD_exome_ALL', 'gnomAD_exome_AFR', 'gnomAD_exome_AMR', 'gnomAD_exome_ASJ', 'gnomAD_exome_EAS',
				'gnomAD_exome_FIN', 'gnomAD_exome_NFE', 'gnomAD_exome_OTH', 'gnomAD_exome_SAS',
				'CADD_phred', 'MetaSVM_score', 'MetaSVM_pred', 'M-CAP_score', 'REVEL', 'MVP_rank']

	## deleteriousness thresholds for missense vairants
	thresh = {'CADD': 20.0, 'MCAP': 0.05, 'REVEL': 0.5,'MVP': 0.75}
	
	with open(annof, 'r') as f:
		for line in f:
			tmp = line.strip().split('\t')
			
			## Header line
			if tmp[0] in ['Chr', 'CHR', 'chr']:
				## sanity check: ensure that all columns used downstream are present
				head_overlap = list(set(key_cols).intersection(tmp))
				missing_cols = [c for c in key_cols if not c in head_overlap]
				if not len(head_overlap) == len(key_cols):
					print "## ERROR: ANNOVAR file missing key columns (%s)"%(','.join(missing_cols))
					sys.exit()
				
				## set up index dictionary to be able to get column indexes by name (idx[colname])
				idx = {col:index for index, col in enumerate(tmp)}
			else:
				## create variant key
				id = tmp[idx['Otherinfo']]
				chr, pos, ref, alt = tmp[idx['Chr']], tmp[idx['Start']], tmp[idx['Ref']], tmp[idx['Alt']]
				key = '|'.join([id, chr, pos, ref, alt])

				## get gene level information
				gene = tmp[idx['Gene.refGene']]
				aa_refgene = tmp[idx['AAChange.refGene']]
				#aa_ens = tmp[idx['AAChange.ensGene']]
				vfunc = tmp[idx['Func.refGene']]
				vclass = tmp[idx['ExonicFunc.refGene']]

				## get ExAC population MAF info
				exac_freq = tmp[idx['ExAC_ALL']]
				ex_afr = tmp[idx['ExAC_AFR']]
				ex_amr = tmp[idx['ExAC_AMR']]
				ex_eas = tmp[idx['ExAC_EAS']]
				ex_fin = tmp[idx['ExAC_FIN']]
				ex_nfe = tmp[idx['ExAC_NFE']]
				ex_oth = tmp[idx['ExAC_OTH']]
				ex_sas = tmp[idx['ExAC_SAS']]
				exac_popfreq = [exac_freq, ex_afr, ex_amr, ex_eas, ex_fin, ex_nfe, ex_oth, ex_sas]
				# get gnomAD population MAF info
				g_freq = tmp[idx['gnomAD_exome_ALL']]
				g_afr = tmp[idx['gnomAD_exome_AFR']]
				g_amr = tmp[idx['gnomAD_exome_AMR']]
				g_asj = tmp[idx['gnomAD_exome_ASJ']]
				g_eas = tmp[idx['gnomAD_exome_EAS']]
				g_fin = tmp[idx['gnomAD_exome_FIN']]
				g_nfe = tmp[idx['gnomAD_exome_NFE']]
				g_oth = tmp[idx['gnomAD_exome_OTH']]
				g_sas = tmp[idx['gnomAD_exome_SAS']]
				g_popfreq = [g_freq, g_afr, g_amr, g_asj, g_eas, g_fin, g_nfe, g_oth, g_sas]
				freqs = exac_popfreq + g_popfreq
				maxfreq = 0.0 ## initial maximum population MAF
				for fr in freqs:
					try: # for frequency values that readily convert to floats
						tmpfreq = float(fr)
						if tmpfreq >= maxfreq:
							maxfreq = tmpfreq
					except: # for frequency values like '.' or 'NA' or comma-separated frequency values
						tmpfreqs = fr.split(',')
						for _ in tmpfreqs:
							if _ in ['.', 'NA']:
								continue
							else:
								tmpfreq = float(_)
								if tmpfreq >= maxfreq:
									maxfreq = tmpfreq
				
				## handle pathogenicity scoring
				""" Range of values:
					.
					frameshift substitution
					nonframeshift substitution
					nonsynonymous SNV
					stopgain
					stoploss
					synonymous SNV
					unknown
				"""
				## check if there are unrecognized vclass values
				if not vclass in set(['.', 'frameshift substitution', 'nonframeshift substitution', 'nonsynonymous SNV', 'stopgain', 'stoploss', 'synonymous SNV', 'unknown']):
					print "## ERROR: ANNOVAR variant class unrecognized (%s)"%(vclass)
					sys.exit()
				
				## initialize values
				cadd, cadd_flag, meta, meta_flag, rev, revel_flag, mcap, mcap_flag, mvp, mvp_flag, adhoc_flag = 'NA','NA','NA','NA','NA','NA','NA','NA','NA', 'NA', 'NA'
				## only focus on exonic/splicing variants
				if ('exonic' in vfunc or 'splicing' in vfunc):
					cadd, meta, mcap, rev, mvp = tmp[idx['CADD_phred']], tmp[idx['MetaSVM_score']], tmp[idx['M-CAP_score']], tmp[idx['REVEL']], tmp[idx['MVP_rank']]
					
					## LGD variants - defined as vclass in ['stopgain', 'stoploss'] OR vfunc in ['splicing']
					## note: frameshift substitutions and nonframeshift substitutions look like indels -- ignore for now
					if vclass in ['stopgain', 'stoploss'] or 'splicing' in vfunc:
						cadd_flag, meta_flag, mcap_flag, revel_flag, mvp_flag, adhoc_flag = 'LGD', 'LGD', 'LGD', 'LGD', 'LGD', 'LGD'
					
					## Missense variants
					elif 'nonsynonymous' in vclass:
						cadd_flag, meta_flag, mcap_flag, revel_flag, mvp_flag, adhoc_flag = 'Bmis', 'Bmis', 'Bmis', 'Bmis', 'Bmis', 'Bmis'
						## CADD Dmis
						if ',' in tmp[idx['CADD_phred']]:
							c = tmp[idx['CADD_phred']].split(',')
							for score in c:
								if score == 'NA' or score == '.':
									cadd = '0.0'
								else:
									cadd = score
						if cadd == 'NA' or cadd == '.':
							cadd = '-1.0'
						if float(cadd) >= thresh['CADD']:
							cadd_flag = 'Dmis'
						## meta Dmis
						if 'D' in tmp[idx['MetaSVM_pred']]:
							meta_flag = 'Dmis' 
						## mcap
						if not mcap=='.':
							if float(mcap) >= thresh['MCAP']:
								mcap_flag = 'Dmis'
						## REVEL
						if not rev == '.':
							if float(rev) >= thresh['REVEL']:
								revel_flag = 'Dmis'
								
						## MVP
						if not mvp == '.':
							if float(mvp) >=thresh['MVP']:
								mvp_flag = 'Dmis'
						## adhoc (here we use REVEL)
						adhoc_flag = revel_flag

					## synonymous variants
					elif vclass in ['synonymousSNV','synonymous SNV']:
						cadd_flag, meta_flag, mcap_flag, revel_flag, mvp_flag, adhoc_flag = 'synonymous', 'synonymous', 'synonymous', 'synonymous', 'synonymous', 'synonymous'

					## non-frameshift deletion
					else:
						cadd_flag, meta_flag, mcap_flag, revel_flag, mvp_flag, adhoc_flag = 'other:' + vclass, 'other:' + vclass, 'other:' + vclass, 'other:' + vclass, 'other:' + vclass, 'other:' + vclass

				else:
					cadd_flag, meta_flag, mcap_flag, revel_flag, mvp_flag, adhoc_flag = 'intronic', 'intronic', 'intronic', 'intronic', 'intronic', 'intronic'

				key = '|'.join([id, chr, pos, ref, alt])
				#out[key] = {'gene':gene, 'aa_refgene':aa_refgene, 'aa_ens':aa_ens, 'vtype':vclass, 'vfunc':vfunc, 'maxfreq':str(maxfreq), 'exac_popfreq':';'.join(exac_popfreq), 'g_popfreq':';'.join(g_popfreq), 'cadd':cadd, 'cadd_flag':cadd_flag, 'meta':meta, 'meta_flag':meta_flag,'mcap':mcap, 'mcap_flag':mcap_flag, 'revel':rev, 'revel_flag':revel_flag,'mvp':mvp, 'mvp_flag':mvp_flag, 'vclass':adhoc_flag}
				out[key] = {'gene':gene, 'aa_refgene':aa_refgene, 'vtype':vclass, 'vfunc':vfunc, 'maxfreq':str(maxfreq), 'exac_popfreq':';'.join(exac_popfreq), 'g_popfreq':';'.join(g_popfreq), 'cadd':cadd, 'cadd_flag':cadd_flag, 'meta':meta, 'meta_flag':meta_flag,'mcap':mcap, 'mcap_flag':mcap_flag, 'revel':rev, 'revel_flag':revel_flag,'mvp':mvp, 'mvp_flag':mvp_flag, 'vclass':adhoc_flag}
	
	return(out)

## Function to apply hard filters to each variant site based on Yufeng's suggestion
# Exclusion criteria: Nalt<6, maxpopfreq>1e-3, PV4: p<0.05 for any of the 4 tests
def check_filter_sam(vline, idx):
	flag = True
	tmp = vline.strip().split('\t')

	#chr, pos, ref, alt = tmp[0], tmp[1], tmp[3], tmp[4]
	#info = tmp[7]
	chr, pos, ref, alt = tmp[idx['#CHROM']], tmp[idx['POS']], tmp[idx['REF']], tmp[idx['ALT']]
	info = tmp[idx['INFO']]

	pv4 = '.'
	# p-value cutoffs for samtools PV4 tests
	pv4_cutoff = 1e-3 
	pv4_cutoff_mapQ = 1e-6
	
	for i in info.split(';'):
		if i.startswith("PV4="):
			pv4 = i
			tests = i.split('=')[1].split(',')
			## only apply pv4_cutoff to baseQ bias, taildistbias, strand bias -- skip mapQ
			for t in xrange(len(tests)):
				## for PV4 values other than mapQ
				if not t == 2:
					if float(tests[t]) < pv4_cutoff:
						flag = False
				## for PV4 mapQ
				else:
					if float(tests[t]) < pv4_cutoff_mapQ:
						flag = False
	return(flag, pv4)


#################
##
## MAIN FUNCTION
##
#################

## Usage/argument checking
try:
	arg1, arg2 = sys.argv[1], sys.argv[2]
except:
	print '\n## ERROR: MISSING ARGUMENTS'
	print "## USAGE:\npython filter_denovo.py <INPUT VCF> <INPUT VCF ANNOVAR MULTIANNO> > <FILTERED DENOVO FILE>\n"
	sys.exit()

## Read ANNOVAR multianno file and create dictionary
annof = sys.argv[2]
annod = parse_annovar(annof)

## filter variants
fullannof = open('fullanno.denovo_nofilt.txt', 'w') # Full, unfiltered, annotated variants file with full annotations for PV4 and which filter the site failed (for QC)
filterd = {'not_snv':[], 'vfunc':[], 'popfreq':[], 'samflag':[], 'muc_hla':[]} # Dictionary tracking how many variants fail each filtering criteria
with open(sys.argv[1],'r') as f1:
	## variant caller
	src = 'samtools'
	
	## population frequency cutoffs
	popfreq_cutoff = 1e-4

	for line in f1:
		tmp = line.strip().split('\t')
		if tmp[0] == '#CHROM':
			idx = {col:index for index, col in enumerate(tmp)}
			## header for filtered file
			head = ['id', 'gene', 'chr', 'pos', 'ref', 'alt', 'refdp', 'altdp', 'aachange', 'maxfreq', 'vfunc', 'vtype', 'cadd', 'cadd_flag', 'meta', 'meta_flag', 'mcap', 'mcap_flag', 'revel', 'revel_flag', 'mvp', 'mvp_flag', 'vclass', 'src', 'adfref', 'adfalt', 'adrref', 'adralt']
			print '\t'.join(head)
			## header for fullanno file
			pv4head = ['id', 'gene', 'chr', 'pos', 'ref', 'alt', 'refdp', 'altdp', 'aachange', 'maxfreq', 'vfunc', 'vtype', 'cadd', 'cadd_flag', 'meta', 'meta_flag', 'mcap', 'mcap_flag', 'revel', 'revel_flag', 'mvp', 'mvp_flag', 'vclass', 'src', 'adfref', 'adfalt', 'adrref', 'adralt', 'strand_bias', 'baseQ_bias','mapQ_bias','tail_dist_bias', 'pre_filter']
			fullannof.write('\t'.join(pv4head)+'\n')
		else:
			gt_format = tmp[idx['FORMAT']].split(':')
			id = tmp[idx['ID']]
			chr, pos, ref, alt = tmp[idx['#CHROM']], tmp[idx['POS']], tmp[idx['REF']], tmp[idx['ALT']]
			
			## handle multiallelic
			if len(alt.split(',')) > 1: 
				pb_ad = tmp[idx['PROBAND']].split(':')[gt_format.index('AD')]
				pb_altidx = map(int, pb_ad.split(',')[1:]).index(max(map(int, pb_ad.split(',')[1:]))) + 1
				refdp = pb_ad.split(',')[0]
				altdp = pb_ad.split(',')[pb_altidx]
				pb_adf = tmp[idx['PROBAND']].split(':')[gt_format.index('ADF')].split(',')
				pb_adr = tmp[idx['PROBAND']].split(':')[gt_format.index('ADR')].split(',')
				dp4 = pb_adf + pb_adr
				adfref, adfalt = pb_adf[0], pb_adf[1]
				adrref, adralt = pb_adr[0], pb_adr[1]
			else:
				pb_ad = tmp[idx['PROBAND']].split(':')[gt_format.index('AD')]
				refdp = pb_ad.split(',')[0]
				altdp = pb_ad.split(',')[1]
				pb_adf = tmp[idx['PROBAND']].split(':')[gt_format.index('ADF')].split(',')
				pb_adr = tmp[idx['PROBAND']].split(':')[gt_format.index('ADR')].split(',')
				dp4 = pb_adf + pb_adr
				adfref, adfalt = pb_adf[0], pb_adf[1]
				adrref, adralt = pb_adr[0], pb_adr[1]

			key = '|'.join([id, chr, pos, ref, alt])
			
			## fetch corresponding ANNOVAR annotation for variant
			try:
				anno = annod[key]
				gene = anno['gene']
				aachange = anno['aa_refgene']
				maxfreq = anno['maxfreq']
				vfunc = anno['vfunc']
				vtype = anno['vtype']
				cadd = anno['cadd']
				cadd_flag = anno['cadd_flag']
				meta = anno['meta']
				meta_flag = anno['meta_flag']
				mcap = anno['mcap']
				mcap_flag = anno['mcap_flag']
				revel = anno['revel']
				revel_flag = anno['revel_flag']
				mvp = anno['mvp']
				mvp_flag = anno['mvp_flag']
				vclass = anno['vclass']	
			except:
				print '## ANNOVAR ERROR'
				break
				aachange, gene, maxfreq, vfunc, vtype, cadd, cadd_flag, meta, meta_flag, mcap, mcap_flag, revel, revel_flag, mvp, mvp_flag, vclass = '.','.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'
			
			## check samtools PV4
			checksam = check_filter_sam(line, idx)
			samflag = checksam[0]
			try:
				pv4 = checksam[1].split('=')[1].split(',')
			except:
				pv4 = ['.', '.', '.', '.']

			## flag for SNV vs. indel
			snv = 1
			if len(ref)>1 or len(alt)>1 or '*' in ref+alt:
				snv = 0

			## variant output
			out = map(str, [id, gene, chr, pos, ref, alt, refdp, altdp, aachange, maxfreq, vfunc, vtype, cadd, cadd_flag, meta, meta_flag, mcap, mcap_flag, revel, revel_flag, mvp, mvp_flag, vclass, src, adfref, adfalt, adrref, adralt])


			## STANDARD FILTERS
			if snv == 1 and ('exonic' in vfunc or 'splicing' in vfunc ) and not ('ncRNA' in vfunc):
				if samflag:
					if float(maxfreq) < popfreq_cutoff:
						if not ('MUC' in gene or 'HLA' in gene):
							print '\t'.join(out)
		
			## annotate sites not passing filters 
			flag = []
			if not snv == 1:
				flag.append('not_snv')
				filterd['not_snv'].append(out)
			if not ('exonic' in vfunc or 'splicing' in vfunc ) or ('ncRNA' in vfunc):
				flag.append('vfunc')
				filterd['vfunc'].append(out)
			if not samflag:
				flag.append('samflag')
				filterd['samflag'].append(out)
			if not maxfreq == '.':
				if not float(maxfreq) < popfreq_cutoff:
					flag.append('popfreq_gt_1e-4')
					filterd['popfreq'].append(out)
			if ('MUC' in gene or 'HLA' in gene):
				flag.append('muc_hla')
				filterd['muc_hla'].append(out)
			## write out filtered sites to fullanno file
			anno_out = '\t'.join(out) + '\t' +'\t'.join(pv4) +'\t' + ','.join(flag)
			fullannof.write(anno_out + '\n')		

	f1.close()
	fullannof.close()

filterlog = open('filter_log.txt','w')
for k in filterd:
	filterlog.write(k + '\t' + str(len(filterd[k]))+'\n')
filterlog.close()
