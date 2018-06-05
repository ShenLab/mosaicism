## Usage: python filter_denovo.py <INPUT VCF> <INPUT VCF ANNOVAR MULTIANNO> > <FILTERED DENOVO FILE>
## Purpose: Filter for exonic SNVs with population frequency < 1e-4 and passing samtools PV4 filter
## Dependencies: Requires ANNOVAR-annotated VCF file 
import sys
import subprocess
import signal
import gzip

f = open(sys.argv[1],'r')
annof = sys.argv[2]

## Function: parse ANNOVAR output for annotations and pathogenicity predictions
def parse_annovar(annof):
	out = {}
	idx = {}
	file = open(annof, 'r')
	for line in file:
		tmp = line.strip().split('\t')
		if tmp[0] == "Chr":
			for i in xrange(len(tmp)):
				idx[tmp[i]] = i
		else:
			id = tmp[idx['Otherinfo']]
			chr = tmp[idx['Chr']]
			pos = tmp[idx['Start']]
			ref = tmp[idx['Ref']]
			alt = tmp[idx['Alt']]
			gene = tmp[idx['Gene.refGene']]
			aachangelist = []
			tmpaachange = tmp[idx['AAChange.refGene']]
			for a in tmpaachange.split(','):
				for b in a.split(':'):
					if b.startswith('p.'):
						aachangelist.append(b)
			if len(aachangelist) >= 1:
				aachange = ','.join(aachangelist)
			else:
				aachange = '.'

			vfunc = tmp[idx['Func.refGene']]
			vclass = tmp[idx['ExonicFunc.refGene']]
			## handle population information
			exac_freq = tmp[idx['ExAC_ALL']]
			ex_afr = tmp[idx['ExAC_AFR']]
			ex_amr = tmp[idx['ExAC_AMR']]
			ex_eas = tmp[idx['ExAC_EAS']]
			ex_fin = tmp[idx['ExAC_FIN']]
			ex_nfe = tmp[idx['ExAC_NFE']]
			ex_oth = tmp[idx['ExAC_OTH']]
			ex_sas = tmp[idx['ExAC_SAS']]
			exac_popfreq = ';'.join([exac_freq, ex_afr, ex_amr, ex_eas, ex_fin, ex_nfe, ex_oth, ex_sas])
			# get gnomad info
			g_freq = tmp[idx['gnomAD_exome_ALL']]
			g_afr = tmp[idx['gnomAD_exome_AFR']]
			g_amr = tmp[idx['gnomAD_exome_AMR']]
			g_asj = tmp[idx['gnomAD_exome_ASJ']]
			g_eas = tmp[idx['gnomAD_exome_EAS']]
			g_fin = tmp[idx['gnomAD_exome_FIN']]
			g_nfe = tmp[idx['gnomAD_exome_NFE']]
			g_oth = tmp[idx['gnomAD_exome_OTH']]
			g_sas = tmp[idx['gnomAD_exome_SAS']]
			g_popfreq = ';'.join([g_freq, g_afr, g_amr, g_asj, g_eas, g_fin, g_nfe, g_oth, g_sas])
			allfreqs = [g_freq, g_afr, g_amr, g_asj, g_eas, g_fin, g_nfe, g_oth, g_sas]
			maxfreq = 0.0
			for fr in allfreqs:
					try:
						freq = float(fr)
						if freq >= maxfreq:
							maxfreq = freq
					except:
						if not ',' in fr:
							if ('.' in fr or 'NA' in fr): # handle . and NA cases
								continue
						elif ',' in fr: # handle x,x cases
							tmpfr = fr.split(',')
							for tmpf in tmpfr:
								if not tmpf in ['.', 'NA']:
									freq = float(tmpf)
									if freq >= maxfreq:
										maxfreq = freq
					continue
			
			## handle pathogenicity scoring
			cadd, cadd_flag, meta, meta_flag, rev, revel_flag, mcap, mcap_flag, mvp, mvp_flag, adhoc_flag = 'NA','NA','NA','NA','NA','NA','NA','NA','NA', 'NA', 'NA'
			if ("exonic" in vfunc or "splicing" in vfunc):
				cadd = tmp[idx['CADD_phred']]
				meta = tmp[idx['MetaSVM_score']]
				mcap = tmp[idx['M-CAP_score']]
				rev = tmp[idx['REVEL']]
				mvp = tmp[idx['MVP_rank']]
				if "stopgain" in vclass or "stoploss" in vclass or "frameshiftdeletion" in vclass or "frameshiftinsertion" in vclass or "splicing" in vfunc:
					if not 'non' in vclass:
						meta_flag = 'LGD'
						cadd_flag = 'LGD'
						mcap_flag = 'LGD'
						revel_flag = 'LGD'
						mvp_flag = 'LGD'
						adhoc_flag = 'LGD'
				elif "nonsynonymous" in vclass:
					meta_flag = 'Bmis'
					cadd_flag = 'Bmis'
					mcap_flag = 'Bmis'
					revel_flag = 'Bmis'
					mvp_flag = 'Bmis'
					adhoc_flag = 'Bmis'
					## meta Dmis
					if 'D' in tmp[idx['MetaSVM_pred']]:
						meta_flag = 'Dmis'
					## CADD Dmis
					if ',' in tmp[idx['CADD_phred']]:
						c = tmp[idx['CADD_phred']].split(',')
						for score in c:
							if score == "NA" or score == ".":
								cadd = "0.0"
							else:
								cadd = score
					if cadd == "NA" or cadd == ".":
						cadd = "-1.0"
					if float(cadd) >= 20.0:
						cadd_flag = 'Dmis'
					## mcap
					if not mcap==".":
						if float(mcap) >= 0.05:
							mcap_flag = 'Dmis'
					## REVEL
					if not rev == ".":
						if float(rev) >= 0.5:
							revel_flag = 'Dmis'
							adhoc_flag = 'Dmis'
					## MVP
					if not mvp == ".":
						if float(mvp) >= 0.75:
							mvp_flag = 'Dmis'

				elif vclass in ['synonymousSNV','synonymous SNV']:
					meta_flag = 'synonymous'
					cadd_flag = 'synonymous'
					mcap_flag = 'synonymous'
					revel_flag = 'synonymous'
					mvp_flag = 'synonymous'
					adhoc_flag = 'synonymous'
				elif vclass in ['nonframeshiftdeletion', 'nonframeshiftinsertion', 'non-frameshiftdeletion', 'non-frameshiftinsertion']:
					meta_flag = 'non-damaging'
					cadd_flag = 'non-damaging'
					mcap_flag = 'non-damaging'
					revel_flag = 'non-damaging'
					mvp_flag = 'non-damaging'
					adhoc_flag = 'non-damaging'
				else:
					meta_flag = 'other: '+vclass
					cadd_flag = 'other: '+vclass
					mcap_flag = 'other: '+vclass
					revel_flag = 'other: '+vclass
					mvp_flag = 'other: '+vclass
					adhoc_flag = 'other: '+vclass

			else:
				meta_flag = 'intronic'
				cadd_flag = 'intronic'
				mcap_flag = 'intronic'
				revel_flag = 'intronic'
				mvp_flag = 'intronic'
				adhoc_flag = 'intronic'
			key = '|'.join([id, chr, pos, ref, alt])
			out[key] = {'gene':gene, 'aachange':aachange, 'vtype':vclass, 'vfunc':vfunc, 'maxfreq':str(maxfreq), 'exac_popfreq':exac_popfreq, 'g_popfreq':g_popfreq, 'cadd':cadd, 'cadd_flag':cadd_flag, 'meta':meta, 'meta_flag':meta_flag,'mcap':mcap, 'mcap_flag':mcap_flag, 'revel':rev, 'revel_flag':revel_flag,'mvp':mvp, 'mvp_flag':mvp_flag, 'vclass':adhoc_flag}
	file.close()
	return(out)

## Function to apply hard filters to each variant site based on Yufeng's suggestion
# Exclusion criteria: Nalt<6, maxpopfreq>1e-3, PV4: p<0.05 for any of the 4 tests
def check_filter_sam(vline):
	flag = True
	tmp = vline.strip().split('\t')
	chr, pos, ref, alt = tmp[0], tmp[1], tmp[3], tmp[4]
	info = tmp[7]
	pv4 = '.'
	pv4_cutoff = 1e-3
	## handle multiallelic sites --> get correct alt index
	if len(alt.split(',')) > 1:
		pb_ad = tmp[9].split(':')[-1]
		pb_altidx = map(int, pb_ad.split(',')[1:]).index(max(map(int, pb_ad.split(',')[1:])))
		for i in info.split(';'):
			if i.startswith("PV4="):
				pv4 = i
				tests = i.split('=')[1].split(',')
				## only apply pv4_cutoff to baseQ bias, taildistbias, strand bias -- skip mapQ
				for t in xrange(len(tests)):
					if not t == 2:
						if float(tests[t]) < pv4_cutoff:
							flag = False
		pbgeno = tmp[9]
	## biallelic sites
	else:
		for i in info.split(';'):
			if i.startswith("PV4="):
				pv4 = i
				tests = i.split('=')[1].split(',')
				## only apply pv4_cutoff to baseQ bias, taildistbias, strand bias -- skip mapQ
				for t in xrange(len(tests)):
					if not t == 2:
						if float(tests[t]) < pv4_cutoff:
							flag = False
		pbgeno = tmp[9]
	return(flag, pv4)

## MAIN FUNCTION
annod = parse_annovar(annof) # Parse annovar file

fullannof = open('fullanno.denovo_nofilt.txt', 'w') # Full, unfiltered, annotated variants file with full annotations for PV4 and which filter the site failed (for QC)
filterd = {'not_snv':[], 'vfunc':[], 'popfreq':[], 'samflag':[], 'muc_hla':[]} # Dictionary tracking how many variants fail each filtering criteria
done = []
for line in f:
	if line.startswith("#CHROM"):
		head = ['id', 'gene', 'chr', 'pos', 'ref', 'alt', 'refdp', 'altdp', 'aachange', 'maxfreq', 'vfunc', 'vtype', 'cadd', 'cadd_flag', 'meta', 'meta_flag', 'mcap', 'mcap_flag', 'revel', 'revel_flag', 'mvp', 'mvp_flag', 'vclass', 'src', 'adfref', 'adfalt', 'adrref', 'adralt']
		print '\t'.join(head)
		pv4head = ['id', 'gene', 'chr', 'pos', 'ref', 'alt', 'refdp', 'altdp', 'aachange', 'maxfreq', 'vfunc', 'vtype', 'cadd', 'cadd_flag', 'meta', 'meta_flag', 'mcap', 'mcap_flag', 'revel', 'revel_flag', 'mvp', 'mvp_flag', 'vclass', 'src', 'adfref', 'adfalt', 'adrref', 'adralt', 'strand_bias', 'baseQ_bias','mapQ_bias','tail_dist_bias', 'pre_filter']
		fullannof.write('\t'.join(pv4head)+'\n')
	else:
		tmp = line.strip().split('\t')
		gt_format = tmp[8].split(':')
		if len(tmp) == 12:
			id = tmp[2]
			chr, pos, ref, alt = tmp[0], tmp[1], tmp[3], tmp[4]
			if len(alt.split(',')) > 1: ## handle multiallelic
				pb_ad = tmp[9].split(':')[gt_format.index('AD')]
				pb_altidx = map(int, pb_ad.split(',')[1:]).index(max(map(int, pb_ad.split(',')[1:]))) + 1
				refdp = pb_ad.split(',')[0]
				altdp = pb_ad.split(',')[pb_altidx]
				pb_adf = tmp[9].split(':')[gt_format.index('ADF')].split(',')
				pb_adr = tmp[9].split(':')[gt_format.index('ADR')].split(',')
				dp4 = pb_adf + pb_adr
				adfref, adfalt = pb_adf[0], pb_adf[1]
				adrref, adralt = pb_adr[0], pb_adr[1]
			else:
				pb_ad = tmp[9].split(':')[gt_format.index('AD')]
				refdp = pb_ad.split(',')[0]
				altdp = pb_ad.split(',')[1]
				pb_adf = tmp[9].split(':')[gt_format.index('ADF')].split(',')
				pb_adr = tmp[9].split(':')[gt_format.index('ADR')].split(',')
				dp4 = pb_adf + pb_adr
				adfref, adfalt = pb_adf[0], pb_adf[1]
				adrref, adralt = pb_adr[0], pb_adr[1]
		else: # samples with many family members
			id = tmp[2]
			chr, pos, ref, alt = tmp[0], tmp[1], tmp[3], tmp[4]
			## FIGURE OUT WHICH GT BELONGS TO PROBAND
			for gt in tmp[9:]: 
				ad = gt.split(':')[gt_format.index('AD')]
				if ad.split(',')[1] != '0':
					pb_ad = tmp[9].split(':')[gt_format.index('AD')]
					refdp = pb_ad.split(',')[0]
					altdp = pb_ad.split(',')[1]
					pb_adf = tmp[9].split(':')[gt_format.index('ADF')].split(',')
					pb_adr = tmp[9].split(':')[gt_format.index('ADR')].split(',')
					dp4 = pb_adf + pb_adr
					adfref, adfalt = pb_adf[0], pb_adf[1]
					adrref, adralt = pb_adr[0], pb_adr[1]
		
		key = '|'.join([id, chr, pos, ref, alt])
		try:
			anno = annod[key]
			gene = anno['gene']
			aachange = anno['aachange']
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
			aachange, gene, maxfreq, vfunc, vtype, cadd, cadd_flag, meta, meta_flag, mcap, mcap_flag, revel, revel_flag, mvp, mvp_flag, vclass = '.','.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'
		
		src = 'samtools'
		# handle samtools parsing and filtering using PV4
		checksam = check_filter_sam(line)
		samflag = checksam[0]
		try:
			pv4 = checksam[1].split('=')[1].split(',')
		except:
			pv4 = ['.', '.', '.', '.']
		
		# flag for SNV vs. indel
		snv = 1
		if len(ref)>1 or len(alt)>1 or '*' in ref+alt:
			snv = 0
		
		out = map(str, [id, gene, chr, pos, ref, alt, refdp, altdp, aachange, maxfreq, vfunc, vtype, cadd, cadd_flag, meta, meta_flag, mcap, mcap_flag, revel, revel_flag, mvp, mvp_flag, vclass, src, adfref, adfalt, adrref, adralt])

		## STANDARD FILTERS
		if snv == 1 and ('exonic' in vfunc or 'splicing' in vfunc ) and not ('ncRNA' in vfunc):
			if samflag:
				if float(maxfreq) < 1e-4:
					if not ('MUC' in gene or 'HLA' in gene):
						out = map(str, [id, gene, chr, pos, ref, alt, refdp, altdp, aachange, maxfreq, vfunc, vtype, cadd, cadd_flag, meta, meta_flag, mcap, mcap_flag, revel, revel_flag, mvp, mvp_flag, vclass, src, adfref, adfalt, adrref, adralt])
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
		if not float(maxfreq) < 1e-4:
			flag.append('popfreq_gt_1e-4')
			filterd['popfreq'].append(out)
		if ('MUC' in gene or 'HLA' in gene):
			flag.append('muc_hla')
			filterd['muc_hla'].append(out)
		anno_out = '\t'.join(out) + '\t' +'\t'.join(pv4) +'\t' + ','.join(flag)
		fullannof.write(anno_out + '\n')		

f.close()
fullannof.close()

filterlog = open('filter_log.txt','w')
for k in filterd:
	filterlog.write(k + '\t' + str(len(filterd[k]))+'\n')
filterlog.close()
