## Usage: python annotate_variants.py <variants file> 
## Purpose: annotate variants with heart exp,  pli, etc for downstream analysis
import sys

f = open(sys.argv[1],'r') # variants file to annotate

en_f = open('/home/local/ARCS/alh2194/reference_files/ensembl_grch37_hg19.txt', 'r')
homsy_f = open('/home/local/ARCS/alh2194/reference_files/homsy_hhe_hbe.txt', 'r')
hi_f = open('/home/local/ARCS/alh2194/reference_files/gnomAD_constraint.txt', 'r')
em_ndd_f = open('/home/local/ARCS/alh2194/mosaic/jin_et_al/syndromic_analysis/jin_em_ndd.txt', 'r')
neuropheno_f = open('/home/local/ARCS/alh2194/reference_files/pcgc_phenotypes/jin_et_al.neuro_pheno.txt', 'r')
#age_f = open('/home/local/ARCS/alh2194/reference_files/pcgc_phenotypes/jin_et_al.ages.txt', 'r')
age_f = open('/home/local/ARCS/alh2194/reference_files/pcgc_phenotypes/jin_2530ids_age_information.converted.txt', 'r')
#union_f = open('/home/local/ARCS/alh2194/reference_files/CHD_risk_union_set.txt', 'r')
#union_f = open('/home/local/ARCS/alh2194/reference_files/CHD_riskgenes_Jin_S2.txt', 'r')
union_f = open('/home/local/ARCS/alh2194/reference_files/chdrisk_hhe_jin.txt', 'r')
#union_f = open('/home/local/ARCS/alh2194/reference_files/test_set.txt', 'r')

en = {}
for line in en_f:
	tmp = line.strip().split('\t')
	if not tmp[0].startswith("Ensembl"):
		enid = tmp[0]
		gene = tmp[2]
		en[enid] = gene
en_f.close()

homsy = {}
for line in homsy_f:
	tmp = line.strip().split('\t')
	if not tmp[0] == "Gene":
		gene, heart, brain = tmp[0], tmp[2], tmp[3]
		homsy[gene] = [heart, brain]
homsy_f.close()

hi = {}
for line in hi_f:
	tmp = line.strip().split('\t')
	if tmp[0] == 'gene':
		hifidx = {col:index for index, col in enumerate(tmp)}
	else:
		gene = tmp[hifidx['gene']]
		pli = tmp[hifidx['pLI']]
		if tmp[hifidx['canonical']] == 'true':
			hi[gene] = pli
hi_f.close()

ndd = {}
em = {}
for line in em_ndd_f:
	tmp = line.strip().split('\t')
	if not tmp[0] == "Blinded ID":
		em[tmp[0]] = tmp[1]
		ndd[tmp[0]] = tmp[2]
em_ndd_f.close()

phenod = {}
for line in neuropheno_f:
	tmp = line.strip().split('\t')
	if not tmp[0] == "id":
		id = tmp[0]
		sex = tmp[1]
		pheno = tmp[2]
		phenod[id] = [sex, pheno]
neuropheno_f.close()

aged = {}
for line in age_f:
	tmp = line.strip().split('\t')
	aged[tmp[0]] = tmp[-1]
age_f.close()

uniond = {}
for line in union_f:
	tmp = line.strip()
	uniond[tmp] = 'yes'
union_f.close()

for line in f:
	tmp = line.strip().split('\t')
	if tmp[0] == "id":
		head = '\t'.join(tmp) + '\t' + '\t'.join(['hexp',  'pli', 'age', 'sex', 'NDD_diagnosis', 'extracardiac', 'pheno', 'CHD_risk', 'CHD_subtype'])
		print head
		try:
			idx = tmp.index('gene')
		except:
			idx = tmp.index('Gene')
	else:
		## gene-level annotations
		gene = tmp[idx]
		if ',' in gene:
			gene = tmp[idx].split(',')[0]
		hexp, pli, = 'NA', 'NA'
		if gene in homsy:
			hexp = homsy[gene][0]
		if gene in hi:
			pli = hi[gene]
		#else:
		#	pli = '-1.0'
		try:
			chd_risk = uniond[gene]
		except:
			chd_risk = 'no'

		
		## sample-level annotations
		id = tmp[0]
		if '_' in tmp[0]:
			id = tmp[0].split('_')[0]
		ageout = '.'
		nddout = '.'
		emout = '.'
		phenoout = '.'
		if id in aged:
			ageout = aged[id]
		if id in ndd:
			nddout = ndd[id]
		if id in em:
			emout = em[id]
		if id in phenod:
			tmpphen = phenod[id]
			sex = tmpphen[0]
			phenoout = tmpphen[1]
		subtype = '.'
		if nddout == 'Yes' or emout == 'Yes':
			subtype = 'syndromic'
		elif nddout == 'No' and emout == 'No':
			subtype = 'isolated'
		elif (nddout=='Unknown' and emout=='No') or (nddout=='No' and emout=='Unknown') or (nddout=='Unknown' and emout=='Unknown'):
			subtype = 'unknown'
		print '\t'.join(tmp) + '\t' + '\t'.join([hexp, pli, ageout, sex, nddout, emout, phenoout, str(chd_risk), subtype])
f.close()
