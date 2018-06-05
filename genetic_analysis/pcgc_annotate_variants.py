## Usage: python annotate_variants.py <variants file> 
## Purpose: annotate variants with heart exp,  pli, etc for downstream analysis
import sys

f = open(sys.argv[1],'r') # variants file to annotate

en_f = open('/home/local/ARCS/alh2194/reference_files/ensembl_grch37_hg19.txt', 'r')
homsy_f = open('/home/local/ARCS/alh2194/reference_files/homsy_hhe_hbe.txt', 'r')
hi_f = open('/home/local/ARCS/alh2194/reference_files/exac_pli_combined.max.txt', 'r')
em_ndd_f = open('/home/local/ARCS/alh2194/mosaic/jin_et_al/syndromic_analysis/jin_em_ndd.txt', 'r')
neuropheno_f = open('/home/local/ARCS/alh2194/reference_files/pcgc_phenotypes/jin_et_al.neuro_pheno.txt', 'r')
age_f = open('/home/local/ARCS/alh2194/reference_files/pcgc_phenotypes/jin_et_al.ages.txt', 'r')

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
	if not tmp[0] == "gene":
		gene = tmp[0]
		pli = tmp[-1]
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
	aged[tmp[0]] = tmp[1]
age_f.close()


for line in f:
	tmp = line.strip().split('\t')
	if tmp[0] == "id":
		head = '\t'.join(tmp) + '\t' + '\t'.join(['hexp',  'pli', 'age', 'sex', 'NDD_diagnosis', 'extracardiac', 'pheno' ])
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
		else:
			pli = '-1.0'
		
		## sample-level annotations
		id = tmp[0]
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
		print '\t'.join(tmp) + '\t' + '\t'.join([hexp, pli, ageout, sex, nddout, emout, phenoout])
f.close()
