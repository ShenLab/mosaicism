## Purpose: Given list of gene symbols, fetch for each gene
## 		1. full gene name
##		2. OMIM MIM #
##		3. any associated OMIM disease phenotypes
## Usage: python gene_to_pheno.py <genes list> <genemap2.txt> 
## NOTE: genemap2.txt file can be downloaded from https://www.omim.org/downloads/
import sys

f1 = open(sys.argv[1],'r') # text file containing gene symbols
omimf = open(sys.argv[2],'r') # omim genemap2.txt file

omimd = {}
idx = {}
for line in omimf:
	tmp = line.strip().split('\t')
	if line.startswith('#'):
		if 'Chromosome' in tmp[0]:
			idx = {col:index for index, col in enumerate(tmp)}
	else:
		mim = tmp[idx['Mim Number']]
		
		try:
			pheno = tmp[idx['Phenotypes']]
			if not len(pheno) > 1:
				pheno = '.'
		except:
			pheno = '.'

		genes = tmp[idx['Gene Symbols']]
		gname = tmp[idx['Gene Name']]
		for g in genes.split(','):
			omimd[g.strip()] = [mim, pheno, gname]
omimf.close()


print '\t'.join(['Gene Symbol', 'Gene Name', 'MIM Number', 'Phenotypes'])

for line in f1:
	original = line.strip()
	gene = original.split(' ')[0].strip()
	## gene symbols that contain single gene names
	if not ';' in gene:
		try:
			entry = omimd[gene]
			mim = entry[0]
			phen = entry[1]
			gname = entry[2]
		except:
			mim, phen, gname = '.', '.', '.'
		print '\t'.join([original, gname, mim, phen])
	## handle cases of ';'-separated gene symbols (ex: CBS;CBSL)
	else:
		mim, phen, gname = [], [], [], []
		for g in gene.split(';'):
			try:
				entry = omimd[g]
				mim.append(entry[0])
				phen.append(entry[1])
				gname.append(entry[2])
			except:
				mim.append('.')
				phen.append('.')
				gname.append('.')
		print original + '\t' + '\t'.join([';'.join(gname), ';'.join(mim), ';'.join(phen)])

f1.close()

