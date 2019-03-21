## Purpose: given a list of tissue BAM files,
##      1. split into trio BAM lists (proband, father, mother)
##      2. write out shell scripts for each trio to qsub samtools calling
## Usage: python tissue.bamlist_to_sh.py <LIST OF BAMS>
import sys

# /share/data/CHD/Exome/PCGC_tissue/scripts
# capture kit: /share/data/CHD/Exome/PCGC_tissue/src/xgen-exome-research-panel-targets.formatted.bed

bamd = {}
with open(sys.argv[2],'r') as f2: # bams1113.txt
	for line in f2:
		fullpath = line.strip()
		fname = fullpath.split('/')[-1]
		id = fname.split('.')[0]
		if id.endswith('-01') or id.endswith('-02'):
			fid = '-'.join(id.split('-')[0:2])
			if not fid in bamd:
				bamd[fid] = []
			bamd[fid].append(fullpath)

all_trio_blist_f = open('all.trio_bams.list', 'w') # master list of all bams in these trios
blistdir = '/share/data/CHD/Exome/PCGC_tissue/trio_bam_lists/'
scriptdir = '/share/data/CHD/Exome/PCGC_tissue/scripts/'
with open(sys.argv[1],'r') as f: # bams.list
	for line in f:
		fullpath = line.strip()
		fname = fullpath.split('/')[-1]
		id = fname.split('.')[0]
		fid = id.split('_')[0]

		if fid in bamd:
			
			## write out trio BAM lists
			bams = bamd[fid]
			bams.append(fullpath)
			outf = open(blistdir+fid+'.list', 'w')
			outf.write('\n'.join(bams))
			all_trio_blist_f.write('\n'.join(bams)+'\n')
			outf.close()

			## write out trio samtools calling scripts
			shf = open(scriptdir+id+'.sh', 'w')
			outsh = []
			outsh.append('#!/bin/bash')
			outsh.append('SAM=/share/data/CHD/Exome/PCGC_tissue/src/ExmVC.7b.SamtoolsCall.sh')
			outsh.append('BAM='+blistdir+fid+'.list')
			outsh.append('REF=/home/yufengshen/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.biocluster.sh')
			outsh.append('TGT=/share/data/CHD/Exome/PCGC_tissue/src/xgen-exome-research-panel-targets.formatted.bed')
			outsh.append('LOG=tmp.'+id+'.log')
			outsh.append('OUTNAME='+id)
			outsh.append('qsub -N sam_'+id+' $SAM -j 1 -i $BAM -r $REF -t $TGT -l $LOG -n $OUTNAME')
			shf.write('\n'.join(outsh))
			shf.close()
		
all_trio_blist_f.close()