## Purpose: given a list of BAM files,
##      1. split into trio BAM lists (proband, father, mother)
##      2. write out shell scripts for each trio to qsub samtools calling
## Usage: python bamlist2triobamlist.py <LIST OF BAMS>
import sys

f = open(sys.argv[1],'r')

bamdir = '/share/data/CHD/Exome/PCGC_1115/fix_missing_1113/bamlists/'
scriptdir = '/share/shenlab/WES/PCGC_1115/CHD_NimbleGenV2/samtools/'

fam = {}
blist = {}
for line in f:
        tmp = line.strip().split('/')
        fullpath = line.strip()
        fname = tmp[-1]
        id = fname.split('.')[0].split('_')[0]
        fid = '-'.join(fname.split('.')[0].split('_')[0].split('-')[0:2])
        if len(id.split('-')) > 3:
            if id.split('-')[2] not in ['01', '02']:
                        id = '-'.join(id.split('-')[0:2])
            else:
                        id = '-'.join(id.split('-')[0:3])
        if id.startswith("G"):
            if not (id.endswith("A") or id.endswith("B")):
                fid = id
            else:
                fid = id[0:-1]
        if not fid in fam:
                fam[fid] = [id]
                blist[fid] = [fullpath]
        else:
            fam[fid].append(id)
            blist[fid].append(fullpath)
f.close()

## write out trio bam lists
for fam_id in blist.keys():
    bams = blist[fam_id]
    outf = open(bamdir+fam_id+'.list', 'w')
    outf.write('\n'.join(bams))
    outf.close()

## write out qsub scripts
for fam_id in fam.keys():
    outf=open(scriptdir+'tmp.'+fam_id+'.sh', 'w')
    out = []
    out.append('#!/bin/bash')
    out.append('SAM=/share/shenlab/WES/PCGC_1115/CHD_NimbleGenV2/src/ExmVC.7b.SamtoolsCall.sh')
    out.append('BAM=/share/shenlab/WES/PCGC_1115/CHD_NimbleGenV2/trio_bam_lists/'+fam_id+'.list')
    out.append('REF=/home/yufengshen/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.biocluster.sh')
    out.append('TGT=/share/shenlab/WES/PCGC_1115/CHD_NimbleGenV2/src/SeqCap_EZ_Exome_v2.hg19.targets.bed')
    out.append('LOG=tmp.'+fam_id+'.trio.log')
    out.append('OUTNAME='+fam_id+'.trio')
    out.append('qsub -N sam_'+fam_id+' $SAM -j 1 -i $BAM -r $REF -t $TGT -l $LOG -n $OUTNAME')
    outf.write('\n'.join(out))
outf.close()
