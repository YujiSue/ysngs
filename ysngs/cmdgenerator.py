import os
from ysngs import common

def getSRACmd(args):
    cmd = 'fasterq-dump ' + args['id'] + ' -O ' + args['outdir']
    if 'thread' in args:
        cmd += ' -e ' + str(args['thread'])
    return cmd   
def getFaiCmd(args):
    return 'samtools faidx ' + args['input']
def getBam2FqCmd(args):
    if args['type'] == 'single':
        return 'samtools fastq ' + args['input'] + ' > ' + args['output'][0]
    else:
        cmd = 'samtools collate -u -O ' + args['input'] + ' | samtools fastq '
        cmd += ' -1 ' + args['output'][0]
        cmd += ' -2 ' + args['output'][1]
        cmd += ' -0 /dev/null -s /dev/null -n'
        return cmd 

def getSam2BamCmd(args):
    cmd = 'samtools view'
    if 'thread' in args:
        cmd += ' -@ '+str(args['thread'])
    cmd += ' -b -o '+ args['output'] + ' ' + args['input']
    return cmd
def getBamSortCmd(args):
    cmd = 'samtools sort -l1 -T tmp'
    if 'thread' in args:
        cmd += ' -@ '+str(args['thread'])
    if 'ram' in args:
        cmd += ' -m ' + str(1000 if args['ram'] > 1 else int(args['ram']*1000)) + 'M'
    cmd += ' -O bam -o ' + args['output'] + ' ' + args['input']
    return cmd
def getBaiCmd(args):
    return 'samtools index ' + args['input']

# Pre-QC
def getQCCmd(args):
    return 'fastqc -o ' + args['outdir'] + ' ' + args['input'] + ' &'
def getCutCmd(args):
    cmd = 'cutadapt'
    if args['site'] == '5p':
        cmd += ' -g'
    elif args['site'] == '3p':
        cmd += ' -a'
    elif args['site'] == 'both':
        cmd += ' -b'
    cmd += ' ' + args['adapter'] + ' ' + args['input'] + ' > ' + args['output']
    return cmd
def getFPCmd(args):
    cmd = 'fastp -i ' + args['input'] + ' -o ' + args['output']
    if 'qual' in args:
        cmd += ' -q ' + str(args['qual'])
    if 'len' in args:
        cmd += ' -l ' + str(args['len'])
    if 'complex' in args:
        cmd += ' -y -Y ' + str(args['complex'])
    if 'thread' in args:
        cmd += ' -w ' + str(args['thread'])
    return cmd
# Mapping
def getBWAiCmd(args):
    return 'bwa index -p ' + args['output'] + ' ' + args['input']
def getBWACmd(args):
    cmd = 'bwa mem '
    if 'thread' in args:
        cmd += '-t '+str(args['thread'])
    if 'sr' in args and args['sr']:
        cmd += ' -M'
    if 'rgroup' in args:
        cmd += ' -R "@RG\\t' + args['rgroup'] + '"'
    cmd += ' ' + args['reference']
    if args['type'] == 'paired':
        cmd += ' ' + args['input'][0] + ' ' + args['input'][1]
    else:
        cmd += ' ' + args['input'][0]
    cmd += ' > ' + args['output']
    return cmd
def getBowtiCmd(args):
    cmd = 'bowtie2-build'
    if 'thread' in args:
        cmd += ' --threads ' + str(args['thread'])
    return cmd + ' -f ' + args['input'] + ' ' + args['output']
def getBowtCmd(args):
    cmd = 'bowtie2'
    if 'rgroup' in args:
        cmd += ' --rg-id ' + args['rgroup']['ID']
        for key in args['rgroup']:
            if key != 'ID':
                cmd += ' --rg ' + key + ':' + args['rgroup'][key]
    if args['type'] == 'single':
        cmd += ' -U ' + args['input'][0]
    elif args['type'] == 'paired':
        cmd += ' -1 ' + args['input'][0] + ' -2 ' + args['input'][1]
    if 'thread' in args:
        cmd += ' -p ' + str(args['thread'])
    cmd += ' -x ' + args['reference'] + ' -S '+ args['output']
    return cmd

def getSTARiCmd(args):
    cmd = 'STAR --runMode genomeGenerate'
    cmd += ' --genomeDir ' + args['reference']
    cmd += ' --genomeFastaFiles ' + args['refseq']
    if 'gtf' in args:
        cmd += ' --sjdbGTFfile ' + args['gtf']
    if 'thread' in args:
        cmd += ' --runThreadN ' + str(args['thread'])
    return cmd
def getSTARCmd(args):
    cmd = 'STAR --outSAMtype BAM SortedByCoordinate'
    cmd += ' --genomeDir ' + args['reference']
    cmd += ' --readFilesIn'
    for f in args['input']:
        cmd += ' ' + f 
    if 'thread' in args:
        cmd += ' --runThreadN ' + str(args['thread'])
    cmd += ' --outFileNamePrefix ' + args['output']
    return cmd

def getMarkDPCmd(args):
    cmd = 'java -jar picard.jar MarkDuplicates'
    cmd += ' -I ' + args['input']
    cmd += ' -O ' + args['output']
    cmd += ' -M ' + args['metric']
    return cmd
# Varaint Call
def getBCFVCCmd(args):
    cmd = 'bcftools mpileup -Ou'
    cmd += ' -f ' + args['reference']
    cmd += ' ' + args['input']
    cmd += ' | bcftools call -vm -Oz'
    cmd += ' -o ' + args['output']
    return cmd
def getGATKRefiCmd(args):
    cmd = 'gatk CreateSequenceDictionary' \
    + ' -R ' + args['input'] \
    + ' -O ' + args['output']
    return cmd
def getGATKVCCmd(args):
    cmd = 'gatk'
    if 'ram' in args:
        cmd += ' --java-options "-Xmx' + str(args['ram'])+'g"'
    cmd += ' HaplotypeCaller'
    if 'target' in args:
        cmd += ' -L ' + args['target']
    cmd += ' -R ' + args['reference']
    cmd += ' -I ' + args['input']
    cmd += ' -O ' + args['output']
    cmd += ' -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation'
    return cmd
def getGATKGenCmd(args):
    cmd = 'gatk'
    if 'ram' in args:
        cmd += ' --java-options "-Xmx' + str(args['ram'])+'g"'
    cmd += ' GenotypeGVCFs'
    cmd += ' -R ' + args['reference']
    cmd += ' -V ' + args['input']
    cmd += ' -O ' + args['output']
    return cmd

def getTVCCmd(args):
    cmd = 'python2 TVC/bin/variant_caller_pipeline.py' + \
    ' --input-bam ' + args['input'] + \
    ' --reference-fasta ' + args['reference'] + \
    ' --parameters-file ' + args['param'] + \
    ' --error-motifs ' + args['motif'] + \
    ' --generate-gvcf'
    if 'target' in args:
        cmd += ' --region-bed ' +  args['target']
    if 'hotspot' in args:
        cmd += ' --hotspot-vcf ' +  args['hotspot']
    if 'control' in args:
        cmd += ' --normal-bam ' +  args['control']
    if 'thread' in args:
        cmd += ' --num-threads ' + str(args['thread'])
    cmd += ' --output-dir ' + args['outdir']
    return cmd

def getIVCCmd(args):
    cmd = ''
    if args['mode'] == 'germ':
        cmd += 'configureStrelkaGermlineWorkflow.py'
        cmd += ' --bam ' + args['input']
        cmd += ' --referenceFasta ' + args['reference']
        cmd += ' --indelCandidates ' + args['output'] + '_indels'
        cmd += ' --runDir ' + args['outdir']
    elif args['mode'] == 'somat':
        cmd += 'configureStrelkaSomaticWorkflow.py'
        cmd += ' --normalBam ' + args['control']
        cmd += ' --tumorBam ' + args['input']
        cmd += ' --referenceFasta ' + args['reference']
        if 'variant' in args:
            cmd += ' --indelCandidates ' + args['variant']
        cmd += ' --runDir ' + args['outdir']
    return cmd
def getGDVCmd(args):
    idir, iname = os.path.split(args['input'])
    odir, oname = os.path.split(args['output'])
    rdir, rname = os.path.split(args['reference'])
    tdir = ''
    tname = ''
    if 'target' in args:
        tdir, tname = os.path.split(args['target'])
    cmd = 'sudo docker run'
    if 'gpu' in args:
        cmd += ' --gpus ' + str(args['gpu'])
    cmd += ' -v "' + rdir + '":"/REF_DIR"' 
    cmd += ' -v "' + idir + '":"/INPUT_DIR"' 
    cmd += ' -v "' + odir + '":"/OUTPUT_DIR"' 
    if tdir != '':
        cmd += ' -v "' + tdir + '":"/TARGET_DIR"' 
    cmd +=' google/deepvariant:"' + args['ver'] + '"' + \
        ' /opt/deepvariant/bin/run_deepvariant --model_type=WGS' 
    cmd += ' --ref ' + os.path.join('/REF_DIR', rname)
    cmd += ' --reads ' + os.path.join('/INPUT_DIR', iname)
    if 'target' in args and os.path.exists(args['target']):
        cmd += ' --regions ' + os.path.join('/TARGET_DIR', tname)
    if 'processor' in args and 0 < args['thread']:
        cmd += ' --num_shards ' + str(args['thread'])
    cmd += ' --output_gvcf ' + os.path.join('/OUTPUT_DIR', oname+'.g')
    cmd += ' --output_vcf ' + os.path.join('/OUTPUT_DIR', oname)
    return cmd

# CNVnator
def getMakeRootFileCmd(args):
    cmd = 'cnvnator -root ' + args['output'] + ' ' + args['input']
    return cmd
def getCNVnatHistCmd(args):
    cmd = 'cnvnator -root ' + args['input'] + ' -his ' + str(args['bin']) + ' -fasta ' + args['reference']
    return cmd
def getCNVNHistCmd(args):
    cmd = 'cnvnator -root ' + args['input'] + ' -his ' + str(args['bin'])
    return cmd
# ./cnvnator -root /media/ionDataArchive/sample.root -stat 10



# HLA-typing
def getHLAHDCmd(args):
    cmd = 'hlahd.sh'
    if 'thread' in args:
        cmd += ' -t ' + str(args['thread'])
    if 'min' in args:
        cmd += ' -m ' + str(args['min'])
    if 'trim' in args:
        cmd += ' -c ' + str(args['trim'])
    cmd += ' -f ' + args['freq']
    cmd += ' '.join(args['input'])
    cmd += ' ' + args['split']
    cmd += ' ' + args['dict']
    cmd += ' ' + args['sample']
    cmd += ' ' + args['outdir']

# Varaint detection
def getSutokuCmd(args):
    cmd = ''
    return cmd

# Variant effect prediction
def getEvepCmd(args):
    cmd = 'vep'
    cmd += ' -i ' + args['input']
    cmd += ' -o ' + args['output']
    if 'off' in args and args['off']:
        cmd += ' --offline --cache'
    if 'species' in args:
        cmd += ' --species ' + args['species']
    return cmd

def getPolyphen2Cmd(args):
    cmd = ''
    return cmd

def getSiftCmd(args):
    cmd = 'java -jar SIFT4G_Annotator.jar -c'
    cmd += ' -i ' + args['input'] # vcf
    cmd += ' -d ' + args['db'] # sift database
    cmd += ' -r ' + args['outdir']
    cmd += ' -t'
    return cmd

# Read count
## HTseq
def getHTSCmd(args):
    cmd = 'htseq-count -r pos -t exon -f bam'
    if 'qual' in args:
        cmd += ' -a ' + str(args['qual'])
    if 'thread' in args:
        cmd += ' -n ' + str(args['thread'])
    cmd += ' ' + args['input']
    cmd += ' ' + args['annotation']
    cmd += ' > ' + args['output']
    return cmd
## FeatureCounts
def getFCCmd(args):
    cmd = 'featureCounts'
    if args['type'] == 'paired':
        cmd += ' -p'
        if 'thread' in args:
            cmd += ' -T ' + str(args['thread'])
    elif args['type'] == 'long':
        cmd += ' -L'
    else:
        if 'thread' in args:
            cmd += ' -T ' + str(args['thread'])
    if 'strand' in args:
        cmd += ' -s '
    cmd += ' -t exon'
    cmd += ' -g gene_id'
    cmd += ' -a ' + args['annotation']
    cmd += ' -o ' + args['output']
    cmd += ' ' + ' '.join(args['input'])
    # featureCounts -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results.bam
    # featureCounts -p -t exon -g gene_id -a annotation.gtf -o counts.txt *.bam
    # featureCounts -T 8 -t exon -g gene_id -a annotation.gtf -o counts.txt input1.bam input2.bam input3.bam
    # -s 0 / 1 / 2 Strand specific
    #   
    return cmd
## Cufflinks
def getCuffLCmd(args):
    cmd = 'cufflinks --no-update-check'
    if 'platform' in args:
        if args['platform'] == 'ion':
            cmd += ' --library-type fr-secondstrand'
        else:
            cmd += ' --library-type fr-unstranded'
    if 'novel' in args:
        cmd += ' -g ' + args['annotation']
    else:
        cmd += ' -G ' + args['annotation']
    if 'mask' in args and os.path.exists(args['mask']):
        cmd += ' -M ' + args['mask']
    if 'thread' in args:
        cmd += ' -p ' + str(args['thread'])
    cmd += ' -o ' + args['output'] + ' ' + args['input']
    return cmd
## RSEM
def getRSEMCmd(args):
    cmd = ''

    return cmd

# DGE
def getCuffDCmd(args):
    cmd = 'cuffdiff --no-update-check'
    if 'platform' in args:
        if args['platform'] == 'ion':
            cmd += ' --library-type fr-secondstrand'
        else:
            cmd += ' --library-type fr-unstranded'
    if 'mask' in args and os.path.exists(args['mask']):
        cmd += ' -M ' + args['mask']
    if 'mincount' in args:
        cmd += ' -c ' + str(args['mincount'])
    if 'thread' in args:
        cmd += ' -p ' + str(args['thread'])
    if 'control' in args and args['control']:
        cmd += ' -g ' + args['control']
    cmd += ' -u -b ' + args['reference']['label'] + ' -o ' + opts['output']
    for label in args['labels']:
        cmd += label + ','
    cmd[-1] = ' '
    cmd += opts['input']
    for group in opts['groups']:
        cmd += ' '
        for reads in group:
            cmd += reads + ','
    cmd = cmd[:-1]
    return cmd

def getCuffMCmd(args):
    cmd = 'cuffmerge --no-update-check'
    if 'thread' in args:
        cmd += ' -p ' + str(args['thread'])
    cmd += ' -s ' + args['reference'] + ' ' + args['input']
    return cmd




def getCmd(app, cfg, opts):
    if app == 'sra':
        return getSRACmd(cfg, opts)
    elif app == 'bwai':
        return getBWAiCmd(cfg, opts)
    elif app == 'bwa':
        return getBWACmd(cfg, opts)
    elif app == 'tvc':
        return getTVCCmd(cfg, opts)
    elif app == 'ivc':
        return getIVCCmd(cfg, opts)
    elif app == 'gdv':
        return getGDVCmd(cfg, opts)


