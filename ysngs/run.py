import os
import subprocess
from subprocess import PIPE 
import common

class apprun:
  def __init__(self, config):
    self.cfg = config
  
  def printMsgLine(self, msg = ''):
    lines = msg.splitlines()
    for l in lines:
      print(' >', l)

  def execCmd(self, cmd, verbose = False):
    print('Run: >', cmd)
    proc = subprocess.run(cmd, stdout=PIPE, stderr=PIPE, text=True, shell=True)
    if proc.returncode:
      print('Error:')
      self.printMsgLine(proc.stderr)
      print('')
      return False
    else:
      if verbose:
        self.printMsgLine(proc.stdout)
      else:
        self.printMsgLine('Completed.')
        print('')
      return True
    
  def downloadSRA(self, srid, output = '.', option = {'thread':8}):
    cmd = 'fasterq-dump ' + srid + ' -O ' + output
    if option['thread']:
      cmd += ' -e ' + str(option['thread'])
    return self.execCmd(cmd)

  def runCutter(self, adaptor = '', site = '', input = '', output = ''):
    cmd = 'cutadapt'
    if site == '5p':
      cmd += ' -g'
    elif site == '3p':
      cmd += ' -a'
    elif site == 'both':
      cmd += ' -b'
    cmd += ' ' + adaptor + ' ' + input + ' > ' + output
    return self.execCmd(cmd)

  def runFastQFilter(self, input = '', output = '', param = {}):
    common.addPath(self.cfg.APPS_DIR)
    cmd = 'fastp -i ' + input + ' -o ' + output
    if 'min_qual' in param:
      cmd += ' -q ' + str(param['min_qual'])
    if 'min_len' in param:
      cmd += ' -l ' + str(param['min_len'])
    if 'min_complex':
      cmd += ' -y -Y ' + str(param['min_complex'])
    if 'thread' in param:
      cmd += ' -w ' + str(param['thread'])
    return self.execCmd(cmd)
  
  def makeRGText(self, rginfo):
    rg = 'ID:' + rginfo['ID']
    for key in rginfo:
      if key != 'ID':
        rg += '\\t' + key + ':' + rginfo[key]
    return rg

  def hasFai(self, ref):
    return os.path.exists(ref+'.fai')

  def makeFai(self, ref):
    cmd = 'samtools faidx ' + ref
    return self.execCmd(cmd)

  def hasBWARefIndex(self, name):
    return os.path.exists(os.path.join(self.cfg.REFERENCE_DIR, name+'.pac'))

  def makeBWARefIndex(self, refpath, name):
    cmd = 'bwa index -p ' + os.path.join(self.cfg.REFERENCE_DIR, name) + ' ' + refpath
    return self.execCmd(cmd)
  
  def runBWA(self, seqtype='single', input=[], ref='', output='', \
             option={'thread':8, 'checksr':True, 'refpath':None, 'addRG':False, 'rgroup':{}}):
    common.addPath(self.cfg.APPS_DIR)
    if not self.hasFai(option['refpath']):
      res = self.makeFai(option['refpath'])
      if not res:
        print(' Reference index (.fai) construction error.')
        return
    if not self.hasBWARefIndex(ref):
      res = self.makeBWARefIndex(option['refpath'], ref)
      if not res:
        print(' Reference transform error.')
        return
    cmd = 'bwa mem '
    if 'thread' in option:
      cmd += '-t '+str(option['thread'])
    if 'checksr' in option:
      cmd += ' -M'
    if 'addRG' in option and option['addRG']:
      cmd += ' -R "@RG\\t' + self.makeRGText(option['rgroup']) + '"'
    cmd += ' '+ref
    cmd += ' ' + input[0]
    if len(input) == 2 and os.path.exists(input[1]) and seqtype == 'paired':
      cmd += ' ' + input[1]
    cmd += ' > ' + output
    os.chdir(self.cfg.REFERENCE_DIR)
    return self.execCmd(cmd)

  def hasBowtRefIndex(self, refname):
    return os.path.exists(os.path.join(self.cfg.REFERENCE_DIR, refname + '.bt2'))

  def makeBowtRefIndex(self, refpath, refname, thread = 8):
    os.chdir(self.cfg.REFERENCE_DIR)
    cmd = 'bowtie2-build' + ' --threads ' + str(thread) + ' -f ' + refpath + ' ' + refname
    return self.execCmd(cmd)
  
  def runBowtie2(self, seqtype='single', input=[], ref='', output='', \
                 option={'thread':8, 'checksr':True, 'refpath':None, 'addRG':False, 'rgroup':{}}):
    if not self.hasFai(option['refpath']):
      res = self.makeFai(option['refpath'])
      if not res:
        print(' Reference index (.fai) construction error.')
        return
    if not self.hasBowtRefIndex(ref):
      res = self.makeBowtRefIndex(option['refpath'], ref, option['thread'])
      os.chdir(self.cfg.WORK_SPACE)
      if not res:
        print(' Reference index construction error.')
        return
    cmd = 'bowtie2'
    if 'addRG' in option and option['addRG']:
      cmd += ' --rg-id ' + option['rgroup']['ID']
      for key in option['rgroup']:
        if key != 'ID':
          cmd += ' --rg ' + key + ':' + option['rgroup'][key]
    if seqtype == 'single' and os.path.exists(input[0]):
      cmd += ' -U '+input[0]
    elif seqtype == 'paired-end' and len(input) == 2 and \
          os.path.exists(input[0]) and os.path.exists(input[1]):
            cmd += ' -1 '+input[0] + ' -2 '+input[1]
    if option['thread']:
      cmd += ' -p '+str(option['thread'])
    cmd += ' -x ' + ref + ' -S '+ output
    os.chdir(self.cfg.REFERENCE_DIR)
    return self.execCmd(cmd)

  def hasSTARRefIndex(self, refdir):
    return os.path.exists(refdir)

  def makeSTARRefIndex(self, refpath='', refname='', \
                    option={'thread':8, 'annotate':True, 'gtf':None}):
    cmd = 'STAR --runMode genomeGenerate' + \
      ' --genomeDir ' + os.path.join(self.cfg.REFERENCE_DIR, refname) + ' --genomeFastaFiles ' + refpath
    if option['annotate'] and os.path.exists(option['gtf']):
      cmd += ' --sjdbGTFfile ' + option['gtf']
    if option['thread']:
      cmd += ' --runThreadN ' + option['thread']
    return self.execCmd(cmd)

  def runSTAR(self, seqtype='single', input=[], ref='', output='', \
              option={'thread':8, 'annotate':True, 'gtf':None, 'refpath':None}):
    if not self.hasFai(option['refpath']):
      res = self.makeFai(option['refpath'])
      if not res:
        print(' Reference index (.fai) construction error.')
        return
    if not self.hasSTARRefIndex(os.path.join(self.cfg.REFERENCE_DIR, ref)):
      self.makeSTARRefIndex(refpath = option['refpath'], refname = ref, option = option)
    
    cmd = 'STAR --runMode genomeGenerate' + \
      ' --genomeDir ' + self.cfg.REFERENCE_DIR + ' --genomeFastaFiles ' + ref
    
    if option['annotate'] and os.path.exists(option['gtf']):
      cmd += ' --sjdbGTFfile ' + option['gtf']
    if option['thread']:
      cmd += ' --runThreadN ' + option['thread']
    return self.execCmd(cmd)

  def runSamtool2Fq(self, seqtype='single', input='', outdir='', outname = ''):
    if seqtype == 'single':
      cmd = 'samtools fastq -0 /dev/null ' + input + ' > ' + os.path.join(outdir, outname + '.fq')
    else:
      cmd = 'samtools collate -u -O ' + input + ' | samtools fastq '
      cmd += ' -1 ' + os.path.join(outdir, outname + '_1.fq')
      cmd += ' -2 ' + os.path.join(outdir, outname + '_2.fq')
      cmd += ' -0 /dev/null -s /dev/null -n'
    return self.execCmd(cmd)

  def runSamtool2BAM(self, input='',output='',option={}):
    cmd = 'samtools view'
    if 'thread' in option:
      cmd += ' -@ '+str(option['thread'])
    cmd += ' -b -o '+ os.path.join(self.cfg.OUT_DIR, output)+' ' + input
    return self.execCmd(cmd)

  def runSamtoolSort(self, input='', output='', option = {}):
    cmd = 'samtools sort -l1 -T tmp' + \
      ' -@ '+str(option['thread']) + ' -m '+str(option['ram']*1000)+'M' + \
      ' -O bam -o '+output + ' ' + input
    return self.execCmd(cmd)

  def runSamtoolIndex(self, input=''):
    cmd = 'samtools index '+ input
    return self.execCmd(cmd)

  def runPicardMD(self, input='', output='', metric = ''):
    common.addPath(self.cfg.APPS_DIR)
    cmd = 'java -jar '+os.path.join(self.cfg.APPS_DIR,'picard.jar') + ' MarkDuplicates'
    cmd += ' -I ' + input
    cmd += ' -O ' + output
    cmd += ' -M ' + metric
    return self.execCmd(cmd)

  def runTVC(self, input = '', output = '', ref = '', 
             option = { 'param' : '', 'motif' : '', 'thread':8, 'target':None, 'hotspot': None }):
    common.addPath(os.path.join(self.cfg.APPS_DIR,'TVC', 'bin'))
    cmd = 'python2 ' + os.path.join(self.cfg.APPS_DIR, 'TVC', 'bin', 'variant_caller_pipeline.py') + \
      ' --input-bam '+input + \
      ' --reference-fasta ' + ref + \
      ' --parameters-file ' + option['param'] + \
      ' --error-motifs ' + option['motif'] + \
      ' --generate-gvcf'
    if 'target' in option and option['target']:
      cmd += ' --region-bed ' +  option['target']
    if 'hotspot' in option and option['hotspot']:
      cmd += ' --hotspot-vcf ' +  option['hotspot']
    if 'control' in option and option['control']:
      cmd += ' --normal-bam ' +  option['control']
    if 'thread' in option and 1 < option['thread']:
      cmd += ' --num-threads ' + str(option['thread'])
    cmd += ' --output-dir ' + output
    return self.execCmd(cmd)

  def runBCFVarCall(self, input = '', output = '', ref = '', option = {}):
    cmd = 'bcftools mpileup -Ou' + \
      ' -f ' + ref + ' ' + input + \
       ' | bcftools call -vm -Oz -o ' + output + '.vcf.gz'
    return self.execCmd(cmd)
  
  def hasGATKRefDict(self, refpath = ''):
    fname, ext = os.path.splitext(os.path.basename(refpath))
    return os.path.exists(os.path.join(self.cfg.REFERENCE_DIR, fname + '.dict'))

  def makeGATKRefDict(self, refpath = ''):
    cmd = 'gatk CreateSequenceDictionary -R ' + refpath
    return self.execCmd(cmd)

  def hasGATKFeatureIndex(self, feature):
    return os.path.exists(feature + '.idx')
  
  def makeGATKFeatureIndex(self, feature):
    cmd = 'gatk IndexFeatureFile -I ' + feature
    return self.execCmd(cmd)

  def runGATKBRecal(self, input = '', output = '', ref = '', known = '',  \
                    option={'ram':8}):
    common.addPath(os.path.join(self.cfg.APPS_DIR,'gatk'))
    if not self.hasGATKRefDict(ref):
      res = self.makeGATKRefDict(ref)
      if not res:
        print(' Reference dictionary (.dict) construction error.')
    if not self.hasGATKFeatureIndex(known):
      res = self.makeGATKFeatureIndex(known)
      if not res:
        print(' Feature index construction error.')
    gatkcmd = 'gatk'
    if option['ram']:
      gatkcmd += ' --java-options "-Xmx' + str(option['ram'])+'g"'
    cmd = gatkcmd + ' BaseRecalibrator'
    cmd += ' -R ' + ref
    cmd += ' -I ' + input
    cmd += ' --known-sites ' + known
    cmd += ' -O ' + output+'_brecal.table'
    res = self.execCmd(cmd)
    if not res:
      return
    cmd = gatkcmd + ' ApplyBQSR'
    cmd += ' -R ' + ref
    cmd += ' -I ' + input
    cmd += ' -bqsr ' + output+'_brecal.table'
    cmd += ' -O ' + output
    return self.execCmd(cmd)

  def runGATKVarCall(self, input = '', output = '', ref = '', option={}):
    common.addPath(os.path.join(self.cfg.APPS_DIR, 'gatk'))
    if not self.hasGATKRefDict(ref):
      res = self.makeGATKRefDict(ref)
      if not res:
        print(' Reference dictionary (.dict) construction error.')
    cmd = 'gatk'
    if option['ram']:
      cmd += ' --java-options "-Xmx' + str(option['ram'])+'g"'
    cmd += ' HaplotypeCaller'
    if 'target' in option and os.path.exists(option['target']):
      cmd += ' -L ' +  option['target']
    cmd += ' -R ' + ref + \
      ' -I ' + input + \
      ' -O ' + output + \
      '.g.vcf.gz -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation'
    res = self.execCmd(cmd)
    if not res:
      print(' Failed to export gvcf.')
      return
    cmd = 'gatk'
    if option['ram']:
      cmd += ' --java-options "-Xmx' + str(option['ram'])+'g"'
    cmd += ' GenotypeGVCFs'
    cmd += ' -R ' + ref + \
      ' -V ' + output + '.g.vcf.gz' + \
      ' -O ' + output + '.vcf.gz'
    return self.execCmd(cmd)
    
  def runGATKVRecal(self, input = '', output = '', ref = '', \
                    option={'ram':8}):
    gatkcmd = 'gatk'
    if option['ram']:
      gatkcmd += ' --java-options "-Xmx' + str(option['ram'])+'g"'
    cmd = gatkcmd + ' VariantRecalibrator'
    cmd += ' -R ' + ref + \
     ' -V ' + input + \
    ' -O ' + output + '_snp.recal'
    cmd += ' -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP'
    cmd += ' --tranches-file ' + output + '_snp.tranches'
    res = self.execCmd(cmd)
    if not res:
      print(' Variant recalibration has failed.')
      return
    cmd = gatkcmd + ' ApplyVQSR'
    cmd += ' -V ' + input
    cmd += ' --recal-file ' + output + '_snp.recal'
    cmd += ' --tranches-file ' + output + '_snp.tranches'
    cmd += ' -O ' + output + '_snp.vcf'
    cmd += ' -mode SNP --create-output-variant-index true'
    return self.execCmd(cmd)
    
  def runGDVCall(self, input = '', output = '', ref = '', \
                  option={'gpu': False, 'processor':4, 'target': ''}):
    idir, iname = os.path.split(input)
    odir, oname = os.path.split(output)
    rdir, rname = os.path.split(ref)
    tdir = ''
    tname = ''
    if 'target' in option:
      tdir, tname = os.path.split(option['target'])
    cmd = 'sudo docker run'
    if 'gpu' in option and option['gpu']:
      cmd += ' --gpus 1'
    cmd += ' -v "' + rdir + '":"/REF_DIR"' 
    cmd += ' -v "' + idir + '":"/INPUT_DIR"' 
    cmd += ' -v "' + odir + '":"/OUTPUT_DIR"' 
    if tdir != '':
      cmd += ' -v "' + tdir + '":"/TARGET_DIR"' 
    cmd +=' google/deepvariant:"' + self.cfg.SOFTWARE_INFO['GDV']['ver'] + '"' + \
      ' /opt/deepvariant/bin/run_deepvariant --model_type=WGS' 
    cmd += ' --ref ' + os.path.join('/REF_DIR', rname)
    cmd += ' --reads ' + os.path.join('/INPUT_DIR', iname)
    if 'target' in option and os.path.exists(option['target']):
      cmd += ' --regions ' + os.path.join('/TARGET_DIR', tname)
    if 'processor' in option and 0 < option['processor']:
      cmd += ' --num_shards ' + str(option['processor'])
    cmd += ' --output_gvcf ' + os.path.join('/OUTPUT_DIR', oname+'.g')
    cmd += ' --output_vcf ' + os.path.join('/OUTPUT_DIR', oname)
    return self.execCmd(cmd)

  def runHTSeqCount(self, input = '', annotation = '', output = ''):
    cmd = 'htseq-count -r pos -t exon -f ' + input + ' ' + annotation + ' > ' + output
    return self.execCmd(cmd)
#htseq-count -f bam align.bam reference.gtf > count.txt
  def runCufflinks(self, input = '', annotation = '', output = ''):
    cmd = 'cufflinks'

    
    return self.execCmd(cmd)

  #-q -p 12 -m 100 -s 60 -G $GENE_GTF -M $RRNA_MASK \
 # --library-type fr-secondstrand --max-bundle-length 3500000 \
 # -o $OUT_CUFFLINKS --no-update-check \
 # ${fastqName}.STARBowtie2.bam"

  def runMACS2(self, input = '', control = None, output = '', species = '', genome = 0, 
              option = { 'bload' : True, 'lambda' : True, 'p-val' : -1, 'q-val': -1}):
    cmd = 'macs2 callpeak -f BAM -t ' + input
    if control :
      cmd = ' -c ' + control
    if len(species) : 
      cmd += ' -g ' + species
    else :
      cmd += ' -g ' + '{:.1E}'.format(genome)
    if not option['bload']:
      cmd += ' --broad'
    if not option['lamda']:
      cmd += ' --nolambda'
    if 'p-val' in option and -1 < option['p-val']:
      cmd += ' -p ' + '{:.1E}'.format(option['p-val'])
    if 'q-val' in option and -1 < option['q-val']:
      cmd += ' -q ' + '{:.1E}'.format(option['q-val'])
    cmd += ' --outdir ' + self.cfg.OUT_DIR + ' -n ' + output
    return self.execCmd(cmd)

#def run


