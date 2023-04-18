import os
import sys
import json
import subprocess
from subprocess import PIPE 
from ysngs import common
from ysngs import cmdgenerator as cmdgen

def getResult(res, prod = None, err = ''):
  if res[0]:
    return { 'state': 0, 'product': prod}
  else:
    return {'state': 1, 'err': err, 'msg': res[2],  'product': prod}

def runCustomRScipt(cfg, opts):
  res = common.runRScript(os.path.join(cfg.SCRIPT_DIR, opts['script']), None, opts['args'])
  return getResult(res, json.loads(res[1]), 'Failed to run R script.')
  
def runCustomScript(cfg, opts):
  res = common.execCmd(opts['cmd'] + ' ' + os.path.join(cfg.SCRIPT_DIR, opts['script']) + ' ' + opts['args'], verbose = opts['verbose'])
  return getResult(res, json.loads(res[1]), 'Failed to run script.')


# Download genome fasta
def downloadReference(cfg, opts):
  os.chdir(cfg.REFERENCE_DIR)
  res = common.execCmd(cmdgen.curlDownloadCmd({
    'url': opts['reference']['url']
    }), verbose=opts['verbose'])
  file = os.path.split(opts['reference']['url'])[1]
  ext = os.path.splitext(file)[1]
  print('file:', file)
  print('ext:', ext)
  if ext == '.gz':
    common.execCmd('gunzip ' + file)
  file = file[0:-len(ext)]
  print('file:', file)
  opts['reference']['path'] = os.path.join(cfg.REFERENCE_DIR, file)
  os.chdir(cfg.WORK_SPACE)
  return getResult(res, {'reference' : opts['reference']}, 'Failed to download.')

# SRA-toolkit
def downloadFromSRA(cfg, opts):
  os.chdir(cfg.TEMPORAL)
  common.addPath(os.path.join(cfg.APPS_DIR, 'sra'))
  res = common.execCmd(cmdgen.getSRACmd({
    'id': opts['sraid'],
    'outdir': cfg.SAMPLE_DIR
  }), verbose=opts['verbose'])
  os.chdir(cfg.SAMPLE_DIR)
  count = 0
  downloaded = []
  res = common.execCmd('find . -maxdepth 1 -name "'+opts['sraid']+'*.fq | wc -l')
  if int(res[1]) != 0:
    count += int(res[1])
    #res = common.execCmd('find . -maxdepth 1 -name "'+opts['sraid']+'*.fq')
    #downloaded.append(res[1].)
  res = common.execCmd('find . -maxdepth 1 -name "'+opts['sraid']+'*.fastq | wc -l')
  if int(res[1]) != 0:
    count += int(res[1])
  type = ''
  if count == 1:
    type = 'single'
  else:
    type = 'paired'
  os.chdir(cfg.WORK_SPACE)
  return getResult(res, {
    'raw': {
      'type': type,
      'path' : downloaded
    }
    })

def downloadTestData(cfg, ver):
  os.chdir(cfg.TEMPORAL)
  assert common.execCmd('curl -LO https://github.com/samtools/samtools/releases/download/' + ver + '/samtools-' + ver + '.tar.bz2')[0], 'Code download error.'
  assert common.execCmd('tar xvf ./samtools-' + ver + '.tar.bz2')[0], 'Expansion error.'
  assert common.execCmd('cp '+cfg.TEMPORAL+'/samtools-' + ver + '/examples/ex1.fa '+cfg.REFERENCE_DIR)[0]
  assert common.execCmd('gunzip '+cfg.TEMPORAL+'/samtools-' + ver + '/examples/ex1.sam.gz')[0]
  assert common.execCmd('cp '+cfg.TEMPORAL+'/samtools-' + ver + '/examples/ex1.sam ' +cfg.TEST_DIR)[0]
  assert common.execCmd('rm -r ./samtools*')[0]
  os.chdir(cfg.WORK_SPACE)
  assert common.execCmd('samtools view -b -T '+cfg.REFERENCE_DIR+'/ex1.fa -o '+cfg.TEST_DIR+'/ex1.bam '+cfg.TEST_DIR+'/ex1.sam')[0]
  assert common.execCmd('samtools fastq '+cfg.TEST_DIR+'/ex1.bam > '+cfg.TEST_DIR+'/ex1.fq')[0]


# Samtools
## Ref. index
def hasFai(ref):
  return os.path.exists(ref + '.fai')
def makeFai(cfg, opts):
  res = common.execCmd(cmdgen.getFaiCmd(opts), verbose=cfg.verbose)
  return getResult(res, {'faidx': opts['reference']['path'] + '.fai'})
## SAM/BAM => fastq
def runSamtool2Fq(cfg, opts):
  res = common.execCmd(cmdgen.getBam2FqCmd(opts), verbose=opts['verbose'] if 'verbose' in opts else False)
  return getResult(res, opts['output'])
def runSamtool2Fq_(cfg, opts):
  opts['raw']['path'] = []
  os.makedirs(os.path.join(cfg.SAMPLE_DIR, opts['sample']), exist_ok=True)
  if opts['raw']['type'] == 'single':
    opts['raw']['path'].append(os.path.join(cfg.OUT_DIR, opts['sample'], opts['export'] + '.fq'))
  else:
    opts['raw']['path'].append(os.path.join(cfg.OUT_DIR, opts['sample'], opts['export'] + '_1.fq'))
    opts['raw']['path'].append(os.path.join(cfg.OUT_DIR, opts['sample'], opts['export'] + '_2.fq'))
  res = common.execCmd(cmdgen.getBam2FqCmd({
    'type': opts['raw']['type'],
    'input': opts['aligned']['path'],
    'output': opts['raw']['path']
  }), verbose=opts['verbose'])
  return getResult(res, {'raw': opts['raw']})
## SAM => BAM
def runSamtool2BAM(cfg, opts):
  os.makedirs(os.path.join(cfg.SAMPLE_DIR, opts['sample']), exist_ok=True)
  out = os.path.join(cfg.OUT_DIR, opts['sample'], opts['export'])
  res = common.execCmd(cmdgen.getSam2BamCmd({
    'input': opts['aligned']['path'],
    'output': out,
    'thread' : opts['thread'] if 'thread' in opts else 1
  }), verbose=opts['verbose'])
  return getResult(res, {'aligned': { 'path': out }})
## Sort BAM
def runSamtoolSort(cfg, opts):
  res = common.execCmd(cmdgen.getBamSortCmd(cfg, opts), verbose=opts['verbose'])
  return getResult(res, {})
## Make BAM index
def runSamtoolIndex(cfg, opts):
  res = common.execCmd(cmdgen.getBaiCmd(cfg, opts), verbose=opts['verbose'])
  return getResult(res, {})



# Pre-process and quality control of fastq
## FastQC
def runFastQC(cfg, opts):
  common.addPath(os.path.join(cfg.APPS_DIR, 'FastQC'))
  os.makedirs(os.path.join(cfg.OUT_DIR, opts['sample']), exist_ok=True)
  res = common.execCmd(cmdgen.getQCCmd(cfg, opts), verbose=opts['verbose'])
  file = os.path.splitext(os.path.basename(opts['raw']['path']))[0]
  return getResult(res, {'fastqc':  os.path.join(cfg.OUT_DIR, opts['sample'], file + '_fastqc.html')})
## CutAdapt
def runCutter(cfg, opts):
  os.makedirs(os.path.join(cfg.OUT_DIR, opts['sample']), exist_ok=True)
  res = common.execCmd(cmdgen.getCutCmd(cfg, opts), verbose=opts['verbose'])
  return getResult(res, {'fqcut': '_cut.fq'})
## fastp
def runFastQFilter(cfg, opts):
  common.addPath(cfg.APPS_DIR)
  os.makedirs(os.path.join(cfg.OUT_DIR, opts['sample']), exist_ok=True)
  os.chdir(os.path.join(cfg.OUT_DIR, opts['sample']))
  res = common.execCmd(cmdgen.getFPCmd(cfg, opts), verbose=opts['verbose'])
  file = os.path.splitext(os.path.basename(opts['raw']['path']))[0]
  return getResult(res, {'fqfilter': os.path.join(cfg.OUT_DIR, opts['sample'], file + '_filtered.fq')})
  
# Mark duplicate by Picard
def runPicardMD(cfg, opts):
  res = common.execCmd(cmdgen.getMarkDPCmd(cfg, opts), verbose=opts['verbose'])
  return getResult(res, {})

# BWA
def makeBWARG(rgid, rg):
  rgopt = 'ID:' + rgid
  for key in rg:
    rgopt += '\\t' + key + ':' + rg[key]
  return rgopt
def hasBWARefIndex(cfg, opts):
  return os.path.exists(os.path.join(cfg.REFERENCE_DIR, opts['reference']['label']+'.pac'))
def makeBWARefIndex(cfg, opts):
  res = common.execCmd(cmdgen.getBWAiCmd(cfg, opts), verbose=opts['verbose'])
  return getResult(res, {})
def runBWA(cfg, opts):
  common.addPath(cfg.APPS_DIR)
  prod = {'intermediate': []}
  if not hasFai(opts['reference']['path']):
    res = makeFai(cfg, opts)
    if res['state'] == 1:
      return getResult(res)
  os.chdir(cfg.REFERENCE_DIR)
  if not hasBWARefIndex(cfg, opts):
    res = makeBWARefIndex(cfg, opts)
    if res['state'] == 1:
      return getResult(res)
  if 'rgid' in opts:
    opts['rgroup'] = makeBWARG(opts['rgid'], opts['rg'])
  os.makedirs(os.path.join(cfg.OUT_DIR, opts['sample']), exist_ok=True)
  opts['output'] = os.path.join(cfg.OUT_DIR, opts['sample'], opts['export'] + '_mapped.sam')
  res = common.execCmd(cmdgen.getBWACmd(cfg, opts), verbose=opts['verbose'])
  prod['intermediate'].append(opts['output'])
  opts['input'] = opts['output']
  opts['output'] = os.path.join(cfg.OUT_DIR, opts['sample'], opts['export'] + '_mapped.bam')
  runSamtool2BAM(cfg, opts)
  prod['intermediate'].append(opts['output'])
  opts['input'] = opts['output']
  opts['output'] = os.path.join(cfg.OUT_DIR, opts['sample'], opts['export'] + '_sorted.bam')
  runSamtoolSort(cfg, opts)
  opts['input'] = opts['output']
  opts['output'] = os.path.join(cfg.OUT_DIR, opts['sample'], opts['export'] + '.bam')
  if opts['markdp']:
    prod['intermediate'].append(opts['output'])
    opts['metric'] = os.path.join(cfg.OUT_DIR, opts['sample'], opts['export'] + '_metric.txt')
    runPicardMD(cfg, opts)
  else:
    common.execCmd('mv ' + opts['input'] + ' ' + opts['output'], verbose=opts['verbose'])  
  opts['input'] = opts['output']
  prod['aligned'] = opts['output']
  runSamtoolIndex(cfg, opts)
  os.chdir(cfg.WORK_SPACE)
  return getResult(res, prod)

# Bowtie2
def hasBowtRefIndex(cfg, opts):
  return os.path.exists(os.path.join(cfg.REFERENCE_DIR, opts['reference']['label'] + '.bt2'))
def makeBowtRefIndex(cfg, opts):
  os.chdir(cfg.REFERENCE_DIR)
  res = common.execCmd(cmdgen.getBowtiCmd(cfg, opts), verbose=opts['verbose'])
  return getResult(res, {})
def runBowtie2(cfg, opts):
  os.chdir(cfg.WORK_SPACE)
  if not hasFai(input['reference']['path']):
    res = makeFai(cfg, opts)
    if res['state'] == 1:
      return getResult(res, {}, 'Reference index (.fai) construction error.')
  if not hasBowtRefIndex(cfg, opts):
    res = makeBowtRefIndex(cfg, opts)
    if res['state'] == 1:
      return getResult(res, {}, 'Reference index construction error.')
  os.chdir(cfg.REFERENCE_DIR)
  res = common.execCmd(cmdgen.getBowtCmd(cfg, opts), verbose=opts['verbose'])
  return getResult(res, {})

# STAR
def hasSTARRefIndex(refdir):
  return os.path.exists(refdir)
def makeSTARRefIndex(cfg, input, v):
  res = common.execCmd(cmdgen.getSTARiCmd(cfg, input), verbose=v)
  return getResult(res, {})
def runSTAR(cfg, input, v):
  common.addPath(cfg.APPS_DIR)
  if not hasFai(input['reference']['path']):
    res = makeFai(cfg, input, v)
    if res['state'] == 1:
      return getResult(res, {}, 'Reference index (.fai) construction error.')
  if not hasSTARRefIndex(cfg, input, v):
    res = makeSTARRefIndex(cfg, input, v)
    if res['state'] == 1:
      return getResult(res, {}, 'Reference index construction error.')
  res = common.execCmd(cmdgen.getSTARCmd(cfg, input), verbose=v)
  return getResult(res, {})

# HISAT2
def hasHISATRefIndex(ref):
  return os.path.exists(ref + '.1.ht2')
def makeHISATRefIndex(cfg, opts, v):
  res = common.execCmd(cmdgen.getHISATiCmd(cfg, opts), verbose=v)
  return getResult(res, {})
def runHISAT(cfg, opts, v):
  common.addPath(cfg.APPS_DIR)
  if not hasHISATRefIndex(cfg, opts, v):
    res = makeHISATRefIndex(cfg, opts, v)
    if res['state'] == 1:
      return getResult(res, {}, 'Reference index construction error.')
  res = common.execCmd(cmdgen.getHISATCmd(cfg, opts), verbose=v)
  return getResult(res, {})

#cmd = 'hisat2-build ' + reference['path'] + ' ' + reference['label']
#cmd = 'hisat2'
# -x TAIR10 -U ${seqlib}.clean.fastq -S ${seqlib}.sam


# RNA-seq
## htseq-count
def runHTSeqCount(cfg, input, v):
  res = common.execCmd(cmdgen.getHTSCmd(cfg, input), verbose=v)
  return getResult(res, {})

## FeatureCount
def runFeatureCount(cfg, input = '', annotation = '', output = '',  option = {}):
  os.chdir(cfg.WORK_SPACE)
  
  return execCmd(cmd)

## Cufflinks
def runCufflinks(cfg, input, v):
  common.addPath(os.path.join(cfg.APPS_DIR, 'cuff'))
  res = common.execCmd(cmdgen.getCuffLCmd(cfg, input), verbose=v)
  return getResult(res, {})
  
## CuffDiff
def runCuffDiff(cfg, input, v):
  common.addPath(os.path.join(self.cfg.APPS_DIR, 'cuff'))
  res = common.execCmd(cmdgen.getCuffDCmd(cfg, input), verbose=v)
  return getResult(res, {})
## CuffMerge
def runCuffMerge(cfg, input, v):
  common.addPath(os.path.join(self.cfg.APPS_DIR, 'cuff'))
  res = common.execCmd(cmdgen.getCuffMCmd(cfg, input), verbose=v)
  return getResult(res, {})

## CummeR

## EdgeR
def runEdgeR(cfg, script = '', output = '', args = []):
  os.chdir(cfg.WORK_SPACE)
  return runRScript(script, output, args)

def removeIntermediateFiles(cfg, opts):
  res = {}
  for file in opts['intermediate']:
    if os.path.exists(file):
      res = common.execCmd('rm -r ' + file)
  return getResult(res, opts)

class Analyzer:
  def __init__(self, config):
    self.cfg = config
    self.result = {'intermediate': []}
  
  def check(self, proc):
    print(proc)
    if proc['state'] == 0:
      self.result.update(proc['product'])
    else:
      print('Error', '#'+str(proc['err']), proc['msg'])
      sys.exit(1)

  def qc(self, opts):
    result = {}
    if opts['fastqc']:
      proc = runFastQC(self.cfg, opts)
      self.check(proc)
    if opts['fqcut']:
      proc = runCutter(self.cfg, opts)
      self.check(proc)
      opts['raw']['ori'] = opts['raw']['path']
      opts['raw']['path'] = proc['product']['fqcut']
    if opts['fqfilter']:
      proc = runFastQFilter(self.cfg, opts)
      self.check(proc)
      if 'ori' not in opts['raw']:
        opts['raw']['ori'] = opts['raw']['path']
      opts['raw']['path'] = proc['product']['fqfilter']
      self.result['intermediate'].append(proc['product']['fqfilter'])

  def mapping(self, opts):
    mapapp = opts['app']
    if mapapp == 'bwa':
      runBWA(self.cfg, opts)
    elif mapapp == 'bowtie':
      runBowtie2(self.cfg, opts)
    elif mapapp == 'star':
      runSTAR(self.cfg, opts)
    elif mapapp == 'hisat':
      runHISAT(self.cfg, opts)
  
  def varcall(self, opts):
    vcapp = opts['app']
    if vcapp == 'bcf':
      runBCFVarCall(self.cfg, opts)
    elif vcapp == 'gatk':
      runGATKVarCall(self.cfg, opts)
    elif vcapp == 'tvc':
      runTVC(self.cfg, opts)
    elif vcapp == 'ivc':
      runStrelka(self.cfg, opts)
    elif vcapp == 'gdv':
      runGDVCall(self.cfg, opts)
  

  def postvc(self, opts):
    if 'hlatype' in opts:
      runHLAtyping()
  
  def count(cfg, opts):
    cntapp = opts['count']['app'] 

  def transcriptome(cfg, opts):
    degapp = opts['deg']['app']

  def peaks(cfg, opts):
    peakapp = opts['peak']['app']

  def clean(self):
    proc = removeIntermediateFiles(self.cfg, self.result)
    self.check(proc)




      
def runHLAtyping():
  cmd = 'hla-da'

  
# GATK
## Reference dict.
def hasGATKRefDict(ref):
  return os.path.exists(ref + '.dict')
def makeGATKRefDict(cfg, input, v):
  common.addPath(os.path.join(cfg.APPS_DIR, 'gatk'))
  os.chdir(cfg.REFERENCE_DIR)
  res = common.execCmd(cmdgen.getGATKRefiCmd(cfg, input), verbose=v)
  os.chdir(cfg.WORK_SPACE)
  return getResult(res, {})
## Feature index
def hasGATKFeatureIndex(cfg, feature):
  return os.path.exists(feature + '.idx')
def makeGATKFeatureIndex(cfg, feature):
  cmd = 'gatk IndexFeatureFile -I ' + feature
  return execCmd(cmd)
## Base recalibration
def runGATKBRecal(cfg, input = '', output = '', ref = '', known = '', option={}):
  common.addPath(os.path.join(cfg.APPS_DIR,'gatk'))
  os.chdir(cfg.WORK_SPACE)
  if not hasGATKRefDict(ref):
    res = makeGATKRefDict(ref)
    if not res:
      print(' Reference dictionary (.dict) construction error.')
  if not hasGATKFeatureIndex(known):
    res = makeGATKFeatureIndex(known)
    if not res:
      print(' Feature index construction error.')
  gatkcmd = 'gatk'
  if 'ram' in option:
    gatkcmd += ' --java-options "-Xmx' + str(option['ram'])+'g"'
  cmd = gatkcmd + ' BaseRecalibrator'
  cmd += ' -R ' + ref
  cmd += ' -I ' + input
  cmd += ' --known-sites ' + known
  cmd += ' -O ' + output+'_brecal.table'
  res = execCmd(cmd)
  if not res:
    return
  cmd = gatkcmd + ' ApplyBQSR'
  cmd += ' -R ' + ref
  cmd += ' -I ' + input
  cmd += ' -bqsr ' + output+'_brecal.table'
  cmd += ' -O ' + output
  return execCmd(cmd)

# Variant call
## BCFtools
def runBCFVarCall(cfg, opts):
  res = common.execCmd(cmdgen.getBCFVCCmd(cfg, opts), verbose=opts['verbose'])
  return getResult(res, {})
## GATK
def runGATKVarCall(cfg, input, v):
  common.addPath(os.path.join(cfg.APPS_DIR, 'gatk'))
  os.chdir(cfg.REFERENCE_DIR)
  if not hasGATKRefDict(ref):
    res = makeGATKRefDict(cfg, input, v)
    if res['state'] == 1:
      return getResult(res, {}, 'Reference index construction error.')
  res = common.execCmd(cmdgen.getGATKVCCmd(cfg, input), verbose=v)
  if res['state'] == 1:
    return getResult(res, {}, 'Failed to export gvcf.')
  res = common.execCmd(cmdgen.getGATKGenCmd(cfg, input), verbose=v)
  return getResult(res, {}, 'Failed to convert gvcf => vcf.')

## TVC
def runTVC(cfg, input, v):
  common.addPath(os.path.join(cfg.APPS_DIR, 'TVC', 'bin'))
  res = common.execCmd(cmdgen.getTVCCmd(cfg, input), verbose=v)
  return getResult(res, {}, 'Failed to run Torrent VariantCaller.')

## IVC
def runStrelka(cfg, input, v):
  common.addPath(os.path.join(cfg.APPS_DIR, 'strelka'))
  res = common.execCmd(cmdgen.getIVCCmd(cfg, input), verbose=v)
  return getResult(res, {}, 'Failed to run Strelka.')

## GDV
def runGDVCall(cfg, input, v):
  res = common.execCmd(cmdgen.getGDVCmd(cfg, input), verbose=v)
  return getResult(res, {})

## Variant recalibration by GATK
def runGATKVRecal(cfg, input = '', output = '', ref = '', resources = [], option={}):
  common.addPath(os.path.join(cfg.APPS_DIR, 'gatk'))
  os.chdir(cfg.WORK_SPACE)
  gatkcmd = 'gatk'
  if 'ram' in option:
    gatkcmd += ' --java-options "-Xmx' + str(option['ram'])+'g"'
  cmd = gatkcmd + ' VariantRecalibrator'
  cmd += ' -R ' + ref + \
    ' -V ' + input + \
    ' -O ' + output + '_snp.recal'
  for resource in resources:
    cmd += ' --resource ' + resource
  cmd += ' -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP'
  cmd += ' --tranches-file ' + output + '_snp.tranches'
  res = execCmd(cmd)
  if not res:
    print(' Variant(SNP) recalibration has failed.')
    return
  cmd = gatkcmd + ' ApplyVQSR'
  cmd += ' -V ' + input
  cmd += ' --recal-file ' + output + '_snp.recal'
  cmd += ' --tranches-file ' + output + '_snp.tranches'
  cmd += ' -O ' + output + '_snp.vcf'
  cmd += ' -mode SNP -truth-sensitivity-filter-leve 99.5 --create-output-variant-index true'
  res = execCmd(cmd)
  if not res:
    print(' Variant(SNP) recalibration apply has failed.')
    return
  cmd = gatkcmd + ' VariantRecalibrator'
  cmd += ' -R ' + ref + \
  ' -V ' + input + \
  ' -O ' + output + '_indel.recal'
  for resource in resources:
    cmd += ' --resource ' + resource
  cmd += ' -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR  -mode INDEL'
  cmd += ' --tranches-file ' + output + '_indel.tranches'
  res = execCmd(cmd)
  if not res:
    print(' Variant(InDel) recalibration has failed.')
    return
  cmd = gatkcmd + ' ApplyVQSR'
  cmd += ' -V ' + input
  cmd += ' --recal-file ' + output + '_indel.recal'
  cmd += ' --tranches-file ' + output + '_indel.tranches'
  cmd += ' -O ' + output + '_indel.vcf'
  cmd += ' -mode INDEL -truth-sensitivity-filter-leve 99.0 --create-output-variant-index true'
  res = execCmd(cmd)
  if not res:
    print(' Variant(InDel) recalibration apply has failed.')
    return
  cmd = gatkcmd + ' GatherVcfs -R ' + ref
  cmd += ' ' + output + '_indel.vcf ' + output + '_snp.vcf'
  cmd += ' -O ' + output + '_recal.vcf'
  return execCmd(cmd)
  
  




  def runMACS2(self, input = '', control = None, output = '', species = '', genome = 0, 
              option = { 'bload' : True, 'lambda' : True, 'p-val' : -1, 'q-val': -1}):
    os.chdir(self.cfg.WORK_SPACE)
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




