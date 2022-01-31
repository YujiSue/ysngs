import os
import subprocess
from subprocess import PIPE 
import ysngs.common as common

class apprun:
  def __init__(self, config):
    self.cfg = config
  
  def printMsgLine(self, msg = ''):
    lines = msg.splitlines()
    for l in lines:
      print(' >', l)

  def execCmd(self, cmd):
    print('Run: >', cmd)
    proc = subprocess.run(cmd, stdout=PIPE, stderr=PIPE, text=True, shell=True)
    if proc.returncode:
      print('Error:')
      self.printMsgLine(proc.stderr)
      return False
    else:
      self.printMsgLine(proc.stdout)
      return True
    
  def downloadSRA(self, srid, output = '.', option = {'thread':8}):
    cmd = 'fasterq-dump ' + srid + ' -O ' + output
    if option['thread']:
      cmd += ' -e ' + str(option['thread'])
    return self.execCmd(cmd)

  def runCutter(self, adaptor='', input = '', output = ''):
      cmd = 'cutadapt -b ' + adaptor + ' '+input + ' > ' + output
      return self.execCmd(cmd)

  def hasBWARefIndex(self, ref):
    base, ext = os.path.splitext(ref)
    return os.path.exists(base+'.pac')

  def makeBWARefIndex(self, refpath, name):
    cmd = 'bwa index -p ' + os.path.join(self.cfg.REFERENCE_DIR, name) + ' ' + refpath
    return self.execCmd(cmd)
  
  def runBWA(self, seqtype='single', input=[], ref='', output='', \
             option={'thread':8, 'checksr':True, 'refpath':None, 'addRG':False, 'rgroup':''}):
    for f in input:
      if not os.path.exists(f):
        print('  ', f, ' is not found.')
        return False
    common.addPath(self.cfg.APPS_DIR)
    if not self.hasBWARefIndex(ref):
      res = self.makeBWARefIndex(option['refpath'], ref)
      if not res:
        print('Reference index construction error.')
        return
    cmd = 'bwa mem '
    if 'thread' in option:
      cmd += '-t '+str(option['thread'])
    if 'checksr' in option:
      cmd += ' -M'
    if 'addRG' in option and option['addRG']:
      cmd += ' -R "@RG\\t' + option['rgroup'] + '"'
    cmd += ' '+ref
    cmd += ' ' + input[0]
    if len(input) == 2 and os.path.exists(input[1]) and seqtype == 'paired':
      cmd += ' ' + input[1]
    cmd += ' > ' + output
    os.chdir(self.cfg.REFERENCE_DIR)
    return self.execCmd(cmd)

  def hasBowtRefIndex(self, ref):
    return os.path.exists(self.cfg.REFERENCE_DIR+'/'+ref+'.bt2')

  def makeBowtRefIndex(self, refpath, ref):
    cmd = 'bowtie2-build -f ' + refpath + ' ' + ref
    return self.execCmd(cmd)
  
  def runBowtie2(self, seqtype='single', input=[], ref='', output='', \
                 option={'thread':8, 'checksr':True, 'refpath':None}):
    self.cfg.addPath(self.cfg.APPS_DIR+'')
    if not (len(ref) and self.hasBowtRefIndex(ref)):
      res = self.makeBowtRefIndex(option['refpath'], ref)
      if res[0]:
        return res
    cmd = 'bowtie2'
    if seqtype == 'single' and os.path.exists(input[0]):
      cmd += ' -U '+input[0]
    elif seqtype == 'paired-end' and len(input) == 2 and \
          os.path.exists(input[0]) and os.path.exists(input[1]):
            cmd += ' -1 '+input[0] + ' -2 '+input[1]
    if option['thread']:
      cmd += ' -p '+str(option['thread'])
    cmd += ' -x ' + ref + ' -S '+ output + '.sam'
    return self.execCmd(cmd)

  def runSTAR(self, seqtype='single', input=[], ref='', output='', \
              option={'thread':8, 'annotate':True, 'gtf':None}):
    cmd = 'STAR --runMode genomeGenerate' + \
      ' --genomeDir ' + self.cfg.REFERENCE_DIR + ' --genomeFastaFiles ' + ref
    if option['annotate'] and os.path.exists(option['gtf']):
      cmd += ' --sjdbGTFfile ' + option['gtf']
    if option['thread']:
      cmd += ' --runThreadN ' + option['thread']
    return self.execCmd(cmd)

  def runSamtool2Fq(self, input='', output=''):
    cmd = 'samtools fastq -b'
    return self.execCmd(cmd)

  def runSamtool2BAM(self, input='',output='',option={}):
    cmd = 'samtools view'
    if 'thread' in option:
      cmd += ' -@ '+str(option['thread'])
    cmd += ' -b -o '+ os.path.join(self.cfg.OUT_DIR, output)+' ' + input
    return self.execCmd(cmd)

  def runSamtoolSort(self, input='', output='', \
                     option={'thread':8,'ram':'1000M'}):
    cmd = 'samtools sort -l1 -T tmp' + \
      ' -@ '+str(option['thread']) + ' -m '+option['ram'] + \
      ' -O bam -o '+output + ' ' + input
    return self.execCmd(cmd)

  def runSamtoolIndex(self, input=''):
    cmd = 'samtools index '+ input
    return self.execCmd(cmd)

  def runPicardMD(self, input='', output='', metric = ''):
    common.addPath(common.APPS_DIR)
    cmd = 'java -jar '+os.path.join(common.APPS_DIR,'picard.jar') + ' MarkDuplicates'
    cmd += ' -I ' + input
    cmd += ' -O ' + output
    cmd += ' -M ' + metric
    return self.execCmd(cmd)

  def runTVC(self, input = '', outdir = '', ref = '', 
             option = { 'param' : '', 'motif' : ''}):
    common.addPath(os.path.join(common.APPS_DIR,'TVC', 'bin'))
    cmd = 'python2 '+self.cfg.APPS_DIR+'/TSVC/bin/variant_caller_pipeline.py' + \
      ' --input-bam '+input + \
      ' --reference-fasta ' + ref + \
      ' --parameters-file ' + self.cfg.PREFERENCE_DIR + '/' + option['param'] + \
      ' --error-motifs ' + self.cfg.PREFERENCE_DIR + '/' + option['motif'] + \
      ' --output-dir ' + outdir
    return self.execCmd(cmd)

  def runBCFVarCall(self, input = '', output = '', ref = '', option = {}):
    cmd = 'bcftools mpileup -Ou' + \
      ' -f ' + ref + ' ' + input + \
       ' | bcftools call -gvm -Oz -o ' + output
    return self.execCmd(cmd)

  def makeGATKRefDict(self, ref = ''):
    cmd = 'gatk CreateSequenceDictionary -R '+ref
    return self.execCmd(cmd)

  def runGATKBRecal(self, input = '', output = '', ref = '', known = '',  \
                    option={'ram':'8g'}):
    common.addPath(os.path.join(common.APPS_DIR,'gatk'))
    gatkcmd = 'gatk'
    if option['ram']:
      gatkcmd += ' --java-options "-Xmx'+option['ram']+'"'
    cmd = gatkcmd + ' BaseRecalibrator'
    cmd += ' -R ' + ref
    cmd += ' -I ' + input
    cmd += ' --known-sites ' + known
    cmd += ' -O ' + output+'_brecal.table'
#-L $i-scattered.interval_list
#$knownSiteArg -O $out
    res = self.execCmd(cmd)
    if 0 < res.returncode:
      return res.stderr
    cmd = gatkcmd + ' ApplyBQSR'
    cmd += ' -R ' + ref
    cmd += ' -I ' + input
    cmd += ' -bqsr ' + output+'_brecal.table'
    cmd += ' -O ' + output
#-L $i-scattered.interval_list -bqsr $bqfile \
#--static-quantized-quals 10 --static-quantized-quals 20 \
#--static-quantized-quals 30 -O $outp
    return self.execCmd(cmd)

  def runGATKVarCall(self, input = '', output = '', ref = '', option={'ram':'8g'}):
    common.addPath(os.path.join(self.cfg.APPS_DIR,'gatk'))
    self.makeGATKRefDict(ref)
    cmd = 'gatk'
    if option['ram']:
      cmd += ' --java-options "-Xmx'+option['ram']+'"'
    cmd += ' HaplotypeCaller'
    cmd += ' -R ' + self.cfg.REFERENCE_DIR + '/' + ref + \
      ' -I ' + input + \
      ' -O ' + output + ' -ERC GVCF -G Standard -G AS_Standard'
    return self.execCmd(cmd)
 
  def runGATKVRecal(self, input = '', output = '', ref = '', \
                    option={'ram':'8g'}):
    gatkcmd = 'gatk'
    if option['ram']:
      gatkcmd += ' --java-options "-Xmx'+option['ram']+'"'
    cmd = gatkcmd + ' BaseRecalibrator'

    cmd += ' -R ' + ref
    cmd += ' -I ' + input
    cmd += ' -O ' + output+'_brecal.table'
#-L $i-scattered.interval_list
#$knownSiteArg -O $out
    res = self.execCmd(cmd)
    if res[0]:
      return res
    cmd = gatkcmd + ' ApplyBQSR'
    cmd += ' -R ' + ref
    cmd += ' -I ' + input
    cmd += ' -bqsr ' + output+'_brecal.table'
    cmd += ' -O ' + output
#-L $i-scattered.interval_list -bqsr $bqfile \
#--static-quantized-quals 10 --static-quantized-quals 20 \
#--static-quantized-quals 30 -O $outp
    return self.execCmd(cmd)

  def runGDVCall(self, input = '', output = '', ref = '', \
                  option={'ram':'8g'}):
    cmd = 'sudo docker run -v "' + self.cfg.WORK_SPACE + '":"/WORKSPACE"' + \
      ' google/deepvariant:"' + self.cfg.SOFTWARE_INFO['GDV']['ver'] + '"' + \
      ' /opt/deepvariant/bin/run_deepvariant --model_type=WGS' 
    cmd += " --ref='" + ref + "'"  
    cmd += " --reads='"
    cmd += " --output_gvcf='" + output + "'"
    return self.execCmd(cmd)

#sudo docker run -v "/home/phys2":"/HOME" 
# google/deepvariant:"${BIN_VERSION}" 
# /opt/deepvariant/bin/run_deepvariant 
# --model_type=WGS 
# --ref='/HOME/Reference/reference.fasta'  
# --reads='/HOME/Test/test.bam'  

def runHTSeqCount():
  return

def runCufflinks():
  return

def runCallPeaks():
  return

#def run


