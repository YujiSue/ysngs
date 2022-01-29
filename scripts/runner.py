import os
import common
import installer
import subprocess
from subprocess import PIPE 

def execCmd(cmd):
  print('  Run: >', cmd)
  return subprocess.run(cmd, stdout=PIPE, stderr=PIPE, text=True, shell=True)

def downloadSRA(srid, output = '.', option = {'thread':8}):
  cmd = 'fasterq-dump ' + srid + ' -O ' + output
  if option['thread']:
    cmd += ' -e ' + str(option['thread'])
  return execCmd(cmd)

def runCutter(adaptor='', input = '', output = ''):
    cmd = 'cutadapt -b ' + adaptor + ' '+input + ' > ' + output
    return execCmd(cmd)

def hasBWARefIndex(ref):
  base, ext = os.path.splitext(ref)
  return os.path.exists(base+'.pac')

def makeBWARefIndex(refpath, name):
  if not os.path.exists(refpath):
    return (1, 'No reference.')
  cmd = 'bwa index -p ' + os.path.join(common.REFERENCE_DIR, name) + ' ' + refpath
  print('Run:\n > ', cmd)
  return execCmd(cmd)
  
def runBWA(seqtype='single', input=[], ref='', output='', \
           option={'thread':8, 'checksr':True, 'refpath':None, 'addRG':False, 'rgroup':''}):
  if not os.path.exists(input):
    return (1, 'No input.')
  common.addPath(common.APPS_DIR)
  if not (len(ref) and hasBWARefIndex(ref)):
    res = makeBWARefIndex(option['refpath'], ref)
    if res[0]:
      return res
  cmd = 'bwa mem '
  if 'thread' in option:
    cmd += '-t '+str(option['thread'])
  if 'checksr' in option:
    cmd += ' -M'
  cmd += ' '+ref
  if os.path.exists(input[0]):
    cmd += ' ' + input[0]
  else:
    return (1, 'File not found. "'+input[0]+'"')
  if len(input) == 2 and os.path.exists(input[1]) and seqtype == 'paired':
    cmd += ' ' + input[1]
  cmd += ' > ' + common.OUT_DIR + '/' + output + '.sam'
  print('Run: ', cmd)
  return execCmd(cmd)

def hasBowtRefIndex(ref):
   return os.path.exists(common.REFERENCE_DIR+'/'+ref+'.bt2')

def makeBowtRefIndex(refpath, ref):
  if not os.path.exists(refpath):
    return (1, 'No reference.')
  cmd = 'bowtie2-build -f ' + refpath + ' ' + ref
  print('Run:\n > ', cmd)
  return execCmd(cmd)
  
def runBowtie2(seqtype='single', input=[], ref='', output='', \
               option={'thread':8, 'checksr':True, 'refpath':None}):
  if len(input) == 0:
    return (1, 'No input.')
  common.addPath(common.APPS_DIR+'')
  if not (len(ref) and hasBowtRefIndex(ref)):
    res = makeBowtRefIndex(option['refpath'], ref)
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
  print('Run:\n > ', cmd)
  return execCmd(cmd)

def runSTAR(seqtype='single', input=[], ref='', output='', \
           option={'thread':8, 'annotate':True, 'gtf':None}):
  
  cmd = 'STAR --runMode genomeGenerate'+\
  ' --genomeDir ' + common.REFERENCE_DIR + ' --genomeFastaFiles ' + ref
  if option['annotate'] and os.path.exists(option['gtf']):
    cmd += ' --sjdbGTFfile ' + option['gtf']
  if option['thread']:
    cmd += ' --runThreadN ' + option['thread']
  print('Run:\n > ', cmd)
  return execCmd(cmd)

def runSamtool2Fq(input='', output=''):
  if not os.path.exists(input):
    return (1, 'File not found.')
  cmd = 'samtools fastq -b'
  return execCmd(cmd)


def runSamtool2BAM(input='',output='',option={}):
  if not os.path.exists(input):
    return (1, 'File not found.')
  cmd = 'samtools view'
  if 'thread' in option:
    cmd += ' -@ '+str(option['thread'])
  cmd += ' -b -o '+ os.path.join(common.OUT_DIR, output)+'.bam ' + input
  return execCmd(cmd)


def runSamtoolSort(input='', output='', \
                   option={'thread':8,'ram':'1000M'}):
  if not os.path.exists(input[0]):
    return (1, 'File not found.')
  cmd = 'samtools sort -l1 -T tmp'+\
  ' -@ '+str(option['thread'])+' -m '+option['ram'] +\
  ' -O bam -o '+output+'.bam'
  print('Run:\n > ', cmd)
  return execCmd(cmd)

def runSamtoolIndex(input=''):
  if not os.path.exists(input):
    return (1, 'File not found.')
  cmd = 'samtools index '+ input
  return execCmd(cmd)

def runPicardMD(input='', output='', metric = ''):
  if not os.path.exists(input):
    return (1, 'File not found.')
  cmd = 'Picard MarkDuplicates'
  cmd += ' -I ' + input
  cmd += ' -O ' + output
  cmd += ' -M ' + metric
  print('Run:\n > ', cmd)
  return execCmd(cmd)

def runTVC(input = '', outdir = '', ref = '', 
           option = { 'param' : '', 'motif' : ''}):
  if not os.path.exists(input):
    return (1, 'File not found.')
  if not os.path.exists(ref):
    return (1, 'Reference not found.')
  common.addPath(common.APPS_DIR+'/TSVC/bin')
  cmd = 'python2 '+common.APPS_DIR+'/TSVC/bin/variant_caller_pipeline.py' + \
  ' --input-bam '+input + \
  ' --reference-fasta ' + ref + \
  ' --parameters-file ' + common.PREFERENCE_DIR + '/' + option['param'] + \
  ' --error-motifs ' + common.PREFERENCE_DIR + '/' + option['motif'] + \
  ' --output-dir ' + outdir
  print('Run:\n > ', cmd)
  return execCmd(cmd)

def runBCFVarCall(input = '', output = '', ref = '', \
                  option = {}):
  if not os.path.exists(input):
    return (1, 'File not found.')
  if not os.path.exists(ref):
    return (1, 'Reference not found.')
  cmd = 'bcftools mpileup -Ou' + \
  ' -f ' + ref + ' ' + input + \
  ' | bcftools call -vm -O z - o ' + output + '.gz'
  return execCmd(cmd)

def runGATKBRecal(input = '', output = '', ref = '', \
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
  res = execCmd(cmd)
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
  return execCmd(cmd)

def runGATKVarCall(input = '', output = '', ref = '', \
                  option={'ram':'8g'}):
  if not os.path.exists(input):
    return (1, 'File not found.')
  if not os.path.exists(ref):
    return (1, 'Reference not found.')
  cmd = 'gatk'
  if option['ram']:
    cmd += ' --java-options "-Xmx'+option['ram']+'"'
  cmd += ' -R ' + common.REFERENCE_DIR + '/' + ref + \
  ' -I ' + input + \
  ' -O ' + output + '.g.vcf.gz -ERC GVCF'
  print('Run: ', cmd)
  return execCmd(cmd)
#   -G Standard \
#   -G AS_Standard
 
def runGATKVRecal(input = '', output = '', ref = '', \
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
  print('Run: ', cmd)
  res = execCmd(cmd)
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
  print('Run:\n > ', cmd)
  return execCmd(cmd)

def runGDVCall(input = '', output = '', ref = '', \
                option={'ram':'8g'}):
  if not os.path.exists(input):
    return (1, 'File not found.')
  if not os.path.exists(ref):
    return (1, 'Reference not found.')
  cmd = 'sudo docker run -v "' + common.WORK_SPACE + '":"/WORKSPACE"' + \
    ' google/deepvariant:"' + installer.SOFTWARE_INFO['GDV']['ver'] + '"' + \
    ' /opt/deepvariant/bin/run_deepvariant --model_type=WGS' 
  cmd += " --ref='" + ref + "'"  
  cmd += " --reads='"

  cmd += " --output_gvcf='" + output + "'"
  return execCmd(cmd)

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


