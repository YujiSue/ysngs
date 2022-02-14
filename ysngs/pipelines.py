from common import config
from run import apprun
import getopt
import os
import sys

def preprocess(pref, app):
  # Cut adapter
  if pref['cutadapt'] and pref['iformat'] == 'fastq':
    for f in pref['inputs']:
      name, ext = os.path.splitext(f)
      ofile = name + '_cut' + ext
      app.runCutter(adaptor = pref['adapter_seq'], site = pref['adapter_site'], input = f, output = ofile)
      pref['mediates'].append(ofile)
      f = ofile
  if pref['fqfliter'] and pref['iformat'] == 'fastq':
    for f in pref['inputs']:
      name, ext = os.path.splitext(f)
      ofile = name + '_filtered' + ext
      app.runFastQFilter(input = f, output = ofile, param = {
        'min_qual': pref['fq_qual'], 'min_len': pref['fq_minlen'], 'min_complex': pref['fq_complex'], 'thread': pref['thread_num']
      })
      pref['mediates'].append(ofile)
      f = ofile

def mapping(pref, app):
  pref['rg_info'] = makeReadGroupInfo(pref)
  
  # BAM => Fastq
  if pref['iformat'] == 'bam':
    app.runSamtool2Fq(seqtype = pref['seqtype'], input = pref['inputs'][0], outdir = app.cfg.OUT_DIR, outname = pref['output_name'] + '_raw')
    if pref['seqtype'] == 'single':
      pref['inputs'][0] = os.path.join(app.cfg.OUT_DIR, pref['output_name'] + '_raw.fq')
    elif pref['seqtype'] == 'paired':
      pref['inputs'][0]= os.path.join(app.cfg.OUT_DIR, pref['output_name'] + '_raw_1.fq')
      pref['inputs'].append(os.path.join(app.cfg.OUT_DIR, pref['output_name'] + '_raw_2.fq'))
    for f in pref['inputs']:
      pref['mediates'].append(f)
  
  # Fastq => SAM
  pref['output'] = os.path.join(app.cfg.OUT_DIR, pref['output_name'] + '_mapped.sam')
  if len(pref['refname']) == 0:
      pref['refname'], ext = os.path.splitext(os.path.basename(pref['reference']))
  
  # BWA-MEM
  if pref['mapper'] == 'BWA':
    app.runBWA(seqtype = pref['seqtype'], input = pref['inputs'], ref = pref['refname'], output = pref['output'], \
             option={'thread':pref['thread_num'], 'checksr':pref['detect_sr'], 'refpath':pref['reference'], 'addRG':pref['add_rg'], 'rgroup':pref['rg_info']})
  # Bowtie2
  elif pref['mapper'] == 'Bowtie2':
    app.runBowtie2(seqtype = pref['seqtype'], input = pref['inputs'], ref = pref['refname'], output = pref['output'], \
             option={'thread':pref['thread_num'], 'checksr':pref['detect_sr'], 'refpath':pref['reference'], 'addRG':pref['add_rg'], 'rgroup':pref['rg_info']})
  # STAR
  
  pref['mediates'].append(pref['output'])
  pref['input'] = pref['output']
  os.chdir(app.cfg.WORK_SPACE)

  # SAM -> BAM
  pref['output'] = os.path.join(app.cfg.OUT_DIR, pref['output_name'] + '_mapped.bam')
  app.runSamtool2BAM(input = pref['input'], output = pref['output'], option={
    'thread':pref['thread_num']
    })
  pref['mediates'].append(pref['output'])

  # Sort BAM
  pref['input'] = pref['output']
  pref['output'] = os.path.join(app.cfg.OUT_DIR, pref['output_name'] + '_sorted.bam')
  app.runSamtoolSort(input = pref['input'], output = pref['output'], option={
    'thread':pref['thread_num'], 'ram': pref['use_ram']/4
    })
  pref['mediates'].append(pref['output'])
  pref['input'] = pref['output']

  # Mark duplicate
  if pref['mark_dup']:
    pref['output'] = os.path.join(app.cfg.OUT_DIR, pref['output_name'] + '_marked.bam')
    pref['metric'] = os.path.join(app.cfg.OUT_DIR, pref['output_name'] + '_metric.txt')
    app.runPicardMD(input = pref['input'], output = pref['output'], metric = pref['metric'])
    pref['mediates'].append(pref['output'])
    pref['input'] = pref['output']

  # Recalibration sequence quality
  if pref['recal_seq']:
    pref['output'] = os.path.join(app.cfg.OUT_DIR, pref['output_name'] + '.bam')
    app.runGATKBRecal(input = pref['input'], output = pref['output'], ref = pref['reference'], known = pref['known_site'], option = {'ram': pref['use_ram'] })
  else:
    pref['output'] = os.path.join(app.cfg.OUT_DIR, pref['output_name'] + '.bam')
    os.system('mv ' + pref['input'] + ' ' + pref['output'])
  pref['input'] = pref['output']

  # Make index
  app.runSamtoolIndex(input = pref['input'])
  
def cleanup(pref):
  for f in pref['mediates']:
    os.system('rm ' + f)

def varcall(pref, app):
  # Preprocessing fastq
  preprocess(pref, app)
  
  # Mapping fastq
  if pref['mapping']:
    mapping(pref, app)
  else:
    pref['input'] = pref['inputs'][0]
  
  # Variant call
  if pref['recal_var']:
    pref['output']  = os.path.join(app.cfg.OUT_DIR, pref['output_name'] + '_raw')
  else:
    pref['output'] = os.path.join(app.cfg.OUT_DIR, pref['output_name'])
  # GATK HaplotypeCaller
  if pref['vcaller'] == 'GATK':
    app.runGATKVarCall(input = pref['input'], output = pref['output'], ref = pref['reference'], option = {'ram': pref['use_ram']})
    pref['input'] = pref['output'] + '.vcf.gz'
  # BCFtools mpileup and call
  elif pref['vcaller'] == 'BCF':
    app.runBCFVarCall(input = pref['input'], output = pref['output'], ref = pref['reference'])
    pref['input'] = pref['output'] + '.vcf.gz'
  # TorrentSuite variant caller
  elif pref['vcaller'] == 'TVC':
    tvcopt = {
      'thread': pref['thread_num'],
      'param' : os.path.join(app.cfg.PREFERENCE_DIR, pref['vcparam']), 
      'motif' : os.path.join(app.cfg.PREFERENCE_DIR, pref['vcmotif']) 
    }
    if pref['targeted_seq']:
      tvcopt['target'] = pref['target']
    if pref['use_hotspot']:
      tvcopt['hotspot'] = pref['known_site']
    if pref['has_control']:
      tvcopt['control'] = pref['control_input']
    app.runTVC(input = pref['input'], output = cfg.OUT_DIR, ref = pref['reference'], option = tvcopt)
  # google deep variant
  elif pref['vcaller'] == 'GDV':
    gdvopt = {
      'gpu': pref['gpgpu'],
      'processor': pref['thread_num']
      }
    if pref['targeted_seq']:
      gdvopt['target'] = pref['target']
    app.runGDVCall(input = pref['input'], output = pref['output'], ref = pref['reference'], option = gdvopt)
  
  # Recalibration variant quality
  if pref['recal_var']:
    pref['mediates'].append(pref['output'] + '.vcf.gz')
    pref['mediates'].append(pref['output'] + '.vcf.gz.tbi')
    pref['output'] = os.path.join(cfg.OUT_DIR, pref['output_name'])
    app.runGATKVRecal(input = pref['input'], output = pref['output'], ref = pref['reference'], option = {'ram': pref['use_ram']})
  
  # Erase intermediate files
  if pref['remove_mediates']:
    cleanup(pref)

def rnaseq(pref, app):
  # Preprocessing fastq
  preprocess(pref, app)
  # Mapping
  if pref['mapping']:
    mapping(pref, app)
  # Count data
  app.runHTSeqCount()
  # 

  # Erase intermediate files
  if pref['remove_mediates']:
    cleanup(pref)

def chipseq(pref, app):
  # Preprocessing fastq
  preprocess(pref, app)  
  # Mapping
  if pref['mapping']:
    mapping(pref, app)
  else:
    pref['input'] = pref['inputs'][0]
  # Call peaks
  macsopt = {}
  app.runMacs2(input = pref['input'], control = None, output = '', species = '', genome = 0, 
              option = { 'bload' : True, 'lambda' : True, 'p-val' : -1, 'q-val': -1})
  # peaks => sequence list
  #app.runMoirei()
  # Find motifs

  # Erase intermediate files
  if pref['remove_mediates']:
    cleanup(pref)

def makeReadGroupInfo(pref):
  return { 
    'ID': pref['read_group_id'], 
    'SM': pref['sample_name'],
    'PL': pref['platform']
  }


def init(pref):
  try:
    pref['root'] = os.getenv('HOME')
    pref['mode'] = sys.argv[1]
    pref['cmd_only'] = False
    # Adapter
    pref['cutadapt'] = False
    pref['adapter_site'] = 'both'
    pref['adapter_seq'] = ''

    # FqQC
    pref['fqfliter'] = False
    pref['fq_qual'] = 15
    pref['fq_minlen'] = 20
    pref['fq_complex'] = 30
    
    # Mapping
    pref['mapping'] = False
    pref['mapper'] = 'BWA'
    if pref['mode'] == 'rnaseq':
       pref['mapper'] = 'STAR'

    # Target seq
    pref['targeted_seq'] = False
    pref['target'] = ''

    # Hotspot/Knownsite
    pref['use_hotspot'] = False
    pref['known_site'] = ''

    # Control data
    pref['has_control'] = False
    pref['control_input'] = ''

    # Read Group
    pref['add_rg'] = False
    pref['read_group_id'] = ''
    pref['sample_name'] = ''
    pref['platform'] = ''

    # Reference
    pref['refdir'] = '',
    pref['reference'] = ''
    pref['refname'] = ''

    # Call
    pref['vcaller'] = 'GATK'
    pref['vcparam'] = ''
    pref['vcmotif'] = ''

    # Input
    pref['iformat'] = 'fastq'
    pref['seqtype'] = 'single'
    pref['inputs'] = []
    pref['paramdir'] = '',
    
    # Output
    pref['outdir'] = '',
    pref['output_name'] = ''

    # Others
    pref['mediates'] = []
    pref['gpgpu'] = False
    pref['remove_mediates'] = False
    pref['mark_dup'] = False
    pref['detect_sr'] = False
    pref['recal_seq'] = False
    pref['recal_var'] = False
    pref['thread_num'] = 16
    pref['use_ram'] = 32
    pref['tmp'] = '',
    pref['verbose'] = False
    opts, args = getopt.getopt(sys.argv[2:], 
    "h:aqmthcgrpdv", 
    [
      "root=", "cmd_only",
      "cutadapt", "adapter_site=", "adapter_seq=",
      "fqfilter", "fq_qual=", "fq_minlen=", "fq_complex=",
      "mapping", "mapper=", 
      "targeted_seq", "target=",
      "use_hotspot", "known_site=",
      "has_control", "control_input=",
      "add_rg", "read_group_id=", "sample_name=", "platform=",
      "refdir=", "reference=", "refname=",
      "vcaller=", "vcparam=", "vcmotif=",
      "iformat=", "paired", "input=", "paramdir=", "outdir=", "output=",
      "remove_mediates", "gpgpu", "mark_dup", "detect_sr", "recal_seq", "recal_var",
      "thread_num=", "use_ram=", "verbose", "tmp="
    ])
  except getopt.GetoptError as err:
    print(err)
    sys.exit(2)
  
  for opt, arg in opts:
    if opt == '--root':
      pref['root'] = arg
    if opt == '--cmd_only':
      pref['cmd_only'] = True
    elif opt in ('-a', '--cutadapt'):
      pref['cutadapt'] = True
    elif opt == '--adapter_site':
      pref['adapter_site'] = arg
    elif opt == '--adapter_seq':
      pref['adapter_seq'] = arg
    elif opt in ('-q', '--fqfilter'):
      pref['fqfliter'] = True
    elif opt == '--fq_qual':
      pref['fq_qual'] = arg
    elif opt == '--fq_minlen':
      pref['fq_minlen'] = arg
    elif opt == '--fq_complex':
      pref['fq_complex'] = arg
    elif opt in ('-m', '--mapping'):
      pref['mapping'] = True
    elif opt == '--mapper':
      pref['mapper'] = arg
    elif opt in ('-t', '--targeted_seq'):
      pref['targeted_seq'] = True
    elif opt == '--target':
      pref['target'] = arg
    elif opt in ('-h', '--use_hotspot'):
      pref['use_hotspot'] = True
    elif opt == '--known_site':
      pref['known_site'] = arg
    elif opt in ('-c', '--has_control'):
      pref['has_control'] = True
    elif opt == '--control_input':
      pref['control_input'] = arg
    elif opt in ('-g', '--add_rg'):
      pref['add_rg'] = True
    elif opt == '--read_group_id':
      pref['read_group_id'] = arg
    elif opt == '--sample_name':
      pref['sample_name'] = arg
    elif opt == '--platform':
      pref['platform'] = arg
    elif opt == '--refdir':
      pref['refdir'] = arg
    elif opt == '--reference':
      pref['reference'] = arg
    elif opt == '--refname':
      pref['refname'] = arg
    elif opt == '--vcaller':
      pref['vcaller'] = arg
    elif opt == '--vcparam':
      pref['vcparam'] = arg
    elif opt == '--vcmotif':
      pref['vcmotif'] = arg
    elif opt == '--iformat':
      pref['iformat'] = arg
    elif opt == '--paired':
      pref['seqtype'] = 'paired'
    elif opt == '--input':
      pref['inputs'] = arg.split(':')
    elif opt == '--paramdir':
      pref['paramdir'] = arg
    elif opt == '--output':
      pref['output_name'] = arg
    elif opt == '--outdir':
      pref['outdir'] = arg
    elif opt in ('-r', '--remove_mediates'):
      pref['remove_mediates'] = True
    elif opt in ('-p', '--mark_dup'):
      pref['mark_dup'] = True
    elif opt in ('-d', '--detect_sr'):
      pref['detect_sr'] = True
    elif opt == '--recal_seq':
      pref['recal_seq'] = True
    elif opt == '--recal_var':
      pref['recal_var'] = True
    elif opt == '--thread_num':
      pref['thread_num'] = arg
    elif opt == '--use_ram':
      pref['use_ram'] = arg
    elif opt == '--gpgpu':
      pref['gpgpu'] = True
    elif opt in ('-v', '--verbose'):
      pref['verbose'] = True
    elif opt == '--tmp':
      pref['tmp'] = arg
    else:
      print('Option: "', opt, '" is not defined.')

if __name__ == "__main__":
  pref = {}
  init(pref)
  cfg = config()
  cfg.setDir({
    'ws': pref['root'],
    'app': os.path.join(pref['root'],'MyApp'),
    'ref': pref['refdir'],
    'pref': pref['paramdir'],
    'sample': os.path.join(pref['root'],'Sample'),
    'tmp': pref['tmp'],
    'test': os.path.join(pref['root'],'test'),
    'out':  pref['outdir']
  })
  cfg.makeDirs()
  app = apprun(cfg)
  pref['reference'] = os.path.join(pref['refdir'], pref['reference'])
  if pref['verbose']:
    print(pref)
  if pref['mode'] == 'varcall':
    varcall(pref, app)
  elif pref['mode'] == 'rnaseq': 
    rnaseq(pref, app)
  elif pref['mode'] == 'chipseq': 
    chipseq(pref, app)
  else:
    print(pref['mode'],'is not defined.')
