import os

def addPath(path):
  if not(path in os.environ['PATH'].split(':')):
    os.environ['PATH'] += ':'+path

class config :
  def __init__(self):
    self.SOFTWARE_INFO = {
      'SRA': { 'ver' : '2.11.3' },
      'fastp': { 'ver' : '0.23.2' },
      'fastQC': {'ver' : '0.11.9' },
      'BWA': { 'ver':'0.7.17' },
      'TVC': { 'ver':'5.12.1' },
      'SAMTools': {'ver':'1.14' },
      'BCFTools': {'ver':'1.14' },
      'Picard': {'ver':'2.26.9' },
      'GATK': {'ver':'4.2.4.0' },
      'GDV' : {'ver':'1.3.0'},
      'Bowtie2': {'ver':'2.4.4' },
      'STAR': {'ver':'2.7.9a' },
      'Cuff': {'ver':'2.2.1' },
      'BiocManager': {'ver':'3.14' },
      'MACS': {'ver':'2.2.7.1' },
      'MEME':{'ver':'5.4.1' }
    }
    self.WORK_SPACE = os.environ.get("HOME")
    self.setDefault()

  def setWorkSpace(self, path):
    self.WORK_SPACE = path
    
  def setDefault(self):
    self.APPS_DIR = os.path.join(self.WORK_SPACE, 'MyApp')
    self.REFERENCE_DIR = os.path.join(self.WORK_SPACE, 'Reference')
    self.PREFERENCE_DIR = os.path.join(self.WORK_SPACE, 'Preference')
    self.SAMPLE_DIR = os.path.join(self.WORK_SPACE, 'Sample')
    self.SCRIPT_DIR = os.path.join(self.WORK_SPACE, 'Script')
    self.TEMPORAL = os.path.join(self.WORK_SPACE, 'temp')
    self.TEST_DIR = os.path.join(self.WORK_SPACE, 'Test')
    self.OUT_DIR = os.path.join(self.WORK_SPACE, 'Output')

  def setDir(self, dirs = {}):
    if 'ws' in dirs:
      self.WORK_SPACE = dirs['ws']
    if 'app' in dirs:
      self.APPS_DIR = dirs['app']
    if 'ref' in dirs:
      self.REFERENCE_DIR = dirs['ref']
    if 'pref' in dirs:
      self.PREFERENCE_DIR = dirs['pref']
    if 'sample' in dirs:
      self.SAMPLE_DIR = dirs['sample']
    if 'script' in dirs:
      self.SCRIPT_DIR = dirs['script']
    if 'tmp' in dirs:
      self.TEMPORAL = dirs['tmp']
    if 'test' in dirs:
      self.TEST_DIR = dirs['test']
    if 'out' in dirs:
      self.OUT_DIR = dirs['out']

  def makeDirs(self):
    os.makedirs(self.WORK_SPACE,exist_ok=True)
    os.makedirs(self.APPS_DIR,exist_ok=True)
    os.makedirs(self.REFERENCE_DIR,exist_ok=True)
    os.makedirs(self.PREFERENCE_DIR,exist_ok=True)
    os.makedirs(self.SAMPLE_DIR,exist_ok=True)
    os.makedirs(self.SCRIPT_DIR,exist_ok=True)
    os.makedirs(self.TEMPORAL,exist_ok=True)
    os.makedirs(self.TEST_DIR,exist_ok=True)
    os.makedirs(self.OUT_DIR,exist_ok=True)
