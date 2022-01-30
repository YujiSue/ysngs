import os

class config :
  def __init__(self):
    self.SOFTWARE_INFO = {
      'SRA': { 'ver' : '2.11.3' },
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
      'MACS': {'ver':'2.2.7.1' },
      'MEME':{'ver':'5.4.1' }
    }
    self.WORK_SPACE = os.environ.get("HOME")
    self.setDefault(self)

  def setWorkSpace(self, path):
    self.WORK_SPACE = path
    
  def setDefault(self):
    self.APPS_DIR = self.WORK_SPACE+'/MyApp'
    self.REFERENCE_DIR = self.WORK_SPACE+'/Reference'
    self.PREFERENCE_DIR = self.WORK_SPACE+'/Preference'
    self.SAMPLE_DIR = self.WORK_SPACE+'/Sample'
    self.TEMPORAL = self.WORK_SPACE + '/temp'
    self.TEST_DIR = self.WORK_SPACE + '/Test'
    self.OUT_DIR = self.WORK_SPACE + '/Output'

  def setDir(self, **kwargs):
    if kwargs['ws'] and os.path.exists(kwargs['ws']):
      self.WORK_SPACE = kwargs['ws']
    if kwargs['app'] and os.path.exists(kwargs['app']):
      self.APPS_DIR = kwargs['app']
    if kwargs['ref'] and os.path.exists(kwargs['ref']):
      self.REFERENCE_DIR = kwargs['ref']
    if kwargs['pref'] and os.path.exists(kwargs['pref']):
      self.PREFERENCE_DIR = kwargs['pref']
    if kwargs['sample'] and os.path.exists(kwargs['sample']):
      self.SAMPLE_DIR = kwargs['sample']
    if kwargs['tmp'] and os.path.exists(kwargs['tmp']):
      self.TEMPORAL = kwargs['tmp']
    if kwargs['test'] and os.path.exists(kwargs['test']):
      self.TEST_DIR = kwargs['test']
    if kwargs['out'] and os.path.exists(kwargs['out']):
      self.OUT_DIR = kwargs['out']

  def makeDirs(self):
    os.makedirs(self.WORK_SPACE,exist_ok=True)
    os.makedirs(self.APPS_DIR,exist_ok=True)
    os.makedirs(self.REFERENCE_DIR,exist_ok=True)
    os.makedirs(self.PREFERENCE_DIR,exist_ok=True)
    os.makedirs(self.SAMPLE_DIR,exist_ok=True)
    os.makedirs(self.TEMPORAL,exist_ok=True)
    os.makedirs(self.TEST_DIR,exist_ok=True)
    os.makedirs(self.OUT_DIR,exist_ok=True)
