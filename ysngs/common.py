import os
import subprocess

def addPath(path):
  if not(path in os.environ['PATH'].split(':')):
    os.environ['PATH'] += ':'+path

def runScript(script, output = None, args = []):
  cmd = "source '"+script+"'"
  if len(args):
    for arg in args:
      cmd += ' ' + str(arg)
  if output:
    cmd += ' > ' + output
  return execCmd(cmd, showcmd=True, verbose=False)

def runRScript(script, output = None, args = [], showcmd = False):
  cmd = 'R --no-save --slave --vanilla'
  if len(args):
    cmd += ' --args'
    for arg in args:
      cmd += ' ' + str(arg)
  cmd += ' < ' + script
  if output and os.path.exists(output):
    cmd += ' > ' + output
  if showcmd:
    print('Run: >', cmd)
  proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
  return [proc.returncode, proc.stdout, proc.stderr]

def download(url, output=None):
  cmd = 'curl -L'
  if output :
    cmd += " -o '" + output + "'"
  else:
    cmd += 'O'
  cmd += " '" + url + "'"
  return execCmd(cmd, showcmd=False, verbose=False)

def execCmd(cmd, showcmd = True, verbose = False):
  if showcmd:
    print('Run: >', cmd)
  proc = subprocess.Popen(cmd, shell=True, stderr=(subprocess.STDOUT if verbose==True else subprocess.PIPE), stdout=subprocess.PIPE, text=True)
  if verbose:
    while proc.poll() is None:
      line = proc.stdout.readline().strip()
      if line:
        print(line)
    return [proc.returncode==0, proc.stderr, proc.stderr]
  else:
    ret = proc.communicate()
    return [proc.returncode==0, ret[0], ret[1]]

def hasMultipleObjects(v):
  return (type(v) is list)

class Config:
  def __init__(self):
    self.WORK_SPACE = os.environ.get("HOME")
    self.setDefault()

  def setWorkSpace(self, path):
    self.WORK_SPACE = path
    
  def setDefault(self):
    self.APPS_DIR = os.path.join(self.WORK_SPACE, 'MyApp')
    self.REFERENCE_DIR = os.path.join(self.WORK_SPACE, 'Reference')
    self.DB_DIR = os.path.join(self.WORK_SPACE, 'DB')
    self.PREFERENCE_DIR = os.path.join(self.WORK_SPACE, 'Preference')
    self.SAMPLE_DIR = os.path.join(self.WORK_SPACE, 'Sample')
    self.OUT_DIR = os.path.join(self.WORK_SPACE, 'Output')
    self.SCRIPT_DIR = os.path.join(self.WORK_SPACE, 'Script')
    self.TEMPORAL = os.path.join(self.WORK_SPACE, 'temp')
    self.TEST_DIR = os.path.join(self.WORK_SPACE, 'Test')
    self.thread = 1
    self.verbose = False
    
  def setDir(self, dirs = {}):
    if 'ws' in dirs:
      self.WORK_SPACE = dirs['ws']
    if 'app' in dirs:
      self.APPS_DIR = dirs['app']
    if 'ref' in dirs:
      self.REFERENCE_DIR = dirs['ref']
    if 'db' in dirs:
      self.DB_DIR = dirs['db']
    if 'pref' in dirs:
      self.PREFERENCE_DIR = dirs['pref']
    if 'sample' in dirs:
      self.SAMPLE_DIR = dirs['sample']
    if 'out' in dirs:
      self.OUT_DIR = dirs['out']
    if 'script' in dirs:
      self.SCRIPT_DIR = dirs['script']
    if 'tmp' in dirs:
      self.TEMPORAL = dirs['tmp']
    if 'test' in dirs:
      self.TEST_DIR = dirs['test']

  def makeDirs(self):
    os.makedirs(self.WORK_SPACE,exist_ok=True)
    os.makedirs(self.APPS_DIR,exist_ok=True)
    os.makedirs(self.REFERENCE_DIR,exist_ok=True)
    os.makedirs(self.DB_DIR,exist_ok=True)
    os.makedirs(self.PREFERENCE_DIR,exist_ok=True)
    os.makedirs(self.SAMPLE_DIR,exist_ok=True)
    os.makedirs(self.OUT_DIR,exist_ok=True)
    os.makedirs(self.SCRIPT_DIR,exist_ok=True)
    os.makedirs(self.TEMPORAL,exist_ok=True)
    os.makedirs(self.TEST_DIR,exist_ok=True)
  
  def setMaxThread(self, n):
    self.thread = n
  
  def setVerbose(self, b):
    self.verbose = b