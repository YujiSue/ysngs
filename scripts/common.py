import os
import subprocess
from subprocess import PIPE

WORK_SPACE = os.environ.get("HOME")
APPS_DIR = WORK_SPACE+'/MyApp'
REFERENCE_DIR = WORK_SPACE+'/Reference'
PREFERENCE_DIR = WORK_SPACE+'/Preference'
SAMPLE_DIR = WORK_SPACE+'/Sample'
TEMPORAL = WORK_SPACE + '/temp'
TEST_DIR = WORK_SPACE + '/Test'
OUT_DIR = WORK_SPACE + '/Output'

def setWorkSpace(path):
  global WORK_SPACE
  WORK_SPACE = path

def init():
  global APPS_DIR
  global REFERENCE_DIR
  global PREFERENCE_DIR
  global SAMPLE_DIR
  global TEMPORAL
  global TEST_DIR
  global OUT_DIR
  APPS_DIR = WORK_SPACE+'/MyApp'
  REFERENCE_DIR = WORK_SPACE+'/Reference'
  PREFERENCE_DIR = WORK_SPACE+'/Preference'
  SAMPLE_DIR = WORK_SPACE+'/Sample'
  TEMPORAL = WORK_SPACE + '/temp'
  TEST_DIR = WORK_SPACE + '/Test'
  OUT_DIR = WORK_SPACE + '/Output'

def setDir(**kwargs):
  global WORK_SPACE
  global APPS_DIR
  global REFERENCE_DIR
  global PREFERENCE_DIR
  global SAMPLE_DIR
  global TEMPORAL
  global TEST_DIR
  global OUT_DIR
  if kwargs['ws'] and os.path.exists(kwargs['ws']):
    WORK_SPACE = kwargs['ws']
  if kwargs['app'] and os.path.exists(kwargs['app']):
    APPS_DIR = kwargs['app']
  if kwargs['ref'] and os.path.exists(kwargs['ref']):
    REFERENCE_DIR = kwargs['ref']
  if kwargs['pref'] and os.path.exists(kwargs['pref']):
    PREFERENCE_DIR = kwargs['pref']
  if kwargs['sample'] and os.path.exists(kwargs['sample']):
    SAMPLE_DIR = kwargs['sample']
  if kwargs['tmp'] and os.path.exists(kwargs['tmp']):
    TEMPORAL = kwargs['tmp']
  if kwargs['test'] and os.path.exists(kwargs['test']):
    TEST_DIR = kwargs['test']
  if kwargs['out'] and os.path.exists(kwargs['out']):
    OUT_DIR = kwargs['out']
  
def makeDirs():
  os.makedirs(WORK_SPACE,exist_ok=True)
  os.makedirs(APPS_DIR,exist_ok=True)
  os.makedirs(REFERENCE_DIR,exist_ok=True)
  os.makedirs(PREFERENCE_DIR,exist_ok=True)
  os.makedirs(SAMPLE_DIR,exist_ok=True)
  os.makedirs(TEMPORAL,exist_ok=True)
  os.makedirs(TEST_DIR,exist_ok=True)
  os.makedirs(OUT_DIR,exist_ok=True)

def addPath(path):
  if not(path in os.environ['PATH'].split(':')):
    os.environ['PATH'] += ':'+path

def execCmd(cmd):
  print('Run: >', cmd)
  #subprocess.Popen(cmd, shell=True)
  return

def errCheck(result):
  if result[0]:
    print(result[1])
    return True
  return False

