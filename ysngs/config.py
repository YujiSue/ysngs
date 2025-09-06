import json
import os
import platform
from ysngs import common
# Config class
class Config :
  def __init__(self):
    # Get/make the workspace path
    pl = platform.system()
    cfgdir = ''
    if pl == 'Linux':
      cfgdir = '/etc/hayami'
    elif pl == 'Darwin':
      cfgdir = '/etc/hayami'
    elif pl == 'Windows':
      cfgdir = 'C:\\ProgramData'
    if not os.path.exists(cfgdir):
      os.makedirs(cfgdir, exist_ok=True)
    cfgpath = os.path.join(cfgdir, 'ysngs.cfg.json')
    self.cfg = {}
    if os.path.exists(cfgpath):
      self.cfg = json.load(open(cfgpath))
    else:
      ws = os.getcwd()
      self.cfg = {
          'workspace': ws,
          'dirs': {
              'app': os.path.join(ws, 'MyApp'),
              'ref': os.path.join(ws, 'Reference'),
              'pref': os.path.join(ws, 'Preference'),
              'db': os.path.join(ws, 'DB'),
              'data': os.path.join(ws, 'Data'),
              'script': os.path.join(ws, 'Script'),
              'tmp': os.path.join(ws, 'Temp'),
              'log': os.path.join(ws, 'Log')
          }
      }
      json.dump(self.cfg, open(cfgpath, 'w'))
  # 
  def makeDirs(self):
    os.makedirs(self.cfg['workspace'], exist_ok=True)
    os.makedirs(self.cfg['dirs']['app'], exist_ok=True)
    os.makedirs(self.cfg['dirs']['ref'], exist_ok=True)
    os.makedirs(self.cfg['dirs']['pref'], exist_ok=True)
    os.makedirs(self.cfg['dirs']['db'], exist_ok=True)
    os.makedirs(self.cfg['dirs']['data'], exist_ok=True)
    os.makedirs(self.cfg['dirs']['script'], exist_ok=True)
    os.makedirs(self.cfg['dirs']['tmp'], exist_ok=True)
    os.makedirs(self.cfg['dirs']['log'], exist_ok=True)

  def setEnvs(self):
    common.setEnv('HYM_WS', self.cfg['workspace'])
    common.setEnv('HYM_APP', self.cfg['dirs']['app'])
    common.setEnv('HYM_REF', self.cfg['dirs']['ref'])
    common.setEnv('HYM_PREF', self.cfg['dirs']['pref'])
    common.setEnv('HYM_DB', self.cfg['dirs']['db'])
    common.setEnv('HYM_DATA', self.cfg['dirs']['data'])
    common.setEnv('HYM_SCRIPT', self.cfg['dirs']['script'])
    common.setEnv('HYM_TEMP', self.cfg['dirs']['tmp'])
    common.setEnv('HYM_LOG', self.cfg['dirs']['log'])
    