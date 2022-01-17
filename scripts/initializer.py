import os
import subprocess

class Setup:
  def __init__(self):
    self.WORK_SPACE = os.environ.get("HOME")
    self.APPS = WORK_SPACE+'/MyApp'
    self.TEMPORAL = WORK_SPACE+'/Downloads'
    self.REFERENCE_DIR = WORK_SPACE+'/Reference'
    self.PREFERENCE_DIR = WORK_SPACE+'/Preference'
    self.TEST_DIR = WORK_SPACE+'/Test'
    self.TEST_OUT_DIR = WORK_SPACE+'/Test/output'
  
  def setPath(options={}):
    return
  
  

