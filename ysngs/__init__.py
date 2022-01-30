import os
from .common import config
from .install import installer
from .run import apprun

def addPath(path):
  if not(path in os.environ['PATH'].split(':')):
    os.environ['PATH'] += ':'+path

