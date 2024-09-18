import os
import importlib
import subprocess
from ysngs import installer
def checkEnv(key):
  return key in os.environ

def addEnv(key,val):
  if not(val in os.environ[key].split(':')):
    os.environ[key] += ':'+val

def setEnv(key,val):
  os.environ[key] = val

def addPath(path):
  addEnv('PATH', path)

def curlDownload(url, output=None, expand=False, showcmd=False, verbose=False):
  res = execCmd(f"curl -L {('-o ' + output) if output else '-O'} '{url}'")
  if expand:
    ext = os.path.splitext(output if output else url)[1]
    file = output if output else os.path.split(url)[1]
    if ext == '.zip':
      cmd = f"unzip {file}"
    elif ext == '.gz':
      if file.endswith('.tar.gz'):
        cmd = f"tar -zvxf {file}"
      else:
        cmd = f"gunzip {file}"
    return execCmd(cmd, showcmd=showcmd, verbose=verbose)
  else:
    return res

def gitClone(url, showcmd=False, verbose=False):
  cmd = f"git clone '{url}'"
  return execCmd(cmd, showcmd=showcmd, verbose=verbose)

def execCmd(cmd, showcmd = True, verbose = False):
  if showcmd:
    print('Run: >', cmd)
  proc = subprocess.Popen(cmd, shell=True, stderr=(subprocess.STDOUT if verbose==True else subprocess.PIPE), stdout=subprocess.PIPE, text=True)
  if verbose:
    while proc.poll() is None:
      while True: 
        line = proc.stdout.readline()
        if line:
          print(line, end='')
        else:
          break
    while True: 
      line = proc.stdout.readline()
      if line:
        print(line, end='')
      else:
        break
    return [proc.returncode==0, proc.stderr.strip() if proc.stderr else None, proc.stderr.strip() if proc.stderr else None]
  else:
    ret = proc.communicate()
    return [proc.returncode==0, ret[0].strip() if ret[0] else None, ret[1].strip() if ret[1] else None]


def execFunc(module, name, *args, **kwargs):
    try:
        func = getattr(module, name)
        return func(*args, **kwargs)
    except (ImportError, AttributeError) as e:
        print(f"Error: {e}")

def runScript(script, output=None, args=[], showcmd=False):
  cmd = "source '"+script+"'"
  if len(args):
    for arg in args:
      cmd += ' ' + str(arg)
  if output:
    cmd += ' > ' + output
  return execCmd(cmd, showcmd=showcmd, verbose=False)

def runRScript(script, output=None, args=[], showcmd=True):
  cmd = 'R --no-save --slave --vanilla'
  if len(args):
    cmd += ' --args'
    for arg in args:
      cmd += ' ' + str(arg)
  cmd += ' < ' + script
  if output and os.path.exists(output):
    cmd += ' > ' + output
  return execCmd(cmd, showcmd = True, verbose = False)

def hasMultipleObjects(v):
  return (type(v) is list)

def checkConda(name) :
  res = execCmd(f"conda list | grep {name}", showcmd = False)
  return res[0] != ''
