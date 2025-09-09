import os
import json
from ysngs import common
from ysngs import appmanager
# WorkFlow class
class WorkFlow:
  def __init__(self, script, opts = {}):
    self.script = script
    self.apps = appmanager.AppManager()
    os.chdir(os.environ['HYM_WS'])
    if not self.apps.is_installed('cromwell'):
        self.apps.install('cromwell')

  def run(self, input):
    res = common.execCmd(f"java -jar $HYM_APP/cromwell.jar run --inputs {input} {self.script}", verbose=True)
    assert res[0], (res[2] if res[1] == None else res[1])
  
  def prepare(self, path):
    res = common.execCmd(f"java -jar $HYM_APP/womtool.jar inputs {self.script} > {path}", showcmd=False, verbose=True)
    assert res[0], (res[2] if res[1] == None else res[1])

  def check(self):
    res = common.execCmd(f"java -jar $HYM_APP/womtool.jar validate {self.script}", showcmd=False, verbose=True)
    assert res[0], (res[2] if res[1] == None else res[1])

  def graph(self, output, detailed = False):
    assert common.execCmd(f"java -jar $HYM_APP/womtool.jar {'womgraph' if detailed else 'graph'} {self.script} > {output}", showcmd=False, verbose=True)[0]
    img = output.replace('.dot', '.png')
    assert common.execCmd(f"dot -Tpng {output} > {img}", showcmd=False, verbose=True)[0]
    print(f"Workflow image was exported to '{img}'")

  def install(self, app, use_conda = False, use_pip = False):
    self.apps.install(app)
  