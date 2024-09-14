import os
import json
from ysngs import common
from ysngs import installer

class AppManager:
    def __init__(self):
        current = os.path.dirname(os.path.abspath(__file__))
        self.apps = json.load(open(os.path.join(current, 'ngsapp.json')))
    
    def applist(self):
        info = {}
        keys = self.apps.keys()
        for k in keys:
            info[k] = {
                'name': self.apps[k]['name'],
                'version': self.apps[k]['ver']
            }
        return info
    
    def is_installed(self, name):
        return common.execFunc(installer, self.apps[name]['checker'])
    
    def install(self, name):
        return common.execFunc(installer, self.apps[name]['installer'], self.apps[name])
    # def update(self):
    # def uninstall(self):
    
    def version(self, name):
        return common.execFunc(installer, self.apps[name]['verchecker'])
