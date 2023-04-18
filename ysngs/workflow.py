


import os
import json

from ysngs import common
from ysngs import analyzer

class WFNode:
  def __init__(self):
    self.nid = ''
    self.type = 'proc'
    self.label = ''

class WorkFlow:
  def __init__(self):
    self.root = WFNode()
    self.childlen = []


  #def addSubFlow(self, flow, nextto=None):
    
    




class FlowNode:
    def __init__(self, id='', type='', label=''):
        self.id = id
        self.type = type
        self.label = label
        self.children = [] if self.type=='group' else {}

def makeNode(obj):
    return FlowNode()

class Node:
    def __init__(self, nid = '', type = '', content = None):
        self.nid = nid
        self.type = type
        self.content = content
        self.input = {}
        self.next = []
        self.enabled = True

    def setOpt(self, opt):
        self.input.update(opt)

    def setContent(self, c):
        self.content = c
        if self.type == 'obj':
            for n in self.next:
                if n.type == 'proc':
                    n.input[self.nid] = self.content
    
    def exec(self, cfg, v):
        if self.type == 'proc':
            if self.enabled:
                #print('run', self.content)
                res = self.content(cfg, self.input, v)
                if res['state'] == 0:
                    if len(self.next) == 0:
                        return res
                    else:
                        for n in self.next:
                            if n.type == 'proc':
                                n.input.update(res['product'])
                            elif n.nid in res['product']:
                                n.setContent(res['product'][n.nid])
                else:
                    return res
            else:
                for n in self.next:
                    if n.type == 'proc':
                        n.input.update(self.input)
                    elif n.nid in self.input:
                        n.setContent(self.input[n.nid])
        for n in self.next:
            if n.type == 'proc':
                return n.exec(cfg, v)
                
    def linkTo(self, n):
        if self.type == 'obj':
            n.input[self.nid] = self.content
        self.next.append(n)

class WorkFlow:
    def __init__(self, cfg):
        self.cfg = cfg
        self.root = FlowNode(id='_root_', type='group')
        self.root.children.append(FlowNode(id='_init_', type='obj'))
        self.nodes = {
            '_root_': self.root,
            '_init_': self.root.children[0]
            }
        """
        if mode == 'ref':
            self.initRefTemplate()
        elif mode == 'sra':
            self.initSRATemplate()
        elif mode == 'vc':
            self.initVCTemplate()
        elif mode == 'rna-seq':
            self.initRSTemplate()
        """

    def load(self, f):
        flow = {}
        (name,ext) = os.path.splitext(f)
        if ext == 'json':
            flow = json.load(f)
        if 'nodes' in flow:
            for node in flow['nodes']:
#                if node['id'].startsWith('_'):
#                    raise 
                self.addNode(makeNode(node))
        if 'edges' in flow:
            for edge in flow['edges']:
                self.getNode(edge['from']).linkTo(self.getNode(edge['to']))
        if 'groups' in flow:
            for group in flow['groups']:
                gnode = makeNode(group)
                children = gnode.children.copy()
                gnode.children = []
                for child in children:
                    gnode.children.append(self.getNode(child))
                self.addNode(gnode)

    def run(self):
        return self.root.exec(self.cfg, self.verbose)

    def addNode(self, n):
        self.nodes[n.nid] = n
    
    def getNode(self, nid):
        return self.nodes[nid]


    #def load(p):

    def initRefTemplate(self):
        self.addNode(Node('refurl', 'obj'))
        self.addNode(Node('download', 'proc', analyzer.downloadReference))
        # self.addNode(Node('namelist', 'obj'))
        # self.addNode(Node('rename', 'proc', analyzer.renameReference))
        # self.addNode(Node('index', 'proc', analyzer.indexReference))
        self.nodes['refurl'].linkTo(self.nodes['download'])
        self.root.linkTo(self.nodes['download'])
        # self.nodes['download'].linkTo(self.nodes['rename'])

    def initSRATemplate(self):
        self.addNode(Node('sraid', 'obj'))
        self.addNode(Node('download', 'proc', analyzer.downloadFromSRA))
        self.nodes['sraid'].linkTo(self.nodes['download'])
        self.root.linkTo(self.nodes['download'])

    def initVCTemplate(self):
        self.addNode(Node('reference', 'obj'))
        self.addNode(Node('data', 'obj'))
        self.addNode(Node('target', 'obj'))

        self.addNode(Node('qc', 'proc', analyzer.runFastQC))
        self.addNode(Node('cut', 'proc', analyzer.runCutter))
        self.addNode(Node('filter', 'procm', analyzer.runFastQFilter))
        self.addNode(Node('mapping', 'proc', analyzer.runMapping))
        self.addNode(Node('markdp', 'proc', analyzer.runPicardMD))
        self.addNode(Node('varcall', 'proc', analyzer.runVC))
        
        self.addNode(Node('remove', 'proc', analyzer.removeIntermediateFiles))
        
        self.nodes['cut'].enabled = False
        self.nodes['filter'].enabled = False
        
        self.nodes['data'].linkTo(self.nodes['qc'])
        self.nodes['reference'].linkTo(self.nodes['mapping'])
        
        self.root.linkTo(self.nodes['qc'])
        self.nodes['qc'].linkTo(self.nodes['cut'])
        self.nodes['cut'].linkTo(self.nodes['filter'])
        self.nodes['filter'].linkTo(self.nodes['mapping'])
        self.nodes['mapping'].linkTo(self.nodes['markdp'])
        self.nodes['markdp'].linkTo(self.nodes['varcall'])



    # def initRSTemplate(self):

    