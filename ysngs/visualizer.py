import os
from graphviz import Digraph

class Visualizer:
    def __init__(self, node, prop = {}):
        self.graph = None
        self.root = node
        self.prop = prop
        self.initGraph(node, prop)
        
    def initGraph(self):
        self.graph = Digraph(format=self.prop['format'])
        self.makeNode(self.graph, self.root)

    def makeNode(self, graph, parent):
        for node in parent.children:
            if node.type == 'GROUP':
                self.makeSubGraph(graph, node)
            elif node.type == 'OBJ':
                self.makeObjNode(graph, node)
                self.makeEdge(parent, node)
            elif node.type == 'EXE':
                self.makeExeNode(graph, node)
                self.makeEdge(parent, node)
            else:
                continue
            
    def makeObjNode(self, graph, node):
        graph.node(node.label)
    
    def makeExeNode(self, graph, node):
        graph.node(node.label)
    
    def makeEdge(self, graph, node1, node2):
        graph.edge(node1.label, node2.label)

    def makeSubGraph(self, graph, parent):
        with graph.subgraph(name=parent.label) as sub:
            for node in parent.children:
                self.makeNode(sub, node)
        
    def export(self, path):
        self.graph.render(path)
        
    
        

        


