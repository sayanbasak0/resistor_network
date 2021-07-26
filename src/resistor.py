from operator import inv
import numpy as np
from functools import reduce
from itertools import product
import networkx as nx

class Network:
    def __init__(self):
        self.nxG = nx.Graph()
        self.edges = {}
        self.nodes = {}
    def get_edge(self,node1,node2):
        return tuple([node1,node2]) if node1>node2 else tuple([node2,node1])
    def add_edge(self, node1, node2, resistance):
        if node1==node2:
            return
        edge = self.get_edge(node1,node2)
        if edge in self.edges:
            if self.edges[edge] * resistance==0:
                self.edges[edge] = 0
                self.nxG[node1][node2]['weight'] = "0"
            else:
                resistance = (self.edges[edge] * resistance)/(self.edges[edge] + resistance)
                self.edges[edge] = resistance
                self.nxG[node1][node2]['weight'] = f"{resistance:.3f}"
        else:
            if node1 not in self.nodes:
                self.nodes[node1] = {}
            self.nodes[node1][node2] = 1
            if node2 not in self.nodes:
                self.nodes[node2] = {}
            self.nodes[node2][node1] = 1
            self.edges[edge] = resistance
            self.nxG.add_edge(node1,node2,weight=f"{resistance:.3f}")
    
    def pop_edge(self, node1, node2):
        edge = self.get_edge(node1,node2)
        resistance = self.edges[edge]
        del(self.edges[edge])
        del(self.nodes[node1][node2])
        del(self.nodes[node2][node1])
        self.nxG.remove_edge(node1,node2)
        return resistance
        
    def star_mesh(self, node):
        if node not in self.nodes:
            return
        if len(self.nodes[node])==0:
            del(self.nodes[node])
            return
        ends = list(self.nodes[node].keys())
        nends = len(ends)
        ress = list(map(lambda x: self.pop_edge(node,x),ends))
        minres,minend = min(zip(ress,ends),key=lambda x: x[0])
        if minres==0:
            new_edges = list(map(lambda i:(minend,i[0],i[1]),zip(ends,ress)))
            list(map(lambda i: self.add_edge(i[0],i[1],i[2]), new_edges))
        else:
            inv_prod = reduce(lambda a,b: a+b, map(lambda x:1/x,ress))
            new_edges = list(filter(lambda i:i[0]>i[1],list(product(range(nends),range(nends)))))
            list(map(lambda i: self.add_edge(ends[i[0]],ends[i[1]],ress[i[0]]*ress[i[1]]*inv_prod), new_edges))
        del(self.nodes[node])
        self.nxG.remove_node(node)
    
    def wye_delta(self, node):
        if len(self.nodes[node].values())==3:
            self.nxG.star_mesh(node)
        
    def delta_wye(self, node1, node2, node3, newnode=None):
        if (node1 not in self.nodes) or (node2 not in self.nodes) or (node3 not in self.nodes):
            return
        if type(newnode)==type(None):
            newnode = ((node1[0]+node2[0]+node3[0])/3,(node1[1]+node2[1]+node3[1])/3)
        nodes = [node1,node2,node3]
        edges = [(node2,node3),(node3,node1),(node1,node2)]
        ress = list(map(lambda x: self.pop_edge(x[0],x[1]),edges))
        inv_sum = 1/reduce(lambda a,b: a+b, ress)
        for i in range(3):
            res_mult = 1.0
            for j in range(3):
                if i!=j:
                    res_mult *= ress[j]
            self.add_edge(nodes[i], newnode, resistance=res_mult*inv_sum)
        return newnode

    def shift_node(self,node,newnode):
        if newnode in self.nodes:
            return
        if node not in self.nodes:
            return
        if len(self.nodes[node])==0:
            del(self.nodes[node])
            return
        ends = list(self.nodes[node].keys())
        ress = list(map(lambda x: self.pop_edge(node,x),ends))
        new_edges = list(map(lambda i:(newnode,i[0],i[1]),zip(ends,ress)))
        list(map(lambda i: self.add_edge(i[0],i[1],i[2]), new_edges))
        del(self.nodes[node])
        self.nxG.remove_node(node)
    def merge_short(self, edge, node):
        if edge not in self.edges:
            return False
        if self.edges[edge]!=0:
            return False
        if len(self.nodes[node])==0:
            del(self.nodes[node])
            return False
        if node not in edge:
            for nodei in edge:
                ends1 = list(self.nodes[nodei].keys())
                ress1 = list(map(lambda x: self.pop_edge(nodei,x),ends1))
                new_edges1 = list(map(lambda i:(node,i[0],i[1]),zip(ends1,ress1)))
                list(map(lambda i: self.add_edge(i[0],i[1],i[2]), new_edges1))
                del(self.nodes[nodei])
                self.nxG.remove_node(nodei)
            return True
        else:
            nodei = edge[1-edge.index(node)]
            endsi = list(self.nodes[nodei].keys())
            ressi = list(map(lambda x: self.pop_edge(nodei,x),endsi))
            new_edgesi = list(map(lambda i:(node,i[0],i[1]),zip(endsi,ressi)))
            list(map(lambda i: self.add_edge(i[0],i[1],i[2]), new_edgesi))
            del(self.nodes[nodei])
            self.nxG.remove_node(nodei)
            return True
            