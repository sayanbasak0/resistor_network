from functools import reduce
from itertools import product
import networkx as nx
import matplotlib.pyplot as plt

class Network(nx.Graph):
    def __init__(self, staticnodes=None):
        nx.Graph.__init__(self)
        
        if type(staticnodes)==list:
            self.staticnodes = staticnodes.copy()
        else:
            self.staticnodes = []
        
    def add_edge_resistors_from(self, edge_resistances):
        for node1,node2,resistance in edge_resistances:
            self.add_edge(node1, node2, resistance=resistance)
    def add_edge_resistor(self, node1, node2, resistance):
        self.add_edge(node1, node2, resistance=resistance)

    def update_edge(self, node1, node2, resistance):
        if node1==node2:
            return
        if self.has_edge(node1, node2):
            parallel_resistance = self[node1][node2]['resistance']
            if parallel_resistance * resistance==0:
                self[node1][node2]['resistance'] = 0
            else:
                resistance = (parallel_resistance * resistance)/(parallel_resistance + resistance)
                self[node1][node2]['resistance'] = resistance
        else:
            self.add_edge(node1, node2, resistance = resistance)

    def pop_edge(self, node1, node2):
        resistance = self[node1][node2]['resistance']
        self.remove_edge(node1, node2)
        return resistance

    def star_mesh(self, node):
        if node in self.staticnodes:
            return []
        if not self.has_node(node):
            return []
        neighbors = list(self.neighbors(node))
        nneighbors = len(neighbors)
        if nneighbors<=1:
            self.remove_node(node)
            return neighbors
        resistances = list(map(lambda x: self.pop_edge(node,x),neighbors))
        min_resistance,min_neighbor = min(zip(resistances,neighbors),key=lambda x: x[0])
        if min_resistance==0:
            new_edges = list(map(lambda i:(min_neighbor,i[0],i[1]),zip(neighbors,resistances)))
            list(map(lambda i: self.update_edge(i[0],i[1],i[2]), new_edges))
        else:
            harmonic_sum = reduce(lambda a,b: a+b, map(lambda x:1/x,resistances))
            new_edges = list(filter(lambda i:i[0]>i[1],list(product(range(nneighbors),range(nneighbors)))))
            list(map(lambda i: self.update_edge(neighbors[i[0]],neighbors[i[1]],resistances[i[0]]*resistances[i[1]]*harmonic_sum), new_edges))
        self.remove_node(node)
        return neighbors

    def wye_delta(self, node):
        if node in self.staticnodes:
            return
        if self.has_node(node):
            if len(self.neighbors(node))==3:
                return self.star_mesh(node)

    def delta_wye(self, node1, node2, node3, newnode):
        if newnode in self.nodes:
            return
        nodes = [node1,node2,node3]
        edges = [(node2,node3),(node3,node1),(node1,node2)]
        for edge in edges:
            if not self.has_edge(*edge):
                return
        ress = list(map(lambda x: self.pop_edge(x[0],x[1]),edges))
        sum_r = reduce(lambda a,b: a+b, ress)
        if sum_r==0:
            for i in range(3):
                self.update_edge(nodes[i], newnode, resistance=0)
        for i in range(3):
            res_mult = 1.0
            for j in range(3):
                if i!=j:
                    res_mult *= ress[j]
            self.update_edge(nodes[i], newnode, resistance=res_mult/sum_r)
        return newnode

    def shift_node(self, fromnode, tonode):
        if fromnode in self.staticnodes:
            return
        if self.has_node(tonode):
            return
        if not self.has_node(fromnode):
            return
        nx.relabel_nodes(self,mapping={fromnode:tonode},copy=False)

    def merge_short(self, keepnode, mergenode):
        if mergenode in self.staticnodes:
            return False
        if not self.has_edge(keepnode,mergenode):
            return False
        if self[keepnode][mergenode]['resistance']!=0:
            return False
        neighbors = list(self.neighbors(mergenode))
        if len(neighbors)==1:
            self.remove_node(mergenode)
            return True
        resistances = list(map(lambda x: self.pop_edge(mergenode,x),neighbors))
        new_edges = list(map(lambda i:(keepnode,i[0],i[1]),zip(neighbors,resistances)))
        list(map(lambda i: self.update_edge(i[0],i[1],i[2]), new_edges))
        self.remove_node(mergenode)
        return True
    
    def draw(self, get_fig_axis=False, fig=None, ax=None, figsize=(10,10), figtitle='', annotate_nodes=False, annotate_resistors=False, resistor_fmt='0.2f'):
        if type(ax)==type(None):
            fig,ax = plt.subplots(figsize=figsize)
            fig.canvas.draw()
            plt.pause(0.01)
        pos = nx.kamada_kawai_layout(self)
        ax.cla()
        ax.set_title(figtitle)
        nx.draw(self, pos=pos, node_size=1, ax=ax, with_labels=annotate_nodes)
        nx.draw_networkx_nodes(self, pos=pos, nodelist=self.staticnodes, node_size=40, ax=ax, node_color='orange')
        
        if annotate_resistors:
            resistances = nx.get_edge_attributes(self,'resistance')
            resistances = {edge:('{:'+resistor_fmt+'}').format(resistance) for edge,resistance in resistances.items()}
            nx.draw_networkx_edge_labels(self, pos=pos, ax=ax, edge_labels=resistances)
        fig.canvas.draw()
        plt.pause(0.00000001)
        if get_fig_axis:
            return fig,ax
        
    def evaluate_equivalent(self):
        for node in list(self.nodes):
            self.star_mesh(node)
        return nx.get_edge_attributes(self,'resistance')
