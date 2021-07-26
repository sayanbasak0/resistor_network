from operator import inv
from symbol import return_stmt, star_expr
import numpy as np
from functools import reduce
from itertools import product
import networkx as nx
import matplotlib.pyplot as plt
import time

class Network(nx.Graph):
    def __init__(self, staticnodes=None, figsize=(10,10), rad=0.3):
        nx.Graph.__init__(self)
        self.figsize = figsize
        self.fig = None
        self.ax = None
        self.pos = None
        if type(staticnodes)==list:
            self.staticnodes = staticnodes.copy()
        else:
            self.staticnodes = []
        self.rad = rad
        self.lx,self.ly = None,None
    def add_edge_resistors_from(self, edge_resistances):
        for node1,node2,resistance in edge_resistances:
            self.add_edge(node1, node2, resistance=resistance)
    def add_edge_resistor(self, node1, node2, resistance):
        self.add_edge(node1, node2, resistance=resistance)

    def update_edge(self, node1, node2, resistance):
        if node1==node2:
            return
        if self.has_edge(node1,node2):
            parallel_resistance = self[node1][node2]['resistance']
            if parallel_resistance * resistance==0:
                self[node1][node2]['resistance'] = 0
            else:
                resistance = (parallel_resistance * resistance)/(parallel_resistance + resistance)
                self[node1][node2]['resistance'] = resistance
        else:
            self.add_edge(node1,node2,resistance = resistance)
    
    def pop_edge(self, node1, node2):
        resistance = self[node1][node2]['resistance']
        self.remove_edge(node1,node2)
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
        min_resistance,min_edge = min(zip(resistances,neighbors),key=lambda x: x[0])
        if min_resistance==0:
            new_edges = list(map(lambda i:(min_edge,i[0],i[1]),zip(neighbors,resistances)))
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

    def shift_node(self,fromnode,tonode):
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
        if len(neighbors)==0:
            self.remove_node(mergenode)
            return False
        resistances = list(map(lambda x: self.pop_edge(mergenode,x),neighbors))
        new_edges = list(map(lambda i:(keepnode,i[0],i[1]),zip(neighbors,resistances)))
        list(map(lambda i: self.update_edge(i[0],i[1],i[2]), new_edges))
        self.remove_node(mergenode)
        return True
    
    def xfy(self,x,l):
        if self.rad==0:
            return 0
        d = self.rad*l
        R = (l**2 + d**2)/(2*d)
        y = (R**2 - (l-x)**2)**0.5 - R + d
        return y
    def node_to_pos(self,xx,yy):
        if not self.ly:
            self.lx,self.ly = np.max(self.nodes,axis=0)
        return (xx+self.xfy(yy,self.ly/2),yy+self.xfy(xx,self.lx/2))
    def draw(self, tViz=0.0001, annotate_resistors=True, annotate_nodes=False , rad=0.3, title=''):
        if type(self.ax)==type(None):
            self.fig,self.ax = plt.subplots(figsize=self.figsize)
        if type(self.pos)==type(None) or any([node not in self.pos.keys() for node in self.nodes]):
            try:
                self.pos = {(xx,yy):self.node_to_pos(xx,yy) for xx,yy in self.nodes}
            except:
                self.pos = nx.kamada_kawai_layout(self)
        self.ax.cla()
        self.ax.set_title(title)
        nx.draw(self, pos=self.pos, node_size=1, ax=self.ax, with_labels=annotate_nodes)
        nx.draw_networkx_nodes(self, pos=self.pos, nodelist=self.staticnodes, node_size=40, ax=self.ax)
        
        if annotate_resistors:
            resistances = nx.get_edge_attributes(self,'resistance')
            nx.draw_networkx_edge_labels(self, pos=self.pos, ax=self.ax, edge_labels=resistances)
        self.fig.canvas.draw()
        plt.pause(tViz)
    def equivalent_resistance(self):
        for node in list(self.nodes):
            self.star_mesh(node)
        return nx.get_edge_attributes(self,'resistance')

class Lattice2DGrid:
    def __init__(self, latticeData, 
        Rfunc,
        terminal_start='left', # 'left','right','up','down'
        terminal_end='right', # 'left','right','up','down'
        rad=0.3,
        plot_initial=False
        ):

        self.terminal_start = terminal_start
        self.terminal_end = terminal_end
        if terminal_start==terminal_end:
            print(f"{terminal_start}<->{terminal_end}: Shorted measurement terminals!")
            return
        
        self.Xlen,self.Ylen = latticeData.shape
        self.Xlen,self.Ylen = self.Xlen+1,self.Ylen+1
        
        if terminal_start=='left':
            self.start_node = (-1,self.Ylen/2)
        if terminal_start=='right':
            self.start_node = (self.Xlen+1,self.Ylen/2)
        if terminal_start=='down':
            self.start_node = (self.Xlen/2,-1)
        if terminal_start=='up':
            self.start_node = (self.Xlen/2,self.Ylen+1)
        if terminal_end=='left':
            self.end_node = (-1,self.Ylen/2)
        if terminal_end=='right':
            self.end_node = (self.Xlen+1,self.Ylen/2)
        if terminal_end=='down':
            self.end_node = (self.Xlen/2,-1)
        if terminal_end=='up':
            self.end_node = (self.Xlen/2,self.Ylen+1)
        
        self.ntwrk = Network(staticnodes=[self.start_node,self.end_node], rad=rad)
        for x in range(self.Xlen):
            for y in range(self.Ylen):
                if y>0:
                    if x==0:
                        resistance_x = Rfunc(latticeData[0,y-1], 0)
                    elif x==self.Xlen-1:
                        resistance_x = Rfunc(latticeData[-1,y-1], 0)
                    else:
                        resistance_x = Rfunc(latticeData[x-1,y-1], 0)+Rfunc(latticeData[x,y-1], 0)
                    self.ntwrk.add_edge_resistor((x,y),(x+1,y),resistance_x)
                if x>0:
                    if y==0:
                        resistance_y = Rfunc(latticeData[x-1,0], 1)
                    elif y==self.Ylen-1:
                        resistance_y = Rfunc(latticeData[x-1,-1], 1)
                    else:
                        resistance_y = Rfunc(latticeData[x-1,y-1], 1)+Rfunc(latticeData[x-1,y], 1)
                    self.ntwrk.add_edge_resistor((x,y),(x,y+1),resistance_y)
        if terminal_start=='left':
            for y in range(1,self.Ylen):
                self.ntwrk.add_edge_resistor(self.start_node,(0,y),0)
        elif terminal_start=='right':
            for y in range(1,self.Ylen):
                self.ntwrk.add_edge_resistor(self.start_node,(self.Xlen,y),0)
        elif terminal_start=='down':
            for x in range(1,self.Xlen):
                self.ntwrk.add_edge_resistor(self.start_node,(x,0),0)
        elif terminal_start=='up':
            for x in range(1,self.Xlen):
                self.ntwrk.add_edge_resistor(self.start_node,(x,self.Ylen),0)
        if terminal_end=='left':
            for y in range(1,self.Ylen):
                self.ntwrk.add_edge_resistor(self.end_node,(0,y),0)
        elif terminal_end=='right':
            for y in range(1,self.Ylen):
                self.ntwrk.add_edge_resistor(self.end_node,(self.Xlen,y),0)
        elif terminal_end=='down':
            for x in range(1,self.Xlen):
                self.ntwrk.add_edge_resistor(self.end_node,(x,0),0)
        elif terminal_end=='up':
            for x in range(1,self.Xlen):
                self.ntwrk.add_edge_resistor(self.end_node,(x,self.Ylen),0)
        if plot_initial:
            self.ntwrk.draw(title=f"{terminal_start}<->{terminal_end}")
    def BruteForceStarMesh(self, update_plot=False, tViz=0.0001, annotate_resistors=False, annotate_nodes=False, plot_final=False):
        tStart = time.time()
        if self.terminal_start==self.terminal_end:
            print(f"{self.terminal_start}<->{self.terminal_end}: R = 0")
            return
        for node in list(self.ntwrk.nodes):
            self.ntwrk.star_mesh(node)
            if update_plot:
                self.ntwrk.draw(title=f"{self.terminal_start}<->{self.terminal_end}", tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
        if plot_final:
            self.ntwrk.draw(title=f"{self.terminal_start}<->{self.terminal_end}")
        if len(list(self.ntwrk.edges))==1:
            print(f"{self.terminal_start}<->{self.terminal_end}: R = {self.ntwrk[self.start_node][self.end_node]['resistance']}")
        else:
            print(f"{self.terminal_start}<->{self.terminal_end}: Incomplete!")
        print(f"{self.terminal_start}<->{self.terminal_end}: Time Elapsed = {time.time()-tStart:.1f} seconds")
        
    def BondPropagation2D(self, update_plot=False, tViz=0.0001, annotate_resistors=False, annotate_nodes=False, plot_final=False):
        tStart = time.time()
        if self.terminal_start==self.terminal_end:
            print(f"{self.terminal_start}<->{self.terminal_end}: R = 0")
            return
        self.remove_dangling_nodes(update_plot=update_plot, tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
        self.remove_shorted_edges(update_plot=update_plot, tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
        self.eliminate_nodes(update_plot=update_plot, tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
        if plot_final:
            self.ntwrk.draw(title=f"{self.terminal_start}<->{self.terminal_end}")
        if len(list(self.ntwrk.edges))==1:
            print(f"{self.terminal_start}<->{self.terminal_end}: R = {self.ntwrk[self.start_node][self.end_node]['resistance']}")
        else:
            print(f"{self.terminal_start}<->{self.terminal_end}: Incomplete!")
        print(f"{self.terminal_start}<->{self.terminal_end}: Time Elapsed = {time.time()-tStart:.1f} seconds")
        
    def remove_dangling_nodes(self, update_plot=False, tViz=0.0001, annotate_resistors=False, annotate_nodes=False):
        for x in range(1,self.Xlen):
            if self.ntwrk.degree((x,0))==1:
                self.ntwrk.remove_node((x,0))
            if update_plot:
                self.ntwrk.draw(tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
            if self.ntwrk.degree((x,self.Ylen))==1:
                self.ntwrk.remove_node((x,self.Ylen))
            if update_plot:
                self.ntwrk.draw(tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
        for y in range(1,self.Ylen):
            if self.ntwrk.degree((0,y))==1:
                self.ntwrk.remove_node((0,y))
            if update_plot:
                self.ntwrk.draw(tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
            if self.ntwrk.degree((self.Xlen,y))==1:
                self.ntwrk.remove_node((self.Xlen,y))
            if update_plot:
                self.ntwrk.draw(tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
    def remove_shorted_edges(self, update_plot=False, tViz=0.0001, annotate_resistors=False, annotate_nodes=False):
        for neighbor in list(self.ntwrk.neighbors(self.start_node)):
            self.ntwrk.merge_short(keepnode=self.start_node, mergenode=neighbor)
            if update_plot:
                self.ntwrk.draw(tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
        for neighbor in list(self.ntwrk.neighbors(self.end_node)):
            self.ntwrk.merge_short(keepnode=self.end_node, mergenode=neighbor)
            if update_plot:
                self.ntwrk.draw(tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)

    def eliminate_nodes(self, update_plot=False, tViz=0.0001, annotate_resistors=False, annotate_nodes=False):
        xs,ys = None,None
        for x2 in [1,self.Xlen-1]:
            for y2 in [1,self.Ylen-1]:
                if self.ntwrk.degree((x2,y2))==2:
                    if x2==1:
                        xs = range(1,self.Xlen)
                    else:
                        xs = range(self.Xlen-1,0,-1)
                    if y2==1:
                        ys = range(1,self.Ylen)
                    else:
                        ys = range(self.Ylen-1,0,-1)
                    for x in xs:
                        for y in ys:
                            nodes = self.ntwrk.star_mesh((x,y))
                            if len(nodes)==2:
                                self.eliminate_bond(nodes[0],nodes[1], update_plot=update_plot, tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
        if type(xs)==type(None) and type(ys)==type(None):
            for x in range(1,self.Xlen):
                for y in range(1,self.Ylen):
                    nodes = self.ntwrk.star_mesh((x,y))
                    self.eliminate_bond((x+1,y), (x,y+1), (x+1,y+1), update_plot=update_plot, tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)

                
    def eliminate_bond(self, node1, node2, node3=None, update_plot=False, tViz=0.0001, annotate_resistors=False, annotate_nodes=False):
        if self.ntwrk.has_edge(node1,node2):
            nodes = list(nx.common_neighbors(self.ntwrk,node1,node2))
            if type(node3)!=type(None):
                if node3 in nodes:
                    nodes = [node3]
            if len(nodes)==1:
                node3 = nodes[0]
                bond = self.propagate(node1, node2, node3, update_plot=update_plot, tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
                while len(bond)==2:
                    nodes = list(nx.common_neighbors(self.ntwrk,bond[0],bond[1]))
                    if len(nodes)==2:
                        if nodes[0]==node3:
                            node3 = nodes[1]
                        elif nodes[1]==node3:
                            node3 = nodes[0]
                        else:
                            return
                        bond = self.propagate(bond[0],bond[1], node3, update_plot=update_plot, tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
                    else:
                        return

    def propagate(self, node1,node2,node3, update_plot=False, tViz=0.0001, annotate_resistors=False, annotate_nodes=False):
        if node3 in self.ntwrk.staticnodes:
            return []
        newnode = tuple(map(np.mean,zip(node1,node2,node3)))
        node = self.ntwrk.delta_wye(node1,node2,node3,newnode)
        if update_plot:
            self.ntwrk.pos[node] = self.ntwrk.node_to_pos(*node)
            self.ntwrk.draw(tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
        retnodes = self.ntwrk.star_mesh(node3)
        if update_plot:
            self.ntwrk.draw(tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
        self.ntwrk.shift_node(node,node3)
        if update_plot:
            self.ntwrk.draw(tViz=tViz, annotate_resistors=annotate_resistors, annotate_nodes=annotate_nodes)
        for i in range(len(retnodes)):
            if retnodes[i]==newnode:
                retnodes.pop(i)
        return retnodes
    