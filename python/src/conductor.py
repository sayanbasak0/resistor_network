import numpy as np
from functools import reduce
from itertools import product
import networkx as nx
import matplotlib.pyplot as plt
import time
import imageio

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

class Lattice2DGrid(Network):
    def __init__(self, latticeData, 
        Rfunc,
        terminal_start='left', # 'left','right','up','down'
        terminal_end='right', # 'left','right','up','down'
        ):
        self.terminal_start = terminal_start
        self.terminal_end = terminal_end
        self.figtitle = f"{self.terminal_start}<->{self.terminal_end}"
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
        
        Network.__init__(self, staticnodes=[self.start_node,self.end_node])
        for x in range(self.Xlen):
            for y in range(self.Ylen):
                if y>0:
                    if x==0:
                        resistance_x = Rfunc(latticeData[0,y-1], 0)
                    elif x==self.Xlen-1:
                        resistance_x = Rfunc(latticeData[-1,y-1], 0)
                    else:
                        resistance_x = Rfunc(latticeData[x-1,y-1], 0)+Rfunc(latticeData[x,y-1], 0)
                    self.add_edge_resistor((x,y),(x+1,y),resistance_x)
                if x>0:
                    if y==0:
                        resistance_y = Rfunc(latticeData[x-1,0], 1)
                    elif y==self.Ylen-1:
                        resistance_y = Rfunc(latticeData[x-1,-1], 1)
                    else:
                        resistance_y = Rfunc(latticeData[x-1,y-1], 1)+Rfunc(latticeData[x-1,y], 1)
                    self.add_edge_resistor((x,y),(x,y+1),resistance_y)
        if terminal_start=='left':
            for y in range(1,self.Ylen):
                self.add_edge_resistor(self.start_node,(0,y),0)
        elif terminal_start=='right':
            for y in range(1,self.Ylen):
                self.add_edge_resistor(self.start_node,(self.Xlen,y),0)
        elif terminal_start=='down':
            for x in range(1,self.Xlen):
                self.add_edge_resistor(self.start_node,(x,0),0)
        elif terminal_start=='up':
            for x in range(1,self.Xlen):
                self.add_edge_resistor(self.start_node,(x,self.Ylen),0)
        if terminal_end=='left':
            for y in range(1,self.Ylen):
                self.add_edge_resistor(self.end_node,(0,y),0)
        elif terminal_end=='right':
            for y in range(1,self.Ylen):
                self.add_edge_resistor(self.end_node,(self.Xlen,y),0)
        elif terminal_end=='down':
            for x in range(1,self.Xlen):
                self.add_edge_resistor(self.end_node,(x,0),0)
        elif terminal_end=='up':
            for x in range(1,self.Xlen):
                self.add_edge_resistor(self.end_node,(x,self.Ylen),0)
        
        self.figsize = (10,10)
        self.update_plot = False
        self.t_pause_plot = 0.00000001
        self.giffile = ''
        self.annotate_nodes = False
        self.node_bend_rad = (0.3)
        self.annotate_resistors = False
        self.resistor_fmt = '.2f'
        self.fig = None
        self.ax = None
        self.pos = None
        self.lx,self.ly = None,None

    def BruteForceStarMesh(self, 
        figsize=(10,10), update_plot=False, t_pause_plot=0.00000001, giffile='', 
        annotate_nodes=False, node_bend_rad='0.3',
        annotate_resistors=False, resistor_fmt='.2f'):
        tStart = time.time()
        self.figsize = figsize
        self.node_bend_rad = node_bend_rad
        self.giffile = giffile
        self.t_pause_plot = t_pause_plot
        self.update_plot = update_plot
        self.annotate_nodes = annotate_nodes
        self.annotate_resistors = annotate_resistors
        self.resistor_fmt = resistor_fmt
        if self.terminal_start==self.terminal_end:
            print(f"{self.figtitle}: R = 0")
            return
        if giffile:
            self.gifWriter = imageio.get_writer(giffile, mode='I') 
        self.plot()
        for node in list(self.nodes):
            self.star_mesh(node)
            self.plot()
        self.plot()
        if giffile:
            self.gifWriter.close() 
        if len(list(self.edges))==1:
            print(f"{self.figtitle}: R = {self[self.start_node][self.end_node]['resistance']}")
        else:
            print(f"{self.figtitle}: Incomplete!")
        print(f"{self.figtitle}: Time Elapsed = {time.time()-tStart:.1f} seconds")
        self.figsize = (10,10)
        self.update_plot = True
        self.t_pause_plot = 0.00000001
        self.giffile = ''
        self.annotate_nodes = False
        self.node_bend_rad = (0.3)
        self.annotate_resistors = False
        self.resistor_fmt = '.2f'

    def BondPropagation2D(self, 
        figsize=(10,10), update_plot=False, t_pause_plot=0.00000001, giffile='', 
        annotate_nodes=False, node_bend_rad=0.3, 
        annotate_resistors=False, resistor_fmt='.2f'):
        tStart = time.time()
        self.figsize = figsize
        self.node_bend_rad = node_bend_rad
        self.giffile = giffile
        self.t_pause_plot = t_pause_plot
        self.update_plot = update_plot
        self.annotate_nodes = annotate_nodes
        self.annotate_resistors = annotate_resistors
        self.resistor_fmt = resistor_fmt
        self.fig = None
        self.ax = None
        self.pos = None
        self.lx,self.ly = None,None
        if self.terminal_start==self.terminal_end:
            print(f"{self.figtitle}: R = 0")
            return
        if self.giffile:
            self.gifWriter = imageio.get_writer(self.giffile, mode='I') 
        self.plot()
        self.remove_dangling_nodes()
        self.remove_shorted_edges()
        self.eliminate_intermediate_nodes()
        self.plot()
        if self.giffile:
            self.gifWriter.close() 
        if len(list(self.edges))==1:
            print(f"{self.figtitle}: R = {self[self.start_node][self.end_node]['resistance']}")
        else:
            print(f"{self.figtitle}: Incomplete!")
        print(f"{self.figtitle}: Time Elapsed = {time.time()-tStart:.1f} seconds")
        
    def remove_dangling_nodes(self):
        for x in range(1,self.Xlen):
            if self.degree((x,0))==1:
                self.remove_node((x,0))
                self.plot()
            if self.degree((x,self.Ylen))==1:
                self.remove_node((x,self.Ylen))
                self.plot()
        for y in range(1,self.Ylen):
            if self.degree((0,y))==1:
                self.remove_node((0,y))
                self.plot()
            if self.degree((self.Xlen,y))==1:
                self.remove_node((self.Xlen,y))
                self.plot()

    def remove_shorted_edges(self):
        for neighbor in list(self.neighbors(self.start_node)):
            self.merge_short(keepnode=self.start_node, mergenode=neighbor)
            self.plot()
        for neighbor in list(self.neighbors(self.end_node)):
            self.merge_short(keepnode=self.end_node, mergenode=neighbor)
            self.plot()

    def eliminate_intermediate_nodes(self):
        xs,ys = None,None
        for x2 in [1,self.Xlen-1]:
            for y2 in [1,self.Ylen-1]:
                if self.degree((x2,y2))==2:
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
                            nodes = self.star_mesh((x,y))
                            self.plot()
                            if len(nodes)==2:
                                self.eliminate_bond(nodes[0], nodes[1], None)
        if type(xs)==type(None) and type(ys)==type(None):
            for x in range(1,self.Xlen):
                for y in range(1,self.Ylen):
                    nodes = self.star_mesh((x,y))
                    self.plot()
                    self.eliminate_bond((x+1,y), (x,y+1), (x+1,y+1))

    def eliminate_bond(self, node1, node2, node3=None):
        if self.has_edge(node1,node2):
            nodes = list(nx.common_neighbors(self, node1, node2))
            if type(node3)!=type(None):
                if node3 in nodes:
                    nodes = [node3]
            if len(nodes)==1:
                node3 = nodes[0]
                bond = self.propagate(node1, node2, node3)
                while len(bond)==2:
                    node1,node2 = bond
                    nodes = list(nx.common_neighbors(self, node1, node2))
                    if len(nodes)==2:
                        if nodes[0]==node3:
                            node3 = nodes[1]
                        elif nodes[1]==node3:
                            node3 = nodes[0]
                        else:
                            return
                        bond = self.propagate(node1, node2, node3)
                    else:
                        return

    def propagate(self, node1, node2, node3):
        if node3 in self.staticnodes:
            return []
        newnode = tuple(map(np.mean, zip(node1, node2, node3)))
        node = self.delta_wye(node1, node2, node3, newnode)
        if self.update_plot:
            self.pos[node] = self.node_to_pos(*node)
            self.plot()
            
        retnodes = self.star_mesh(node3)
        self.plot()
        
        self.shift_node(node,node3)
        if self.update_plot:
            del(self.pos[node])
            self.plot()

        for i in range(len(retnodes)):
            if retnodes[i]==newnode:
                retnodes.pop(i)
        return retnodes
    
    def xfy(self, x, l):
        if self.node_bend_rad==0:
            return 0
        d = self.node_bend_rad*l
        R = (l**2 + d**2)/(2*d)
        y = (R**2 - (l-x)**2)**0.5 - R + d
        return y

    def node_to_pos(self, xx, yy):
        if (not self.ly) or (not self.lx):
            self.lx,self.ly = np.max(self.nodes,axis=0)
        return (xx+self.xfy(yy,self.ly/2),yy+self.xfy(xx,self.lx/2))

    def plot(self):
        if not self.update_plot:
            return
        if type(self.ax)==type(None):
            self.fig,self.ax = plt.subplots(figsize=self.figsize)
            self.fig.canvas.draw()
            plt.pause(0.01)
        if type(self.pos)==type(None):
            self.pos = {(xx,yy):self.node_to_pos(xx,yy) for xx,yy in self.nodes}
            
        self.ax.cla()
        self.ax.set_title(self.figtitle)
        nx.draw(self, pos=self.pos, node_size=1, ax=self.ax, with_labels=self.annotate_nodes)
        nx.draw_networkx_nodes(self, pos=self.pos, nodelist=[self.start_node], node_size=40, ax=self.ax, node_color='green')
        nx.draw_networkx_nodes(self, pos=self.pos, nodelist=[self.end_node], node_size=40, ax=self.ax, node_color='red')
        
        if self.annotate_resistors:
            resistances = nx.get_edge_attributes(self,'resistance')
            resistances = {edge:('{:'+self.resistor_fmt+'}').format(resistance) for edge,resistance in resistances.items()}
            nx.draw_networkx_edge_labels(self, pos=self.pos, ax=self.ax, edge_labels=resistances)
        self.fig.canvas.draw()
        plt.pause(self.t_pause_plot)
        if self.giffile:
            fig2frame(self.fig, self.gifWriter)
        

    def draw(self, get_fig_axis=False, fig=None, ax=None, figsize=(10,10), figtitle='', annotate_nodes=False, node_bend_rad=0.3, annotate_resistors=False, resistor_fmt='0.2f'):
        if type(ax)==type(None):
            fig,ax = plt.subplots(figsize=figsize)
            fig.canvas.draw()
            plt.pause(0.01)
        self.node_bend_rad = node_bend_rad
        pos = {(xx,yy):self.node_to_pos(xx,yy) for xx,yy in self.nodes}
        ax.cla()
        ax.set_title(figtitle)
        nx.draw(self, pos=pos, node_size=1, ax=ax, with_labels=annotate_nodes)
        nx.draw_networkx_nodes(self, pos=pos, nodelist=[self.start_node], node_size=40, ax=ax, node_color='green')
        nx.draw_networkx_nodes(self, pos=pos, nodelist=[self.end_node], node_size=40, ax=ax, node_color='red')
        if annotate_resistors:
            resistances = nx.get_edge_attributes(self,'resistance')
            resistances = {edge:('{:'+resistor_fmt+'}').format(resistance) for edge,resistance in resistances.items()}
            nx.draw_networkx_edge_labels(self, pos=pos, ax=ax, edge_labels=resistances)
        fig.canvas.draw()
        plt.pause(0.00000001)
        self.node_bend_rad = 0.3
        if get_fig_axis:
            return fig,ax

def fig2frame(fig, gifWriter):
    """
    Convert a matplotlib figure to a 3D numpy array with RGBA channels and appends frame to gif file writer
    Args:
        fig (matplotlib.pyplot.figure):
            a matplotlib figure
    Returns:
        None
    """
    # Get the RGBA buffer from the figure
    w,h = fig.canvas.get_width_height()
    buf = np.fromstring ( fig.canvas.tostring_argb(), dtype=np.uint8 )
    buf.shape = ( w, h, 4 )

    # canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
    buf = np.roll ( buf, 3, axis = 2 )
    gifWriter.append_data(buf)