// #include "Network.h"
#include "Graph.h"
#include "covernout.h"
#include <cmath>
#include <float.h>
#include <map>
#include <iostream>
#include <algorithm>
#include <vector>

template<typename V, typename R>
class Network: public Graph<V, R> {
    protected:
    // typename std::map<V, R> typedef Vend;
    // typename std::map<V, Vend> typedef Vs;
    // // typename std::set<V> typedef Vend;
    // // typename std::map<V, Vend> typedef Vs;
    // // typename std::pair<V, V> typedef E;
    // // typename std::map<E, R> typedef Es;
    std::vector<V> static_nodes;
    public:
    using Graph<V,R>::Graph;
    using Graph<V,R>::add_node;
    using Graph<V,R>::add_edge;
    using Graph<V,R>::has_node;
    using Graph<V,R>::has_edge;
    using Graph<V,R>::remove_node;
    using Graph<V,R>::remove_edge;
    using Graph<V,R>::show;
    using Graph<V,R>::get_edge_resistance;
    using Graph<V,R>::get_node_neighbors;
    using Graph<V,R>::get_nodes;
    using Graph<V,R>::no_of_nodes;
    using Graph<V,R>::no_of_edges;
    using Graph<V,R>::clear;
    using Graph<V,R>::get_node_degree;
    using Graph<V,R>::common_neighbors;
    
    void add_edge_resistor(V node1, V node2, R resistance) {
        add_edge(node1, node2, resistance);
    }
    
    void update_edge(V node1, V node2, R resistance) {
        if (node1==node2) {
            return;
        }
        if (has_edge(node1, node2)) {
            /* typename */ R old_resistance = get_edge_resistance(node1, node2);
            if (old_resistance==0 || resistance==0) {
                add_edge(node1, node2, 0);
            } else if (std::isinf(old_resistance) && std::isinf(resistance)) {
                remove_edge(node1, node2);
            } else if (std::isinf(old_resistance) && std::isfinite(resistance)) {
                add_edge(node1, node2, resistance);
            } else if (std::isfinite(old_resistance) && std::isinf(resistance)) {
                // add_edge(node1, node2, old_resistance);
            } else if (std::isfinite(old_resistance) && std::isfinite(resistance)) {
                R new_resistance;
                if (old_resistance>resistance) {
                    R ratio = resistance/old_resistance;
                    new_resistance = old_resistance * (ratio/(1.0+ratio));
                } else if (resistance>old_resistance) {
                    R ratio = old_resistance/resistance;
                    new_resistance = resistance * (ratio/(1.0+ratio));
                } else {
                    new_resistance = resistance/2.0;
                }
                add_edge(node1, node2, new_resistance);
            }
        } else if (std::isfinite(resistance)) {
            add_edge(node1, node2, resistance);
        }
    }

    R pop_edge(V node1, V node2) {
        if (has_edge(node1, node2)) {
            /* typename */ R resistance = get_edge_resistance(node1, node2);
            remove_edge(node1, node2);
            return resistance;
        } else {
            return -1;
        }
    }

    std::vector<V> star_mesh(V node) {
        // if node in staticnodes:
        //     return []
        typename std::vector<V> neighbors;
        if (!has_node(node)) {
            return neighbors;
        }
        neighbors = get_node_neighbors(node);
        
        if (neighbors.size()<=1) {
            remove_node(node);
            return neighbors;
        }
        std::vector<R> resistances;
        int short_index = -1;
        for (auto it: range<int>(0,neighbors.size())) {
            R resistance = get_edge_resistance(node,neighbors[it]);
            resistances.push_back(resistance);
            if (resistance==0 && short_index==-1) {
                short_index = it;
            }
        }
        if (short_index>-1) {
            for (auto it: range<int>(0,neighbors.size())) {
                update_edge(neighbors[short_index],neighbors[it],resistances[it]);
            }
        } else {
            for (auto i1: range<int>(0,neighbors.size())) {
                R resistance = 0;
                R r1 = resistances[i1];
                for (auto i2: range<int>(0,neighbors.size())) {
                    R r2 = resistances[i2];
                    if (i1 > i2) {
                        R rmin = std::min(r1,r2);
                        R rmax = std::max(r1,r2);
                        R resistance = rmin+rmax;
                        if (std::isfinite(resistance)) {
                            int i3 = 0;
                            while ( i3 < neighbors.size() ) {
                                if (i3!=i1 && i3!=i2) {
                                    R r3 = resistances[i3];
                                    resistance += rmax*(rmin/r3);
                                    if (std::isinf(resistance)) {
                                        i3 = neighbors.size()+1;
                                    }
                                }
                                i3 += 1;
                            }
                            if (i3 == neighbors.size()) {
                                update_edge(neighbors[i1], neighbors[i2], resistance);
                            }
                        }
                    }
                }
            }
        }
        remove_node(node);
        return neighbors;
    }
    
    std::vector<V> wye_delta(V node) {
        // if node in self.staticnodes:
        //     return
        if (has_node(node)) {
            if (get_node_neighbors(node).size()==3) {
                return star_mesh(node);
            }
        }
        return std::vector<V>();
    }

    V delta_wye(V node1, V node2, V node3, V newnode) {
        if (has_node(newnode)) {
            return newnode;
        }
        if (!has_edge(node1, node2)) {
            return newnode;
        } else {
            if (!has_edge(node2, node3)) {
                return newnode;
            } else {
                if (!has_edge(node3, node1)) {
                    return newnode;
                }
            }
        }
        
        R ress1 = pop_edge(node2, node3);
        R ress2 = pop_edge(node3, node1);
        R ress3 = pop_edge(node1, node2);
        /* typename */ R sum_r = ress1 + ress2 + ress3;
        if (sum_r==0) {
            update_edge(node1, newnode, 0);
            update_edge(node2, newnode, 0);
            update_edge(node3, newnode, 0);
        } else {
            update_edge(node1, newnode, ress2*ress3/sum_r);
            update_edge(node2, newnode, ress3*ress1/sum_r);
            update_edge(node3, newnode, ress1*ress2/sum_r);
        }
        return newnode;
    }
    
    bool shift_node(V fromnode, V tonode) {
        // if fromnode in self.staticnodes:
        //     return
        
        if (has_node(tonode)) {
            return false;
        }
        if (!has_node(fromnode)) {
            return false;
        }
        typename std::vector<V> neighbors;
        neighbors = get_node_neighbors(fromnode);
        /* typename */ R resistance;
        for (auto it : neighbors) {
            resistance = pop_edge(fromnode, it);
            add_edge(tonode, it, resistance);
        }
        remove_node(fromnode);
        return true;
    }

    bool merge_short(V keepnode, V mergenode) {
        // if mergenode in self.staticnodes:
        //     return False
        if (!has_edge(keepnode,mergenode)) {
            return false;
        }
        if (get_edge_resistance(keepnode, mergenode)!=0) {
            return false;
        }
        typename std::vector<V> neighbors;
        neighbors = get_node_neighbors(mergenode);
        if (neighbors.size()==1) {
            remove_node(mergenode);
            return true;
        }
        /* typename */ R resistance;
        for ( auto it : neighbors ) {
            resistance = pop_edge(mergenode, it);
            update_edge(keepnode, it, resistance);
        }

        remove_node(mergenode);
        return true;
    }

    void evaluate_equivalent(V node1, V node2) {
        typename std::vector<V> nodes = get_nodes();
        for (auto it : nodes) {
            if (((it)!=node1) && ((it)!=node2)) {
                // std::cout << " " << it << "\n";
                // std::flush(std::cout);
                star_mesh(it);
            } 
        }
    }
    
};

