// #include "Graph.h"
#include "covernout.h"
#include <map>
#include <iostream>
#include <algorithm>
#include <vector>

template <typename V, typename R, bool directed = false>
class Graph {
    protected:
    bool isDirected = directed;

    typename std::map<V, R> typedef Vend;
    typename std::map<V, Vend> typedef Vs;
    Vs node_map;
    
    public:
    // Graph();
    // ~Graph();
    // void add_node(V node);
    // void add_edge(V node1, V node2, R resistance);
    // bool has_node(V node);
    // bool has_edge(V node1, V node2);
    // void remove_node(V node);
    // void remove_edge(V node1, V node2);
    // void show();
    // R get_edge_resistance(V node1, V node2);
    // Vend get_node_neighbors(V node);
    // Vs get_nodes();
    
    Graph() {
    }
    ~Graph() {
        typename Vs::iterator it;
        it = node_map.begin(); 
        while (node_map.size()>0) {
            node_map.erase(it);
            it = node_map.begin(); 
        }
        // std::cout << "|"  ;
    }
    
    void add_node(V node) {
        if (!has_node(node)) {
            node_map[node] = Vend();
        }
    }
    void add_edge(V node1, V node2, R resistance) {
        typename Vs::iterator it;
        typename Vend::iterator it2;
        it = node_map.find(node1);
        if (it != node_map.end()) {
            it2 = (it->second).find(node2);
            if (it2 != (it->second).end()) {
                (it2->second) = resistance;
            } else {
                (it->second)[node2] = resistance;
            }
        } else {
            node_map[node1] = Vend();
            node_map[node1][node2] = resistance;
        }
        if (!isDirected) {
            it = node_map.find(node2);
            if (it != node_map.end()) {
                it2 = (it->second).find(node1);
                if (it2 != (it->second).end()) {
                    (it2->second) = resistance;
                } else {
                    (it->second)[node1] = resistance;
                }
            } else {
                node_map[node2] = Vend();
                node_map[node2][node1] = resistance;
            }
        }
    }
    bool has_node(V node) {
        typename Vs::iterator it;
        it = node_map.find(node);
        if (it != node_map.end()) {
            return true;
        } else {
            if (isDirected) {
                typename Vend::iterator it2;
                for (it=node_map.begin(); it!=node_map.end(); ++it) {
                    it2 = (it->second).find(node);
                    if (it2 != (it->second).end()) {
                        return true;
                    }
                }
                return false;
            } else {
                return false;
            }
        }
    }
    bool has_edge(V node1, V node2) {
        typename Vs::iterator it;
        it = node_map.find(node1);
        if (it == node_map.end()) {
            return false;
        } else {
            typename Vend::iterator it2;
            it2 = (it->second).find(node2);
            if (it2 == (it->second).end()){
                return false;
            } else {
                return true;
            }
        }
    }
    
    void remove_node(V node) {
        typename Vs::iterator it;
        it = node_map.find(node);
        if (it != node_map.end()) {
            typename std::vector<V> neighbors = get_node_neighbors(node);
            for (auto it2 : neighbors) {
                remove_edge( node, it2 );
            }
            node_map.erase(node);
        }
    }
    void remove_edge(V node1, V node2) {
        node_map[node1].erase(node2);
        if (!isDirected) {
            node_map[node2].erase(node1);
        }
    }
    void show() {
        std::cout << "======\n" 
                  << node_map 
                  << "\n======\n" ;
    }
    R get_edge_resistance(V node1, V node2) {
        return node_map[node1][node2];
    }
    std::vector<V> get_node_neighbors(V node) {
        std::vector<V> nodelist;
        for (auto it : node_map[node]) {
            nodelist.push_back(it.first);
        }
        return nodelist;
    }
    std::vector<V> get_nodes() {
        std::vector<V> nodelist;
        for (auto it : node_map) {
            nodelist.push_back(it.first);
        }
        return nodelist;
    }
    int no_of_nodes() {
        return node_map.size();
    }
    int no_of_edges() {
        int count = 0;
        for (auto it:node_map) {
            count += (it.second).size();
        }
        if ( !isDirected ) {
            return count/2;
        } else {
            return count;
        }

    }
    void clear() {
        typename Vs::iterator it;
        it = node_map.begin(); 
        while (node_map.size()>0) {
            node_map.erase(it);
            it = node_map.begin(); 
        }
    }
    int get_node_degree(V node) {
        if ( has_node( node ) ) {
            return node_map[node].size();
        }
        return 0;
    }
    std::vector<V> common_neighbors(std::vector<V> nodes) {
        std::vector<V> neighbors0;
        if (nodes.size()==0) {
            return neighbors0;
        }
        neighbors0 = get_node_neighbors(nodes[0]);
        if (nodes.size()==1) {
            return neighbors0;
        }
        std::vector<V> neighbors;
        bool inall = true;
        for (auto neigh : neighbors0) {
            inall = true;
            for (int i=1; i<nodes.size(); i++){
                if ( !has_edge(neigh, nodes[i]) ) {
                    inall = false;
                }
            }
            if (inall) {
                neighbors.push_back(neigh);
            }
        }
        return neighbors;
    }
};

