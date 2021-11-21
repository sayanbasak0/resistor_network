// #include "Graph2.h"
#include "covernout.h"
#include <map>
#include <iostream>
#include <algorithm>
#include <vector>
#include <set>

template <typename V, typename R, bool directed = false>
class Graph {
    protected:
    bool isDirected = directed;
    
    typename std::set<V> typedef Vend;
    typename std::map<V, Vend> typedef Vs;
    typename std::pair<V, V> typedef E;
    typename std::map<E, R> typedef Es;
    Vs node_map;
    Es edge_map;
    
    public:
    
    Graph() {
    }
    ~Graph() {
        typename Vs::iterator it;
        it = node_map.begin(); 
        while ( node_map.size() > 0 ) {
            node_map.erase( it );
            it = node_map.begin(); 
        }
        typename Es::iterator it2;
        it2 = edge_map.begin(); 
        while ( edge_map.size() > 0 ) {
            edge_map.erase( it2 );
            it2 = edge_map.begin(); 
        }
        // std::cout << "|"  ;
    }
    
    void add_node(V node) {
        if (!has_node(node)) {
            node_map[node] = Vend();
        }
    }
    void add_edge(V node1, V node2, R resistance) {
        add_node( node1 );
        add_node( node2 );
        if ( !isDirected ) {
            if (node1 > node2) {
                edge_map[ std::make_pair(node1,node2) ] = resistance;
            }
            else {
                edge_map[ std::make_pair(node2,node1) ] = resistance;
            }
            node_map[node1].insert(node2);
            node_map[node2].insert(node1);
        } else {
            edge_map[std::make_pair(node1,node2)] = resistance;
            node_map[node1].insert(node2);
        }
    }
    bool has_node(V node) {
        typename Vs::iterator it;
        it = node_map.find(node);
        if (it != node_map.end()) {
            return true;
        } else {
            if ( isDirected ) {
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
        if ( !isDirected ) {
            if ( node1 > node2 ) {
                typename Es::iterator it;
                it = edge_map.find( std::make_pair(node1, node2) );
                if ( it != edge_map.end() ) {
                    return true;
                } else {
                    return false;
                }
            } else {
                typename Es::iterator it;
                it = edge_map.find( std::make_pair(node2, node1) );
                if ( it != edge_map.end() ) {
                    return true;
                } else {
                    return false;
                }
            }
        } else {
            typename Es::iterator it;
            it = edge_map.find( std::make_pair(node1, node1) );
            if (it != edge_map.end()) {
                return true;
            } else {
                return false;
            }
        }
    }
    
    void remove_node(V node) {
        std::vector<V> neighbors = get_node_neighbors( node );
        for ( auto it : neighbors ) {
            remove_edge( node, it );
        }
        node_map.erase( node );
    }
    void remove_edge(V node1, V node2) {
        if ( !isDirected ) {
            if ( node1 > node2 ) {
                edge_map.erase( std::make_pair(node1, node2) );
            } else {
                edge_map.erase( std::make_pair(node2, node1) );
            }
            node_map[ node1 ].erase( node2 );
            node_map[ node2 ].erase( node1 );
        } else {
            edge_map.erase( std::make_pair(node1, node2) );
            node_map[ node1 ].erase( node2 );
        }
    }
    void show() {
        std::cout << "======\n" 
                  << edge_map 
                  << "\n------\n"
                  << node_map 
                  << "\n======\n" ;
    }
    R get_edge_resistance(V node1, V node2) {
        if ( !isDirected ) {
            if ( node1 > node2 ) {
                return edge_map[ std::make_pair(node1, node2) ];
            } else {
                return edge_map[ std::make_pair(node2, node1) ];
            }
        } else {
            return edge_map[ std::make_pair(node1, node2) ];
        }
    }
    std::vector<V> get_node_neighbors(V node) {
        std::vector<V> nodelist;
        for (auto it : node_map[node]) {
            nodelist.push_back(it);
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
        return edge_map.size();
    }
    void clear() {
        typename Vs::iterator it;
        it = node_map.begin(); 
        while ( node_map.size() > 0 ) {
            node_map.erase( it );
            it = node_map.begin(); 
        }
        typename Es::iterator it2;
        it2 = edge_map.begin(); 
        while ( edge_map.size() > 0 ) {
            edge_map.erase( it2 );
            it2 = edge_map.begin(); 
        }
    }
    int get_node_degree(V node) {
        if ( has_node( node ) ) {
            return node_map[node].size();
        }
        return 0;
    }
    
};

