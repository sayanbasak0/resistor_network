// #include "Resistor2dGrid.h"
#include "Network.h"
#include "covernout.h"
#include <string>
#include <vector>
#include <functional>

typedef std::vector<double> V;

template<typename R>
class Resistor2dGrid: public Network<std::vector<double>,R> {
    protected:
    std::string start_terminal;
    std::string end_terminal;
    std::vector<double> start_node;
    std::vector<double> end_node;
    double Xlen=0;
    double Ylen=0;
    public:
    using Network<std::vector<double>,R>::Network;
    using Network<std::vector<double>,R>::shift_node;
    using Network<std::vector<double>,R>::merge_short;
    using Network<std::vector<double>,R>::star_mesh;
    using Network<std::vector<double>,R>::delta_wye;
    using Network<std::vector<double>,R>::wye_delta;
    using Network<std::vector<double>,R>::add_edge_resistor;
    using Network<std::vector<double>,R>::evaluate_equivalent;
    
    using Network<std::vector<double>,R>::add_edge;
    using Network<std::vector<double>,R>::add_node;
    using Network<std::vector<double>,R>::has_edge;
    using Network<std::vector<double>,R>::has_node;
    using Network<std::vector<double>,R>::remove_edge;
    using Network<std::vector<double>,R>::remove_node;
    using Network<std::vector<double>,R>::show;
    using Network<std::vector<double>,R>::get_edge_resistance;
    using Network<std::vector<double>,R>::get_node_neighbors;
    using Network<std::vector<double>,R>::get_nodes;
    using Network<std::vector<double>,R>::no_of_nodes;
    using Network<std::vector<double>,R>::no_of_edges;
    using Network<std::vector<double>,R>::clear;
    using Network<std::vector<double>,R>::get_node_degree;
    using Network<std::vector<double>,R>::common_neighbors;
    
    void initialize(const std::string& terminal_start, 
                    const std::string& terminal_end, 
                    std::vector<std::vector<R>> latticeData,
                    std::function<R(R,int)> Rfunc
                    ) {
        start_terminal = terminal_start;
        end_terminal = terminal_end;
        Xlen = latticeData.size();
        Ylen = (*latticeData.begin()).size();
        Xlen = Xlen+1;
        Ylen = Ylen+1;
        
        if (terminal_start==terminal_end) {
            std::cout << terminal_start << "<->" << terminal_end 
                      << " : Shorted measurement terminals!" << std::endl; 
            return;
        }
        
        
        if (terminal_start=="left") {
            start_node = {-1,Ylen/2};
        } else {
            if (terminal_start=="right") {
                start_node = {Xlen+1,Ylen/2};
            } else {
                if (terminal_start=="down") {
                    start_node = {Xlen/2,-1};
                } else {
                    if (terminal_start=="up") {
                        start_node = {Xlen/2,Ylen+1};
                    }
                }
            }
        }
        if (terminal_end=="left") {
            end_node = {-1,Ylen/2};
        } else {
            if (terminal_end=="right") {
                end_node = {Xlen+1,Ylen/2};
            } else {
                if (terminal_end=="down") {
                    end_node = {Xlen/2,-1};
                } else {
                    if (terminal_end=="up") {
                        end_node = {Xlen/2,Ylen+1};
                    }
                }
            }
        }
        R resistance_x;
        R resistance_y;
        for (int x=0; x<Xlen; x++) {
            for (int y=0; y<Ylen; y++) {
                if (y>0) {
                    if (x==0) {
                        resistance_x = Rfunc(latticeData[0][y-1], 0);
                    } else {
                        if (x==Xlen-1) {
                            resistance_x = Rfunc(latticeData[x-1][y-1], 0);
                        } else {
                            resistance_x = Rfunc(latticeData[x-1][y-1], 0) + Rfunc(latticeData[x][y-1], 0);
                        }
                    }
                    add_edge_resistor({(double) x,(double) y}, {(double) (x+1),(double) y}, resistance_x);
                }
                if (x>0) {
                    if (y==0) {
                        resistance_y = Rfunc(latticeData[x-1][0], 1);
                    } else {
                        if (y==Ylen-1) {
                            resistance_y = Rfunc(latticeData[x-1][y-1], 1);
                        } else {
                            resistance_y = Rfunc(latticeData[x-1][y-1], 1) + Rfunc(latticeData[x-1][y], 1);
                        }
                    }
                    add_edge_resistor({(double) x,(double) y}, {(double) x,(double) (y+1)}, resistance_y);
                }
            }
        }
        if (terminal_start=="left") {
            for (int y=1; y<Ylen; y++) {
                add_edge_resistor(start_node, {0,(double) y}, 0);
            }
        } else {
            if (terminal_start=="right") {
                for (int y=1; y<Ylen; y++) {
                    add_edge_resistor(start_node, {Xlen,(double) y}, 0);
                }
            } else {
                if (terminal_start=="down") {
                    for (int x=1; x<Xlen; x++) {
                        add_edge_resistor(start_node, {(double) x,0}, 0);
                    }
                } else { 
                    if (terminal_start=="up") {
                        for (int x=1; x<Xlen; x++) {
                            add_edge_resistor(start_node, {(double) x,Ylen}, 0);
                        }
                    }
                }
            }
        }
        if (terminal_end=="left") {
            for (int y=1; y<Ylen; y++) {
                add_edge_resistor(end_node, {0,(double) y}, 0);
            }
        } else {
            if (terminal_end=="right") {
                for (int y=1; y<Ylen; y++) {
                    add_edge_resistor(end_node, {Xlen,(double) y}, 0);
                }
            } else {
                if (terminal_end=="down") {
                    for (int x=1; x<Xlen; x++) {
                        add_edge_resistor(end_node, {(double) x,0}, 0);
                    }
                } else { 
                    if (terminal_end=="up") {
                        for (int x=1; x<Xlen; x++) {
                            add_edge_resistor(end_node, {(double) x,Ylen}, 0);
                        }
                    }
                }
            }
        }
    }
    void BruteForceStarMesh() {
        evaluate_equivalent(start_node, end_node);
        // for (int x=0; x<Xlen+1; x++) {
        //     for (int y=0; y<Xlen+1; y++) {
        //         star_mesh({(double) x, (double) y});
        //     }
        //         std::cout << x 
        //                   << ": " << no_of_nodes()
        //                   << ": " << no_of_edges()
        //                   << std::endl;
        // }
    }
        
    void remove_dangling_nodes() {
        for (int x=1; x<Xlen; x++) {
            if ( get_node_degree( {(double) x,0} ) == 1 ) {
                remove_node( {(double) x,0} );
            }
            if ( get_node_degree( {(double) x,Ylen} ) == 1 ) {
                remove_node( {(double) x,Ylen});
            }
        }
        for (int y=1; y<Ylen; y++) {
            if ( get_node_degree( {0,(double) y} ) == 1 ) {
                remove_node( {0,(double) y} );
            }
            if ( get_node_degree( {Xlen,(double) y} ) == 1 ) {
                remove_node( {Xlen,(double) y} );
            }
        }
    }
    void remove_shorted_edges() {
        for ( auto neighbor : get_node_neighbors(start_node) ) {
            merge_short(start_node, neighbor);
        }
        for ( auto neighbor : get_node_neighbors(end_node) ) {
            merge_short(end_node, neighbor);
        }
    }
    std::vector<std::vector<double>> propagate(std::vector<double> node1, std::vector<double> node2, std::vector<double> node3) {
        std::vector<std::vector<double>> retnodes;
        if (node3==start_node || node3==end_node) {
            return retnodes;
        }
        std::vector<double> newnode({(node1[0]+node2[0]+node3[0])/3,(node1[1]+node2[1]+node3[1])/3});
        std::vector<double> node = delta_wye(node1, node2, node3, newnode);
            
        retnodes = star_mesh(node3);
        
        shift_node(node,node3);
        
        for (auto it=retnodes.begin(); it!=retnodes.end(); ++it) {
            if (*it==newnode) {
                retnodes.erase(it);
                return retnodes;
            }
        }
        return retnodes;
    }
    void eliminate_bond(std::vector<double> node1, std::vector<double> node2, std::vector<double> node3, bool automatic) {
        if ( has_edge(node1,node2) ) {
            std::vector<std::vector<double>> nodes;
            std::vector<std::vector<double>> bond;
            nodes = common_neighbors({node1, node2});
            if (!automatic) {
                for (auto node : nodes) {
                    if (node==node3) {
                        nodes = std::vector<std::vector<double>>({node3});
                    }
                }
            }
            if (nodes.size() == 1){
                node3 = nodes[0];
                bond = propagate(node1, node2, node3);
                while (bond.size()==2) {
                    node1 = bond[0];
                    node2 = bond[1];
                    nodes = common_neighbors({node1, node2});
                    if (nodes.size()==2) {
                        if (nodes[0]==node3) {
                            node3 = nodes[1];
                        } else {
                            if (nodes[1]==node3) {
                                node3 = nodes[0];
                            } else {
                                return;
                            }
                        }
                        bond = propagate(node1, node2, node3);
                    } else {
                        return;
                    }
                }
            }
        }
    }
    void eliminate_intermediate_nodes() {
        std::vector<int> xs;
        std::vector<int> ys;
        std::vector<std::vector<double>> nodes;
        bool done = false;
        for (auto x2 : std::vector<double>({1,Xlen-1}) ) {
            for (auto y2 : std::vector<double>({1,Ylen-1}) ) {
                if ( get_node_degree( {x2,y2} ) == 2 ) {
                    if (x2==1) {
                        xs = range<int>(1,Xlen);
                    } else {
                        xs = range<int>(Xlen-1,0,-1);
                    }
                    if (y2==1) {
                        ys = range<int>(1,Ylen);
                    } else {
                        ys = range<int>(Ylen-1,0,-1);
                    }
                    for (auto x : xs) {
                        for (auto y : ys) {
                            nodes = star_mesh( {(double) x,(double) y} );
                            if (nodes.size() == 2 ) {
                                eliminate_bond(nodes[0], nodes[1], nodes[1], true);
                            }
                        }
                    }
                    done = true;
                }
            }
        }
        if ( !done ) {
            for (auto x : range<int>(1,Xlen)) {
                for (auto y : range<int>(1,Ylen)) {
                    nodes = star_mesh( {(double) x,(double) y} );
                    eliminate_bond( {(double) (x+1),(double) y}, {(double) x,(double) (y+1)}, {(double) (x+1),(double) (y+1)}, false );
                }
            }
        }
    }
        
    void BondPropagation2D() {
        remove_dangling_nodes();
        remove_shorted_edges();
        eliminate_intermediate_nodes();
    }
    void info() {
        std::cout << "---------------- " << std::endl
                  << "[Start] Terminal : " << start_terminal << ", Node : " << start_node << std::endl
                  << "[End] Terminal : " << end_terminal << ", Node : " << end_node << std::endl
                  << "Grid Dimension : " << Xlen << " x " << Ylen << std::endl
                  << "No. of Nodes : " << no_of_nodes() << std::endl
                  << "No. of Edges : " << no_of_edges() << std::endl
                  << "---------------- " << std::endl;
    }
};

