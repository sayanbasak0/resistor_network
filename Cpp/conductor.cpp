#include "Resistor2dGrid.h"
#include "covernout.h"
#include <vector>
#include <map>
#include <cmath>

#define PI 3.141592653589793238462643383279502884L


double Risotropic(double val, int dir) {
    return (cos(PI*val)+3)/2;
}
double Rnematic(double val, int dir) {
    return ((cos(PI*val)+3)/2)*dir + ((-cos(PI*val)+3)/2)*(1-dir) ;
}
int main() {
    // Graph<std::vector<double>,double> ntwrk;
    // ntwrk.show();
    // ntwrk.add_edge({1,1}, {2}, 1.1);
    // ntwrk.show();
    // ntwrk.add_edge({1,1}, {3,1}, 1.2);
    // ntwrk.show();
    // ntwrk.add_edge({3,1}, {4,1}, 1.3);
    // ntwrk.show();
    // ntwrk.add_node({2,2});
    // ntwrk.show();
    // ntwrk.remove_node({2,2});
    // ntwrk.show();
    // ntwrk.remove_edge({1,1}, {2,2});
    // ntwrk.show();
    // ntwrk.remove_edge({1,1}, {3,1});
    // ntwrk.show();
    // ntwrk.remove_edge({3,1}, {4,1});
    // ntwrk.show();
    // ntwrk.remove_edge({3,1}, {4,1});
    // ntwrk.show();
    // ntwrk.get_nodes();
    // ntwrk.show();
    
    // Network<std::vector<double>,double> ntwrk2;
    // ntwrk2.show();
    // ntwrk2.add_edge({1,1}, {2}, 1.1);
    // ntwrk2.show();
    // ntwrk2.add_edge({1,1}, {3,1}, 1.2);
    // ntwrk2.show();
    // ntwrk2.add_edge({3,1}, {4,1}, 1.3);
    // ntwrk2.show();
    // ntwrk2.evaluate_equivalent({3,1}, {4,1});
    // ntwrk2.show();

    chronotimer timer;
    std::vector<std::vector<double>> lat;
    for (int i=0; i<10; i++ ) {
        std::vector<double> lati;
        for (int j=0; j<10; j++ ) {
            lati.push_back(1);
        }
        lat.push_back(lati);
    }
    Resistor2dGrid<double> ntwrk3;
    
    timer.start();

    ntwrk3.initialize("left", "right", lat, Risotropic);
    ntwrk3.info();
    ntwrk3.BruteForceStarMesh();
    ntwrk3.show();
    ntwrk3.clear();
    timer.lap("[Brute-Iso-LR]");
    
    ntwrk3.initialize("left", "right", lat, Risotropic);
    ntwrk3.info();
    ntwrk3.BondPropagation2D();
    ntwrk3.show();
    ntwrk3.clear();
    timer.lap("[BondProp-Iso-LR]");
    
    ntwrk3.initialize("left", "right", lat, Rnematic);
    ntwrk3.info();
    ntwrk3.BruteForceStarMesh();
    ntwrk3.show();
    ntwrk3.clear();
    timer.lap("[Brute-Nem-LR]");
    
    ntwrk3.initialize("left", "right", lat, Rnematic);
    ntwrk3.info();
    ntwrk3.BondPropagation2D();
    ntwrk3.show();
    ntwrk3.clear();
    timer.lap("[BondProp-Nem-LR]");
    
    ntwrk3.initialize("up", "down", lat, Risotropic);
    ntwrk3.info();
    ntwrk3.BruteForceStarMesh();
    ntwrk3.show();
    ntwrk3.clear();
    timer.lap("[Brute-Iso-UD]");
    
    ntwrk3.initialize("up", "down", lat, Risotropic);
    ntwrk3.info();
    ntwrk3.BondPropagation2D();
    ntwrk3.show();
    ntwrk3.clear();
    timer.lap("[BondProp-Iso-UD]");
    
    ntwrk3.initialize("up", "down", lat, Rnematic);
    ntwrk3.info();
    ntwrk3.BruteForceStarMesh();
    ntwrk3.show();
    ntwrk3.clear();
    timer.lap("[Brute-Nem-UD]");
    
    ntwrk3.initialize("up", "down", lat, Rnematic);
    ntwrk3.info();
    ntwrk3.BondPropagation2D();
    ntwrk3.show();
    ntwrk3.clear();
    timer.lap("[BondProp-Nem-UD]");
    
    // for (auto blah : range<int>(100,1000,100)) {
    //     for (auto blah2 : range<int>(blah,1000,100)) {
    //         lat.clear();
    //         for (int i=0; i<blah; i++ ) {
    //             std::vector<double> lati;
    //             for (int j=0; j<blah2; j++ ) {
    //                 lati.push_back(1);
    //             }
    //             lat.push_back(lati);
    //         }
    //         ntwrk3.initialize("up", "down", lat, Risotropic);
    //         ntwrk3.info();
    //         ntwrk3.BondPropagation2D();
    //         ntwrk3.show();
    //         ntwrk3.clear();
    //         timer.lap("[BondProp-Iso-UD"+std::to_string(blah)+"x"+std::to_string(blah2)+"]");
    //     }
    // }
    timer.stop();
    
    return 0;
}