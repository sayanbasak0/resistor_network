#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <chrono>

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    for (auto it=v.begin(); it != v.end(); ++it) {
        if (it != v.begin()) {
            os << ", " ;
        }
        os << *it ;
    }
    os << "]";
    return os;
}

template <typename T, typename S>
std::ostream& operator<<(std::ostream& os, const std::map<T, S>& v)
{
    os << "{";
    for (auto it=v.begin(); it != v.end(); ++it) {
        if (it != v.begin()) {
            os << ", " ;
        }
        os << it->first << ": " << it->second ;
    }
    os << "}";
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& v)
{
    os << "{";
    for (auto it : v) {
        os << it;
        if (it != *v.rbegin())
            os << ", ";
    }
    os << "}";
    return os;
}

template <typename T, typename S>
std::ostream& operator<<(std::ostream& os, const std::pair<T,S>& v)
{
    os << "(";
    os << v.first << ", " 
       << v.second << ")";
    return os;
}

class chronotimer {
    protected:
    std::chrono::high_resolution_clock::time_point tstart;
    std::chrono::high_resolution_clock::time_point tstop;
    public:
    chronotimer() {
        tstart = std::chrono::high_resolution_clock::now();
        tstop = tstart;
    }
    void start() {
        tstart = std::chrono::high_resolution_clock::now();
        tstop = tstart;
    }
    void stop(const std::string& app="") {
        tstop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tstop - tstart);
        std::cout << app << " Time : " << duration.count()/1000000.0 << " seconds\n";
    }
    void lap(const std::string& app="") {
        auto tlap = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tlap - tstop);
        std::cout << app << " Time Lap : " << duration.count()/1000000.0 << " seconds\n";
        tstop = tlap;
    }
    
};

template<typename T>
std::vector<T> range(T start, T stop, T step=1 ) {
    std::vector<T> v;
    T vi = start;
    while (vi*step < stop*step) {
        v.push_back(vi);
        vi += step;
    }
    return v;
}
