#ifndef __GRAPH_HPP__
#define __GRAPH_HPP__

// class which is basically a sparse matrix of node to node connectivity using std::map
#include <map>
#include <array>

typedef short int id;

class Graph {
 public:
    Graph(){}; // private constructor

    // add edge
    void add_edge(int v1, int v2){
        graph[{v1, v2}] = true;
    }

    // remove edge
    void remove_edge(int v1, int v2){
        graph.erase({v1, v2});
    }

    // check if edge exists
    bool operator()(int v1, int v2) {
        if (graph.find({v1, v2}) != graph.end()){
            return graph.at({v1, v2});
        } else {
            return false;
        }
    }

    void operator()(int v1, int v2, bool val){
        graph[{v1, v2}] = val;
    }

    void Print(){
        for (auto it = graph.begin(); it != graph.end(); ++it){
            std::cout << it->first[0] << " " << it->first[1] << " : " << it->second << std::endl;
        }
    }



 private:
    std::map<std::array<int,2>,bool> graph;
};

class GraphCRS {
 public:
    GraphCRS(){}; // private constructor

    id operator()(int v1, int v2){
        if (graph.find(v1) != graph.end()){
            if (graph.at(v1).find(v2) != graph.at(v1).end()){
                return graph.at(v1).at(v2);
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    void operator()(int v1, int v2, id val){
        graph[v1][v2] = val;
    }

    void Print(){
        for (auto it = graph.begin(); it != graph.end(); ++it){
            for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2){
                std::cout << it->first << " " << it2->first << " : " << it2->second << std::endl;
            }
        }
    }

    std::vector<std::array<int,2>> operator()(int v){
        std::vector<std::array<int,2>> neighbors;
        for (auto it = graph[v].begin(); it != graph[v].end(); ++it){
            if (it->second >= 0){
                neighbors.push_back({it->first, it->second});
            }
        }
        return neighbors;
    }

 private:
    std::map<int, std::map<int,id>> graph;
};


#endif //__EDHEGRAPH_HPP__