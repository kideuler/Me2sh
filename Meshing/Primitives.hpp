#ifndef __PRIMITIVES_HPP__
#define __PRIMITIVES_HPP__

#include <array>
#include <vector>
#include <iostream>
#include <stack>
#include <fstream>
#include <cassert>
#include <functional>
#include <queue>
#include <unordered_set>

#include "Graph.hpp"

#define NUM_EDGE_CONNECTIONS 10

// class forwarding
class vertex;
class edge;
class triangle;
class triangulation;

// vertex class
class vertex {
 public:
    // default pposition
    vertex(){}

    vertex(std::array<double, 2> V){
        coord = V;
    }

    void operator=(std::array<double, 2> V){
        coord = V;
    } 
    
    // giving constant so as to compare
    const double& operator[](int v) const {
        return coord[v];
    }
    
    // access to change
    double& operator[](int v){
        return coord[v];
    }

    double distance(const vertex &other){
        double x = other.coord[0] - coord[0];
        double y = other.coord[1] - coord[1];
        return sqrt(x*x + y*y);
    }

 private:
    std::array<double, 2> coord;
    std::array<int, NUM_EDGE_CONNECTIONS> edges;
};


// edge class
class edge {
 public:
    // default vertices
    edge(){vertices = {-1,-1}; triangles = {-1,-1}; lids = {-1,-1};}

    edge(int v1, int v2) {
        vertices = {v1, v2};
        triangles = {-1, -1};
        lids = {-1, -1};
    }

    edge(std::array<int, 2> V) {
        vertices = V;
        triangles = {-1, -1};
        lids = {-1, -1};
    }

    // giving constant to compare
    const int& operator[](int v) const {
        return vertices[v];
    }

    // access to change
    int& operator[](int v){
        return vertices[v];
    }

    // giving constant to compare
    const int& operator()(int t) const {
        return triangles[t];
    }

    // access to change
    int& operator()(int t){
        return triangles[t];
    }

    // giving constant to compare
    const int& L(int t) const {
        return lids[t];
    }

    int& L(int t){
        return lids[t];
    }

    bool is_equal(const edge &other){
        return (vertices[0] == other[0] && vertices[1] == other[1]) || (vertices[0] == other[1] && vertices[1] == other[0]);
    }

    bool& on_boundary(){
        return bdy;
    }



 private:
    std::array<int,2> vertices;
    std::array<int,2> triangles;
    std::array<int,2> lids; // local edges on the triangle
    bool bdy = false;
};

// triangle class
class triangle {
 public:
    triangle(){}

    triangle(std::array<int,3> V){
        vertices = V;
    }



    // giving constant to compare
    const int& operator[](int v) const {
        return vertices[v];
    }

    // access to change
    int& operator[](int v){
        return vertices[v];
    }

    // giving constant to compare
    const int& operator()(int e) const {
        return edges[e];
    }

    // access to change
    int& operator()(int e){
        return edges[e];
    }

    int& O(int e){
        return orientation[e];
    }

 private:
    std::array<int,3> vertices;
    std::array<int,3> edges;
    std::array<int,3> orientation;
};

// triangulation class
class triangulation {

 public:
    triangulation(){}

    triangulation(std::vector<std::array<int, 3> > elems, std::vector<std::array<double, 2> > coords){
        int nv = coords.size();
        this->vertices.resize(nv);
        for (int v = 0; v<nv; ++v){
            this->vertices[v] = coords[v];
        }

        int nelems = elems.size();
        this->triangles.resize(nelems);
        for (int t = 0; t<nelems; ++t){
            this->triangles[t] = elems[t];
        }
        
        int e, ee;
        edge E;
        bool found;
        for (int t = 0; t<nelems; ++t){

            for (int i = 0; i<3; i++){
                e = 0;
                found = false;
                E[0] = triangles[t][i];
                E[1] = triangles[t][(i+1)%3];
                while (e<edges.size() && !found){

                    if (edges[e].is_equal(E)){
                        found = true;
                        ee = e;
                    }
                    e++;
                }

                if (found){
                    edges[ee](1) = t;
                    edges[ee].L(1) = i;
                    triangles[t](i) = ee;
                    triangles[t].O(i) = 1;
                } else {
                    edge E2;
                    E2[0] = E[0];
                    E2[1] = E[1];
                    E2(0) = t;
                    E2.L(0) = i;
                    triangles[t](i) = edges.size();
                    triangles[t].O(i) = 0;
                    edges.push_back(E2);
                }
            }
        }
    }

    // edge length
    double EdgeLength(int e){
        int v1 = edges[e][0];
        int v2 = edges[e][1];

        return vertices[v1].distance(vertices[v2]);
    }

    // triangle area
    double TriangleArea(int t){
        int v1 = triangles[t][0];
        int v2 = triangles[t][1];
        int v3 = triangles[t][2];

        double x1 = vertices[v1][0];
        double y1 = vertices[v1][1];
        double x2 = vertices[v2][0];
        double y2 = vertices[v2][1];
        double x3 = vertices[v3][0];
        double y3 = vertices[v3][1];

        return 0.5 * (x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2));
    }

    // triangle circumradius
    double TriangleCircumradius(int t){

        double A = TriangleArea(t);
        int e1 = triangles[t](0);
        int e2 = triangles[t](1);
        int e3 = triangles[t](2);

        double a = EdgeLength(e1);
        double b = EdgeLength(e2);
        double c = EdgeLength(e3);

        return (a*b*c)/(4*A);
    }

    vertex TriangleCircumcenter(int t){
        int v1 = triangles[t][0];
        int v2 = triangles[t][1];
        int v3 = triangles[t][2];

        double x1 = vertices[v1][0];
        double y1 = vertices[v1][1];
        double x2 = vertices[v2][0];
        double y2 = vertices[v2][1];
        double x3 = vertices[v3][0];
        double y3 = vertices[v3][1];

        double D = 2*(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
        double ux = (x1*x1 + y1*y1)*(y2-y3) + \
            (x2*x2 + y2*y2)*(y3-y1) + \
            (x3*x3 + y3*y3)*(y1-y2);
        double uy = (x1*x1 + y1*y1)*(x3-x2) + \
            (x2*x2 + y2*y2)*(x1-x3) + \
            (x3*x3 + y3*y3)*(x2-x1);

        vertex C;
        C[0] = ux/D; C[1] = uy/D;
        return C;
    }

    int TriangleOpposite_tri(int t, int e){
        int o = triangles[t].O(e);
        int edge = triangles[t](e);
        return edges[edge](1-o);
    }

    int TriangleOpposite_edge(int t, int e){
        int edge = triangles[t](e);
        int o = triangles[t].O(e);
        return edges[edge].L(1-o);
    }

    void EdgeFlip(int e){
        // Get the two triangles sharing the edge
        int t1 = edges[e](0);
        int t2 = edges[e](1);
        int v1e = edges[e][0];
        int v2e = edges[e][1];
        int e1_index = edges[e].L(0);
        int e2_index = edges[e].L(1);

        if (t1 < 0 || t2 < 0){
            return;
        }

        // Get the vertices of the quad connecting the triangles
        int v1 = (triangles[t1][0] != v1e && triangles[t1][0] != v2e) ? triangles[t1][0] : (triangles[t1][1] != v1e && triangles[t1][1] != v2e) ? triangles[t1][1] : triangles[t1][2];
        int v2 = v1e;
        int v3 = (triangles[t2][0] != v1e && triangles[t2][0] != v2e) ? triangles[t2][0] : (triangles[t2][1] != v1e && triangles[t2][1] != v2e) ? triangles[t2][1] : triangles[t2][2];
        int v4 = v2e;

        // Update the triangle vertices
        triangles[t1][0] = v1;
        triangles[t1][1] = v2;
        triangles[t1][2] = v3;

        triangles[t2][0] = v3;
        triangles[t2][1] = v4;
        triangles[t2][2] = v1;

        edges[e][0] = v3;
        edges[e][1] = v1;
        edges[e].L(0) = 2;
        edges[e].L(1) = 2;

        // find quad edges
        int e1 = triangles[t1]((e1_index + 2) % 3);
        int e2 = triangles[t2]((e2_index + 1) % 3);
        int e3 = triangles[t2]((e2_index + 2) % 3);
        int e4 = triangles[t1]((e1_index + 1) % 3);

        // find quad orientation
        int o1 = triangles[t1].O((e1_index + 2) % 3);
        int o2 = triangles[t2].O((e2_index + 1) % 3);
        int o3 = triangles[t2].O((e2_index + 2) % 3);
        int o4 = triangles[t1].O((e1_index + 1) % 3);

        // re-assigning triangles to edges
        triangles[t1](0) = e1;
        triangles[t1].O(0) = o1;
        triangles[t1](1) = e2;
        triangles[t1].O(1) = o2;
        triangles[t1](2) = e;
        triangles[t1].O(2) = 0;

        triangles[t2](0) = e3;
        triangles[t2].O(0) = o3;
        triangles[t2](1) = e4;
        triangles[t2].O(1) = o4;
        triangles[t2](2) = e;
        triangles[t2].O(2) = 1;

        // reassigning edges corresponding triangles
        edges[e1](o1) = t1;
        edges[e1].L(o1) = 0;
        edges[e2](o2) = t1;
        edges[e2].L(o2) = 1;
        edges[e3](o3) = t2;
        edges[e3].L(o3) = 0;
        edges[e4](o4) = t2;
        edges[e4].L(o4) = 1;

    }

    bool EdgeIsDelaunay(int e){
        // Get the two triangles sharing the edge
        int t1 = edges[e](0);
        int t2 = edges[e](1);
        int v1e = edges[e][0];
        int v2e = edges[e][1];

        // If the edge is a boundary edge, it is Delaunay by definition
        if (t1 == -1 || t2 == -1) {
            return true;
        }

        vertex C = TriangleCircumcenter(t1);
        double R = TriangleCircumradius(t1);

        int v4 = triangles[t2][0];
        int v5 = triangles[t2][1];
        int v6 = triangles[t2][2];
        int opp = (v4 != v1e && v4 != v2e) ? v4 : (v5 != v1e && v5 != v2e) ? v5 : v6;
        
        return vertices[opp].distance(C) >= R;
    }

    void TriangleFlip(int t, int e){
        EdgeFlip(triangles[t](e));
    }

    bool TriangleContainsPoint(int t, const vertex &p){
        int v1 = triangles[t][0];
        int v2 = triangles[t][1];
        int v3 = triangles[t][2];

        double v0x = vertices[v3][0] - vertices[v1][0];
        double v0y = vertices[v3][1] - vertices[v1][1];
        double v1x = vertices[v2][0] - vertices[v1][0];
        double v1y = vertices[v2][1] - vertices[v1][1];
        double v2x = p[0] - vertices[v1][0];
        double v2y = p[1] - vertices[v1][1];

        double dot00 = v0x*v0x + v0y*v0y;
        double dot01 = v0x*v1x + v0y*v1y;
        double dot02 = v0x*v2x + v0y*v2y;
        double dot11 = v1x*v1x + v1y*v1y;
        double dot12 = v1x*v2x + v1y*v2y;

        // compute barycentric
        double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
        double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        return (u >= 0.0) && (v >= 0.0) && (u + v <= 1.0);
    }

    int NextTriangle(int currentTriangle, const vertex& point) {
        const triangle& t = triangles[currentTriangle];
        for (int i = 0; i < 3; ++i) {
            int neighbor = TriangleOpposite_tri(currentTriangle,i);
            if (neighbor != -1) {
                const vertex& A = vertices[t[i]];
                const vertex& B = vertices[t[(i + 1) % 3]];
                double crossProduct = (B[0] - A[0]) * (point[1] - A[1]) - (B[1] - A[1]) * (point[0] - A[0]);
                if (crossProduct < 0) {
                    return neighbor;
                }
            }
        }
        return -1;
    }

    int FindTriangleContainingPoint(int tstart, const vertex &p){

        int currentTriangle = tstart;
        while (true) {
            if (TriangleContainsPoint(currentTriangle, p) ){
                return currentTriangle;
            }
            int nextTriangle = NextTriangle(currentTriangle, p);
            if (nextTriangle == -1) {
                break;
            }
            currentTriangle = nextTriangle;
        }
        return -1;  // Point is outside the triangulation
    }

    void InsertPoint(const int &p, int tstart = 0);
    
    void DelaunayTriangulation(std::vector<std::array<double, 2> > coords);

    bool IsQuadConvex(int e){
        int t1 = edges[e](0);
        int t2 = edges[e](1);

        int e1_index = edges[e].L(0);
        int e2_index = edges[e].L(1);

        // Get the vertices of the quad connecting the triangles
        int v1 = triangles[t1][(e1_index + 2) % 3];
        int v2 = triangles[t1][(e1_index)];
        int v3 = triangles[t2][(e2_index + 2) % 3];
        int v4 = triangles[t2][(e2_index)];

        // Helper function to compute the cross product of two vectors
        auto cross = [](const std::array<double, 2>& O, const std::array<double, 2>& A, const std::array<double, 2>& B) {
            return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]);
        };

        // Get the coordinates of the vertices
        std::array<double, 2> p1 = {vertices[v1][0], vertices[v1][1]};
        std::array<double, 2> p2 = {vertices[v2][0], vertices[v2][1]};
        std::array<double, 2> p3 = {vertices[v3][0], vertices[v3][1]};
        std::array<double, 2> p4 = {vertices[v4][0], vertices[v4][1]};

        // Check the orientation of the quadrilateral
        bool convex = (cross(p1, p2, p3) > 0) && (cross(p2, p3, p4) > 0) && (cross(p3, p4, p1) > 0) && (cross(p4, p1, p2) > 0);

        return convex;
    }

    bool IsRayIntersectingEdge(const int vr1, const int vr2, const int e1, const int e2){
        double x1 = vertices[vr1][0];
        double y1 = vertices[vr1][1];
        double x2 = vertices[vr2][0];
        double y2 = vertices[vr2][1];
        double x3 = vertices[e1][0];
        double y3 = vertices[e1][1];
        double x4 = vertices[e2][0];
        double y4 = vertices[e2][1];

        // Compute the direction vectors
        double dx1 = x2 - x1;
        double dy1 = y2 - y1;
        double dx2 = x4 - x3;
        double dy2 = y4 - y3;

        // Compute the denominator of the intersection formula
        double denominator = dx1 * dy2 - dy1 * dx2;

        // If the denominator is zero, the lines are parallel
        if (std::abs(denominator) < 1e-10) {
            return false;
        }

        // Compute the intersection parameters
        double t = ((x3 - x1) * dy2 - (y3 - y1) * dx2) / denominator;
        double u = ((x3 - x1) * dy1 - (y3 - y1) * dx1) / denominator;

        // Check if the intersection parameters are within the valid range
        return (t >= 0) && (u >= 0) && (u <= 1);
    }

    void DeleteExteriorElements(const std::vector<std::array<int, 2>>& segments) {
        
        // finding centroid of all segments
        double x = 0.0;
        double y = 0.0;
        for (auto segment : segments) {
            x += vertices[segment[0]][0] + vertices[segment[1]][0];
            y += vertices[segment[0]][1] + vertices[segment[1]][1];
        }
        x /= 2 * segments.size();
        y /= 2 * segments.size();

        // find the triangle containing the centroid
        int t = FindTriangleContainingPoint(0, vertex({x, y}));

        // find all triangles connected to the centroid triangle which do not cross segment edges
        std::queue<int> q;
        std::vector<bool> visited(triangles.size(), false);
        q.push(t);
        visited[t] = true;
        while (!q.empty()) {
            int currentTriangle = q.front();
            q.pop();
            for (int i = 0; i < 3; ++i) {
                int neighbor = TriangleOpposite_tri(currentTriangle, i);
                if (neighbor != -1 && !visited[neighbor] && !edges_on_boundary[triangles[currentTriangle](i)]) {
                    q.push(neighbor);
                    visited[neighbor] = true;
                }
            }
        }

        // delete all triangles not connected to the centroid triangle
        std::vector<triangle> new_triangles;
        std::vector<edge> new_edges;
        std::vector<bool> new_edges_on_boundary;
        std::vector<vertex> new_vertices;
        std::unordered_set<int> valid_edges;
        for (int i = 0; i < triangles.size(); ++i) {
            if (visited[i]) {
                new_triangles.push_back(triangles[i]);
                for (int j = 0; j < 3; ++j) {
                    int edge = triangles[i](j);
                    if (valid_edges.find(edge) == valid_edges.end()) {
                        valid_edges.insert(edge);
                        new_edges.push_back(edges[edge]);
                        new_edges_on_boundary.push_back(edges_on_boundary[edge]);
                    }
                }
            }
        }

        // update the boundary edges
        edges_on_boundary = new_edges_on_boundary;

        // update the data structures
        edges = new_edges;
        triangles = new_triangles;
    }

    void DelaunayTriangluationConstrained(std::vector<std::array<double, 2> > coords, std::vector<std::array<int,2>> segments){

        // create initial triangulation
        DelaunayTriangulation(coords);
        
        // create the edge graph from the triangulation
        for (auto edge : edges){
            edgeGraph(edge[0], edge[1],true);
            edgeGraph(edge[1], edge[0],true);
        }

        // create the vertex-triangle graph
        for (int i = 0; i<triangles.size(); i++){
            for (int j = 0; j<3; j++){
                vertTriangleGraph(triangles[i][j],i,j);
            }
        }

        // create the constrained edges
        bool connected,inray,flipped;
        std::stack<std::array<int,2>> triedges;
        for (auto segment : segments){
            int v1 = segment[0];
            int v2 = segment[1];
            
            // check if the edge is already present
            if (edgeGraph(v1,v2)){
                continue;
            }

            while (!edgeGraph(v1,v2)){
                std::vector<std::array<int,2>> ngbrs = vertTriangleGraph(v1);
                int kk = 0;
                inray = false;
                while (kk<ngbrs.size() && !inray){
                    int t = ngbrs[kk][0];
                    int v = ngbrs[kk][1];
                    int e = (v+1)%3;
                    int vt0 = triangles[t][v];
                    int vt1 = triangles[t][(v+1)%3];
                    int vt2 = triangles[t][(v+2)%3];
                    int ve1 = edges[triangles[t](e)][0];
                    int ve2 = edges[triangles[t](e)][1];

                    // check if the edge is in the ray casted by the segment
                    if (IsRayIntersectingEdge(v1,v2,ve1,ve2)){
                        inray = true;
                        triedges.push({t,e});
                    }
                    kk++;
                }

                flipped = false;
                while (!triedges.empty() && !flipped){
                    std::array<int,2> triedge = triedges.top();
                    triedges.pop();
                    int t = triedge[0];
                    int e = triedge[1];
                    int o = triangles[t].O(e);
                    int edge = triangles[t](e);
                    if (IsQuadConvex(edge)){
                        // turn off edges in the edge graph and vert-tri graph
                        int e1 = edges[edge][0];
                        int e2 = edges[edge][1];
                        edgeGraph(e1,e2,false);
                        edgeGraph(e2,e1,false);

                        for (int i = 0; i<2; i++){
                            int tg = edges[edge](i);
                            if (tg != -1){
                                for (int j = 0; j<3; j++){
                                    int vg = triangles[tg][j];
                                    vertTriangleGraph(vg,tg,-1);
                                }
                            }
                        }

                        // flip the edge
                        EdgeFlip(edge);

                        // turn on edges in the edge graph and vert-tri graph
                        e1 = edges[edge][0];
                        e2 = edges[edge][1];
                        edgeGraph(e1,e2,true);
                        edgeGraph(e2,e1,true);

                        for (int i = 0; i<2; i++){
                            int tg = edges[edge](i);
                            if (tg != -1){
                                for (int j = 0; j<3; j++){
                                    int vg = triangles[tg][j];
                                    vertTriangleGraph(vg,tg,j);
                                }
                            }
                        }

                        flipped = true;
                    } else {
                        std::cout << "Edge " << triangles[t](e) << " is not convex" << std::endl;
                        // checking other edges of opposite triangle to see if they are in the ray
                        t = edges[edge](1-o);
                        int leid = edges[edge].L(1-o);
                        for (int i = 1; i<=2; i++){
                            int e = triangles[t]((leid+i)%3);
                            int ve1 = edges[e][0];
                            int ve2 = edges[e][1];
                            if (IsRayIntersectingEdge(v1,v2,ve1,ve2)){
                                triedges.push({t,(leid+i)%3});
                            }
                        }
                    }
                }
            }
        }

        edges_on_boundary.resize(edges.size(),false);
        for (auto segment : segments){
            int v1 = segment[0];
            int v2 = segment[1];
            for (int e = 0; e<edges.size(); e++){
                if ((edges[e][0] == v1 && edges[e][1] == v2) || (edges[e][0] == v2 && edges[e][1] == v1)){
                    edges_on_boundary[e] = true;
                }
            }
        }

        // find and delete all elements on the opposite side of the boundary
        DeleteExteriorElements(segments);
    }


    void DelaunayRefinement(double h_target, int npasses=2);

    void Print(){
        std::cout << "\nCOORDINATES" << std::endl;
        for (auto v : vertices){
            std::cout << v[0] << "," << v[1] << std::endl;
        }

        std::cout << "\nEDGES" << std::endl;
        for (auto e : edges){
            std::cout << e[0] << "," << e[1] << std::endl;
        }

        std::cout << "\nTRIANGLES" << std::endl;
        for (auto t : triangles){
            std::cout << t[0] << "," << t[1] << "," << t[2] << std::endl;
        }
    }

    void writeVTK(const std::string &filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open file for writing");
        }

        // Write VTK header
        file << "# vtk DataFile Version 3.0\n";
        file << "Triangulation\n";
        file << "ASCII\n";
        file << "DATASET UNSTRUCTURED_GRID\n";

        // Write points
        file << "POINTS " << vertices.size() << " float\n";
        for (const auto &coord : vertices) {
            file << coord[0] << " " << coord[1] << " 0.0\n";
        }

        // Write cells
        file << "CELLS " << triangles.size() << " " << 4 * triangles.size() << "\n";
        for (const auto &triangle : triangles) {
            file << "3 " << triangle[0] << " " << triangle[1] << " " << triangle[2] << "\n";
        }

        // Write cell types
        file << "CELL_TYPES " << triangles.size() << "\n";
        for (size_t i = 0; i < triangles.size(); ++i) {
            file << "5\n"; // VTK_TRIANGLE
        }

        file.close();
    }

    std::vector<triangle> triangles;
    std::vector<edge> edges;
    std::vector<vertex> vertices;
    std::vector<bool> edges_on_boundary;
    std::stack<int> flip_Stack;
    Graph edgeGraph;
    GraphCRS vertTriangleGraph;
};

#endif //__PRIMITIVES_HPP__