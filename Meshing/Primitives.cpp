#include "Primitives.hpp"


void triangulation::InsertPoint(const int &p, int tstart){

    // insert point into triangle and split triangle
    int ntris = triangles.size();
    int t = FindTriangleContainingPoint(tstart,vertices[p]);

    int v1 = triangles[t][0];
    int v2 = triangles[t][1];
    int v3 = triangles[t][2];
    int e1 = triangles[t](0);
    int e2 = triangles[t](1);
    int e3 = triangles[t](2);
    int o1 = triangles[t].O(0);
    int o2 = triangles[t].O(1);
    int o3 = triangles[t].O(2);

    triangles[t][0] = v1;
    triangles[t][1] = v2;
    triangles[t][2] = p;
    edges[e1](o1) = t;
    edges[e1].L(o1) = 0;

    triangle T2({v2,v3,p});
    triangles.push_back(T2);
    edges[e2](o2) = ntris;
    edges[e2].L(o2) = 0;

    triangle T3({v3,v1,p});
    triangles.push_back(T3);
    edges[e3](o3) = ntris+1;
    edges[e3].L(o3) = 0;

    // Create new edges
    int e4 = edges.size();
    edges.push_back(edge{v2, p});
    edges[e4](0) = t;
    edges[e4].L(0) = 1;
    edges[e4](1) = ntris;
    edges[e4].L(1) = 2;

    int e5 = edges.size();
    edges.push_back(edge{v3, p});
    edges[e5](0) = ntris;
    edges[e5].L(0) = 1;
    edges[e5](1) = ntris + 1;
    edges[e5].L(1) = 2;

    int e6 = edges.size();
    edges.push_back(edge{v1, p});
    edges[e6](0) = ntris + 1;
    edges[e6].L(0) = 1;
    edges[e6](1) = t;
    edges[e6].L(1) = 2;

    // Update the new triangles with the new edges
    triangles[t](0) = e1;
    triangles[t](1) = e4;
    triangles[t](2) = e6;
    triangles[t].O(0) = o1;
    triangles[t].O(1) = 0;
    triangles[t].O(2) = 1;

    triangles[ntris](0) = e2;
    triangles[ntris](1) = e5;
    triangles[ntris](2) = e4;
    triangles[ntris].O(0) = o2;
    triangles[ntris].O(1) = 0;
    triangles[ntris].O(2) = 1;

    triangles[ntris + 1](0) = e3;
    triangles[ntris + 1](1) = e6;
    triangles[ntris + 1](2) = e5;
    triangles[ntris + 1].O(0) = o3;
    triangles[ntris + 1].O(1) = 0;
    triangles[ntris + 1].O(2) = 1;


    flip_Stack.push(e1);
    flip_Stack.push(e2);
    flip_Stack.push(e3);

    

    int E,t1,t2;
    while (!flip_Stack.empty()){
        E = flip_Stack.top(); 
        flip_Stack.pop();
        if (!EdgeIsDelaunay(E)){
            EdgeFlip(E);
            t1 = edges[E](0);
            t2 = edges[E](1);
            
            flip_Stack.push(triangles[t1](0));
            flip_Stack.push(triangles[t1](1));
            flip_Stack.push(triangles[t2](0));
            flip_Stack.push(triangles[t2](1));
        }
    }
}

void triangulation::DelaunayTriangulation(std::vector<std::array<double, 2>> coords) {
    // Clear existing data
    vertices.clear();
    triangles.clear();
    edges.clear();

    // Find the bounding box of the points
    double minX = coords[0][0], maxX = coords[0][0];
    double minY = coords[0][1], maxY = coords[0][1];
    for (const auto& coord : coords) {
        if (coord[0] < minX) minX = coord[0];
        if (coord[0] > maxX) maxX = coord[0];
        if (coord[1] < minY) minY = coord[1];
        if (coord[1] > maxY) maxY = coord[1];
    }
    double dx = maxX - minX;
    double dy = maxY - minY;

    std::srand(std::time(0)); // Seed for random number generation

    // Add all points to the vertices list
    for (auto coord : coords) {
        coord[0] += 1e-4*dx*(static_cast<double>(std::rand()) / RAND_MAX);
        coord[1] += 1e-4*dy*(static_cast<double>(std::rand()) / RAND_MAX);
        vertices.push_back(coord);
    }

    // Create a super triangle that encompasses all the points
    double deltaMax = std::max(dx, dy);
    double midX = (minX + maxX) / 2.0;
    double midY = (minY + maxY) / 2.0;

    vertex v1({midX - 20 * deltaMax, midY - deltaMax});
    vertex v2({midX + 20 * deltaMax, midY - deltaMax});
    vertex v3({midX, midY + 20 * deltaMax});
    
    int nv = vertices.size();
    vertices.push_back(v1);
    vertices.push_back(v2);
    vertices.push_back(v3);

    // Add the super triangle to the triangulation
    triangle tri_0({nv, nv+1, nv+2});
    triangles.push_back(tri_0);
    triangles[0](0) = 0;
    triangles[0](1) = 1;
    triangles[0](2) = 2;
    triangles[0].O(0) = 0;
    triangles[0].O(1) = 0;
    triangles[0].O(2) = 0;

    // Create edges for the super triangle
    edge e1(nv, nv+1);
    edge e2(nv+1, nv+2);
    edge e3(nv+2, nv);
    edges.push_back(e1);
    edges.push_back(e2);
    edges.push_back(e3);
    edges[0](0) = 0;
    edges[0](1) = -1;
    edges[0].L(0) = 0;
    edges[1](0) = 0;
    edges[1](1) = -1;
    edges[1].L(0) = 1;
    edges[2](0) = 0;
    edges[2](1) = -1;
    edges[2].L(0) = 2;


    // Insert each point into the triangulation
    for (int i = 0; i < nv; ++i) {
        InsertPoint(i, triangles.size()-1);
    }

    // Remove triangles connected to the super triangle vertices
    int ntris = triangles.size();
    std::vector<int> map(ntris,-1);
    std::vector<triangle> newTriangles;
    int nt = 0;
    for (int i = 0; i < ntris; ++i) {
        if (triangles[i][0] >= nv || triangles[i][1] >= nv || triangles[i][2] >= nv) {
            continue;
        }
        newTriangles.push_back(triangles[i]);
        map[i] = nt;
        nt++;
    }
    triangles = newTriangles;

    // Update the edges
    int ne = edges.size();
    std::vector<edge> newEdges;
    std::vector<int> mapedge(ne);
    int ne_new = 0;
    for (int i = 0; i < ne; ++i) {
        if (edges[i][0] >= nv || edges[i][1] >= nv) {
            continue;
        }
        edges[i](0) = map[edges[i](0)];
        edges[i](1) = map[edges[i](1)];
        newEdges.push_back(edges[i]);
        mapedge[i] = ne_new;
        ne_new++;
    }
    edges = newEdges;

    // update the triangle edges
    for (int i = 0; i < nt; ++i) {
        for (int j = 0; j < 3; ++j) {
            triangles[i](j) = mapedge[triangles[i](j)];
        }
    }

    // Remove the super triangle vertices from the vertices list
    vertices.resize(vertices.size() - 3);
}


void triangulation::DelaunayRefinement(double h_target, int npasses){ 
    double r_target =(sqrt(3)/3)*h_target;
    std::queue<int> Q;


    for (int np = 0; np<npasses; np++){

        // create queue with all available triangles
        for (int i = 0; i < triangles.size(); ++i) {
            Q.push(i);
        }

        // while queue is not empty refine triangles
        int nv = vertices.size();
        while (!Q.empty()){
            int t = Q.front();
            Q.pop();
            double R = TriangleCircumradius(t);
            if (R > h_target){
                vertex C = TriangleCircumcenter(t);

                int tri = FindTriangleContainingPoint(t,C);

                if (tri == -1){
                    continue;
                }
                
                vertices.push_back(C);
                InsertPoint(nv,t);
                nv++;
                Q.push(tri);
                Q.push(triangles.size()-1);
                Q.push(triangles.size()-2);
            }
        }
    }
}