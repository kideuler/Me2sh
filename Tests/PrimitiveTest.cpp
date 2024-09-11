#include "Primitives.hpp"
#include <nlohmann/json.hpp> // Include the nlohmann/json library

#define TOL 1e-12

const double PI = 3.14159265358979323846;

#define ASSERT_NEAR(a,b)(abs(a-b)<TOL)

std::pair<std::vector<std::array<double, 2> >, std::vector<std::array<int, 3> > > constructTriangulation(int nx, int ny) {
    std::vector<std::array<double, 2> > coords;
    std::vector<std::array<int, 3> > triangles;

    // Create the coordinates
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            coords.push_back({i * 1.0 / nx, j * 1.0 / ny});
        }
    }

    // Create the triangles
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int p0 = j * (nx + 1) + i;
            int p1 = p0 + 1;
            int p2 = p0 + (nx + 1);
            int p3 = p2 + 1;

            // First triangle of the cell
            triangles.push_back({p0, p1, p2});
            // Second triangle of the cell
            triangles.push_back({p1, p3, p2});
        }
    }

    return {coords, triangles};
}

std::vector<std::array<double, 2>> createUnitCirclePoints(int numPoints) {
    
    std::vector<std::array<double, 2>> points;

    for (int i = 0; i < numPoints; ++i) {
        double angle = 2 * PI * i / numPoints;
        points.push_back({std::cos(angle), std::sin(angle)});
    }

    return points;
}

std::vector<std::array<double, 2>> createEllipsePoints(int numPoints, double a, double b) {
    
    std::vector<std::array<double, 2>> points;

    for (int i = 0; i < numPoints; ++i) {
        double angle = 2 * PI * i / numPoints;
        points.push_back({a*std::cos(angle), b*std::sin(angle)});
    }

    return points;
}

std::array<double, 2> computeCentroid(const std::vector<std::array<double, 2>>& points) {
    double cx = 0.0, cy = 0.0;
    for (const auto& point : points) {
        cx += point[0];
        cy += point[1];
    }
    cx /= points.size();
    cy /= points.size();
    return {cx, cy};
}

// Function to generate a random parameterized closed shape
std::vector<std::array<double, 2>> generateRandomClosedShape(int numPoints, double radius = 1.0) {
    std::vector<std::array<double, 2>> points(numPoints);

    // Seed the random number generator
    std::srand(std::time(0));

    // Generate random points within a circle of given radius
    for (int i = 0; i < numPoints; ++i) {
        double angle = 2 * PI * std::rand() / RAND_MAX;
        double r = radius * std::sqrt(static_cast<double>(std::rand()) / RAND_MAX);
        points[i] = {r * std::cos(angle), r * std::sin(angle)};
    }

    // Compute the centroid of the points
    std::array<double, 2> centroid = computeCentroid(points);

    // Sort points by angle around the centroid to form a closed shape
    std::sort(points.begin(), points.end(), [&centroid](const std::array<double, 2>& a, const std::array<double, 2>& b) {
        double angleA = std::atan2(a[1] - centroid[1], a[0] - centroid[0]);
        double angleB = std::atan2(b[1] - centroid[1], b[0] - centroid[0]);
        return angleA < angleB;
    });

    return points;
}

std::pair<std::vector<std::array<double, 2>>, std::vector<std::array<int, 2>>> readLakeSuperiorData(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    int numPoints;
    file >> numPoints;

    std::vector<std::array<double, 2>> points(numPoints);
    for (int i = 0; i < numPoints; ++i) {
        file >> points[i][0] >> points[i][1];
        // Skip the third value (z-coordinate) as it's not needed
        double z;
        file >> z;
    }

    int numSegments;
    file >> numSegments;

    std::vector<std::array<int, 2>> segments(numSegments);
    for (int i = 0; i < numSegments; ++i) {
        file >> segments[i][0] >> segments[i][1];
        // Skip the third value as it's not needed
        int dummy;
        file >> dummy;
    }

    return {points, segments};
}

using json = nlohmann::json;

std::pair<std::vector<std::array<double, 2>>, std::vector<std::array<int, 2>>> readGeoJSON(const std::string& filename) {
    std::ifstream file(filename);
    json geojson;
    file >> geojson;

    std::vector<std::array<double, 2>> coords;
    std::vector<std::array<int, 2>> segments;

    for (const auto& geometry : geojson["geometries"]) {
        if (geometry["type"] == "Polygon") {
            const auto& polygon = geometry["coordinates"][0];
            int start_index = coords.size();
            for (const auto& point : polygon) {
                coords.push_back({point[0], point[1]});
            }
            int end_index = coords.size();
            for (int i = start_index; i < end_index - 1; ++i) {
                segments.push_back({i, i + 1});
            }
            segments.push_back({end_index - 1, start_index});
        }
    }

    return {coords, segments};
}

int TEST1(){
    std::vector<std::array<double,2> > coords(4);
    coords[0][0] = 0.0; coords[0][1] = 0.0;
    coords[1][0] = 1.0; coords[1][1] = 0.0;
    coords[2][0] = 0.0; coords[2][1] = 1.0;
    coords[3][0] = 1.0; coords[2][1] = 1.0;

    std::vector<std::array<int,3> > tri(2);
    tri[0][0] = 0;
    tri[0][1] = 1;
    tri[0][2] = 2;

    tri[1][0] = 1;
    tri[1][1] = 3;
    tri[1][2] = 2;

    triangulation T(tri,coords);
    
    if (!ASSERT_NEAR(T.TriangleArea(0),0.5)){return 1;}

    if (!ASSERT_NEAR(T.TriangleCircumradius(0),sqrt(2)/2)){return 1;}

    vertex C = T.TriangleCircumcenter(0);
    if (!ASSERT_NEAR(C[0],0.5) || !ASSERT_NEAR(C[1],0.5)){return 1;}

    if (!(T.TriangleOpposite_tri(0,1)==1 && T.TriangleOpposite_tri(1,2) == 0)){return 1;}

    vertex p({0.3,0.3});
    if (!T.TriangleContainsPoint(0,p)){return 1;}

    return 0;
}

int TEST2(){
    auto input = constructTriangulation(3, 3);

    triangulation T(input.second,input.first);

    std::cout << T.TriangleOpposite_tri(0,0) << " " << T.TriangleOpposite_tri(0,1) << " " << T.TriangleOpposite_tri(0,2) << std::endl;
    if (!(T.TriangleOpposite_tri(0,0) == -1 && T.TriangleOpposite_tri(0,1) == 1 && T.TriangleOpposite_tri(0,2) == -1)){return 1;}
    

    std::cout << T.TriangleOpposite_tri(8,0) << " " << T.TriangleOpposite_tri(8,1) << " " << T.TriangleOpposite_tri(8,2) << std::endl;
    if (!(T.TriangleOpposite_tri(8,0)==3 && T.TriangleOpposite_tri(8,1)==9 && T.TriangleOpposite_tri(8,2)==7)){return 1;}
    

    int t1 = 8;
    int edge = T.triangles[t1](1);
    int t1v1 = T.triangles[t1][0]; int t1v2 = T.triangles[t1][1]; int t1v3 = T.triangles[t1][2];
    int t1e1 = T.triangles[t1](0); int t1e2 = T.triangles[t1](1); int t1e3 = T.triangles[t1](2);

    int t2 = T.TriangleOpposite_tri(t1,1);
    int eopp = T.TriangleOpposite_edge(t1,1);
    int t2opp = T.triangles[t2][(eopp+1)%3];
    int t2e1 = T.triangles[t2](eopp); int t2e2 = T.triangles[t2]((eopp+1)%3); int t2e3 = T.triangles[t2]((eopp+2)%3);


    T.EdgeFlip(T.triangles[t1](1));
    // check all tris have psoitive area
    for (int i = 0; i<T.triangles.size(); i++){
        if (T.TriangleArea(i)<0){
            std::cerr << "Triangle " << i << " has negative area." << std::endl;
            return 1;
        }
    }

    // checking vertices
    std::cout << T.triangles[t1][0] << " " << T.triangles[t1][1] << " " << T.triangles[t1][2] << std::endl;
    if (!(T.triangles[t1][0]==t1v1 && T.triangles[t1][1]==t1v2 && T.triangles[t1][1]==t2opp)){return 1;}

    // checking edges
    std::cout << T.triangles[t1](0) << " " << T.triangles[t1](1) << " " << T.triangles[t1](2) << std::endl;
    if (!(T.triangles[t1](0) == t1e1 && T.triangles[t1](1) == t2e2 && T.triangles[t1](2) == edge)){return 1;}

    //checking edge -> triangles
    int e1 = T.triangles[t1](0); int e2 = T.triangles[t1](1); int e3 = T.triangles[t1](2);
    int o1 = T.triangles[t1].O(0); int o2 = T.triangles[t1].O(1); int o3 = T.triangles[t1].O(2);
    std::cout << e1 << " " << e2 << " " << e3 <<std::endl;
    if (!(T.edges[e1](o1) == t1 && T.edges[e2](o2) == t1 && T.edges[e3](o3) == t1)){return 1;}

    e1 = T.triangles[t2](0); e2 = T.triangles[t2](1); e3 = T.triangles[t2](2);
    o1 = T.triangles[t2].O(0); o2 = T.triangles[t2].O(1); o3 = T.triangles[t2].O(2);
    std::cout << T.edges[e1](o1) << " " << T.edges[e2](o2) << " " << T.edges[e3](o3) <<std::endl;
    if (!(T.edges[e1](o1) == t2 && T.edges[e2](o2) == t2 && T.edges[e3](o3) == t2)){return 1;}

    vertex p({0.99,0.99});
    int t_in = T.FindTriangleContainingPoint(0,p);
    std::cout << t_in << std::endl;
    if (!(t_in == T.triangles.size()-1)){return 1;}

    return 0;
}

int TEST3(){
    // Create a simple triangle
    std::vector<std::array<double, 2>> coords = {
        {0.0, 0.0}, {1.0, 0.0}, {0.5, 1.0}
    };

    std::vector<std::array<int, 3>> triangles = {
        {0, 1, 2}
    };

    triangulation T(triangles, coords);

    // Insert a point in the middle of the triangle
    vertex p({0.5, 0.5});
    T.vertices.push_back(p);
    T.InsertPoint(3);

    // Verify the new state after inserting the point
    // There should be 3 new triangles
    if (T.triangles.size() != 3) {
        std::cerr << "Number of triangles incorrect after inserting point." << std::endl;
        return 1;
    }

    // Verify the vertices of the new triangles
    std::vector<std::array<int, 3>> expected_triangles = {
        {0, 1, 3}, {1, 2, 3}, {2, 0, 3}
    };

    for (int i = 0; i < 3; ++i) {
        std::cout << T.triangles[i][0] << " " << T.triangles[i][1] << " " << T.triangles[i][2] << std::endl;
        for (int j = 0; j<3; j++){
            if (T.triangles[i][j] != expected_triangles[i][j]) {
                std::cerr << "Triangle vertices incorrect after inserting point." << std::endl;
                return 1;
            }
        }
    }

    // Verify the edges of the new triangles
    std::vector<std::array<int, 2>> expected_edges = {
        {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}
    };

    for (const auto& edge : expected_edges) {
        bool found = false;
        for (const auto& e : T.edges) {
            if ((e[0] == edge[0] && e[1] == edge[1]) || (e[0] == edge[1] && e[1] == edge[0])) {
                found = true;
                break;
            }
        }
        if (!found) {
            std::cerr << "Edge incorrect after inserting point." << std::endl;
            return 1;
        }
    }

    // Verify the edge orientations
    for (int i = 0; i < 3; ++i) {
        if (!(T.triangles[i].O(0) == 0 && T.triangles[i].O(1) == 0 && T.triangles[i].O(2) == 1)) {
            std::cerr << "Triangle orientations incorrect after inserting point." << std::endl;
            return 1;
        }
    }

    std::cout << "InsertPoint test passed." << std::endl;
    return 0;
}

int TEST4(){
    // Generate a set of points within the unit box
    std::vector<std::array<double, 2>> points;
    int num_points = 1000; // Number of points to generate
    std::srand(std::time(0)); // Seed for random number generation

    for (int i = 0; i < num_points; ++i) {
        double x = static_cast<double>(std::rand()) / RAND_MAX;
        double y = static_cast<double>(std::rand()) / RAND_MAX;
        points.push_back({x, y});
    }

    // Create a triangulation object
    triangulation T;

    // Perform Delaunay triangulation
    T.DelaunayTriangulation(points);
    T.writeVTK("TEST.vtk");

    // Verify the resulting triangulation is Delaunay
    bool pass = true;
    for (const auto& tri : T.triangles) {
        for (int i = 0; i < 3; ++i) {
            int e = tri(i);
            if (!T.EdgeIsDelaunay(e)) {
                std::cout << "Edge " << e << " is not Delaunay." << std::endl;
                pass = false;
                return 1;
            }
        }
    }

    if (pass){
        std::cout << "DelaunayTriangulation test passed." << std::endl;
    } else {
        std::cout << "DelaunayTriangulation test failed." << std::endl;
    }
    return 0;
}

int TEST5(){
    int npoints = 500;
    std::vector<std::array<double, 2>> points = generateRandomClosedShape(npoints, 1.0);

    // Create a triangulation object
    triangulation T;

    // Perform Delaunay triangulation
    T.DelaunayTriangulation(points);
    T.writeVTK("TEST_pre.vtk");
    std::cout << "Delaunay Triangulation done." << std::endl;
    double dx = points[1][0] - points[0][0];
    double dy = points[1][1] - points[0][1];
    double h_target =  2*PI/npoints;

    // perform Delaunay refinement
    T.DelaunayRefinement(h_target);
    T.writeVTK("TEST_post.vtk");
    std::cout << "Delaunay Refinement done." << std::endl;
    return 0;
}

int TEST6(){
    auto [points, segments] = readLakeSuperiorData("/Users/jacob/Documents/Me2sh/Tests/data/lake_superrior.dat");

    // Print points and segments for verification
    std::cout << "Points:" << std::endl;
    for (const auto& point : points) {
        std::cout << point[0] << ", " << point[1] << std::endl;
    }

    std::cout << "Segments:" << std::endl;
    for (const auto& segment : segments) {
        std::cout << segment[0] << ", " << segment[1] << std::endl;
    }

    // Create a triangulation object
    triangulation T;
    T.DelaunayTriangulation(points);
    T.writeVTK("TEST6_unconstrained.vtk");

    // Perform Delaunay triangulation
    T.DelaunayTriangluationConstrained(points, segments);
    T.writeVTK("TEST6_constrained.vtk");

    // Verify that the resulting triangulation has all segments in it
    for (const auto& segment : segments) {
        bool found = false;
        for (const auto& edge : T.edges) {
            if ((edge[0] == segment[0] && edge[1] == segment[1]) || (edge[0] == segment[1] && edge[1] == segment[0])) {
                found = true;
                break;
            }
        }
        if (!found) {
            std::cout << "Segment (" << segment[0] << "," << segment[1] << ") not found in triangulation." << std::endl;
            return 1;
        }
    }

    return 0;
}

int main(){
    int pass1 = TEST1();
    int pass2 = TEST2();
    int pass3 = TEST3();
    int pass4 = TEST4();
    int pass5 = TEST5();
    int pass6 = TEST6();

    return pass1 || pass2 || pass3 || pass4 || pass5 || pass6;
}