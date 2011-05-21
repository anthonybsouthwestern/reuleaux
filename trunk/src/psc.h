#ifndef __PSC__
#define __PSC__

#include "common.h"

using namespace std;

 
struct CurveAngle {
  double angle;
  Vector dir;
  Curve * c;
  double curvatureInfo;
};

struct SmallAnglePt {
    int pt;
    double dist;
};


//
// psc.h
// The Piecewise Smooth Complex class
//
// TODO:
//  * Support more output formats including:
//     - .node/.ele used by Triangle
//
//
class PSC {
public:
	
    string name;
	
    Points *pts;      // pointer to points class
    Geometry *geom;    // pointer to geometry routines
    TrianglePool *triPool; // memory allocator for triangles
	
    Triangulation * tri;  // the Delaunay triangulation
    vector <Curve *> curves;  // psc input curves with their 
    // subdivisions into subcurves
    vector <int> inputVertices; // list of input vertices
    vector <int> presplitVertices; // list of presplit vertices
    vector <int> bboxVertices;

    set<int> protectedVertices;

    set<Triangle *> exteriorTriangles;
    set<Triangle *> smallAngleTriangles;

    Configuration *config;
    Stats *stats;
	
    map<int, double> vertex2Radius;
	
    // constructor...
    PSC();           // plain constructor
	
    void OpenLogFile(string filename); // open .log file for output
	
    void ReadCfg(string filename);  // read options from the .cfg file
	
    void ReadPsc(string filename); // read input .psc file
    void ReadPlc(string filename); //read input .plc file

    void PreSplit(); // split all arcs to 90 degrees
	
    void TriangulateInputVertices();  // compute initial Delaunay triangulation
    void RuppertsAlgorithm();         // run Ruppert's algorithm
	
    bool checkSteinerPoint(Point & p);
    bool checkSteinerPoint(Point & p, Triangle * t);
    bool checkSteinerPoint(Point & p, vector<Triangle *> tris);
	
    void PrintInput(string filename); // print the input PSC
    // add function to print PSC in svg format
    void PrintMesh(string filename); // print .eps file
    void PrintSVG(string filename, bool printTriangulation,
                  bool noTriangulationWithArcs, bool printVertices); // print svg image
    void PrintSVG(string filename, bool printTriangulation,
                  bool noTriangulationWithArcs, bool printVertices,
                  Point center, double size); // print svg image
    void PrintNodeEle(string filename);
    void InitBoundingBox();
	
    void SetUpExteriorTriangles();
    void SetUpSmallAngleTriangles();

    void LogResults();
	
    void FixSmallAngles();
    void addSmallAngleCircle(int vertex, double radius);
    double anglebtw(Vector l, Vector r);
    double anglebtw1(Vector l, Vector r);
	
    void updateTriangulation(int p, Circle * newCir);
    
    ofstream *LOGFILE;
	
    bool logfileisopen;
private:
	
    // check encroachment for ruppert's algorithm:
    bool IsEncroached(Subsegment s);
    // DoesEncroach: does point p encroach a subcurve with an endpoint
    //               in nbs (which is a vertex of t)?
    bool DoesEncroach(Point & p, set<int> nbs, set<Subsegment> & s); 
    bool DoesEncroach(Point & p, Triangle * t, set<Subsegment> & s);
    bool DoesEncroach(Point & p, vector<int> nbs, vector<Subsegment> & s); 
    bool DoesEncroach(Point & p, Triangle * t, vector<Subsegment> & s);
	
    bool IsProtected(int p);
    bool IsProtected(Subsegment ss);
    bool IsProtected(Triangle * t);

    bool IsArc(int a, int b);

    bool IsExterior(Triangle * t);
    bool IsInterior(Triangle * t);
	
    bool IsTooClose(Point &p);

    // useful helper functions
    void SetRadiusAndCC(Triangle * s);  
    void PrintBegin(string s);
    void PrintEnd(string s);

    // The STL set class is used as a priority queue...
    // this works since the stl set stores items in a search tree
    // and iterates through the tree using an in order traversal.
    // The segmentQueue.begin() points to the highest priority item.
    set<Subsegment> segmentQueue;
    set<Triangle *> triangleQueue;
};

#endif              // !defined(__PSC__)
