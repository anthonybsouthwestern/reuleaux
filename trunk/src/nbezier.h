#ifndef __NBEZIER__
#define __NBEZIER__

#include "curve.h" 
#include "psc.h"

// see curve.h for inherited class
class NBezier : public Bezier {
//nth degree bezier curve class
public:
	
    NBezier(PSC * c, int a, int b, int degree);
	
    void PreSplit(double tv);
	
    double findT(int p);
    double findT(int p, double tlower, double tupper);
    Point findP(double t);
    int Binom(int n, int m);
    Vector findTang(double t);
    void findAngle(int p, vector<Vector> & angs);
    void addPt(Point & pt, int pos);
    void getOutput(int n, vector<Point> & pts);
    void get2Beziers(vector<Point> & pts);
    Point MidPoint(int l);
    //Control Points  include endpoints
    vector <Point> controlPts;
    int degree;
    map <int, double> pt2t;
    vector <double> splitAtDistTvals;

    void SplitAtDist(int l, double r, vector<Point> & pts);
	
};

#endif
