#ifndef __BEZIER3__
#define __BEZIER3__

#include "bezier.h" 
#include "psc.h"

// see curve.h for inherited class
class Bezier3 : public Bezier {
//Cubic bezier curve class (2 control points)
public:
	
    Bezier3(PSC * c, int a, int b, Point p1, Point p2);
	
    void PreSplit(double tv);
	
    double findT(int p);
    double findT(int p, double tlower, double tupper);
    Point findP(double t);
    Vector findTang(double t);
    void findAngle(int p, vector<Vector> & angs);
	
    vector <Point> controlPts;

    bool QuadraticEquation(  double p,  double q,  double r,  double s,  double t, vector< double> & roots);
    void CubicEquation(double a, double b, double c, double d, vector <double> & roots);
};

#endif
