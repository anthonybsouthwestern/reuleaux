#ifndef __BEZIER__
#define __BEZIER__

#include "curve.h" 
#include "psc.h"

// see curve.h for inherited class
class Bezier : public Curve {
//Bezier curve superclass with some common functions
public:
	
    Bezier(PSC * c);
	
    Point MidPoint(int l);
    void SplitAtDist(int l, double r, vector<Point> & pts);
    void InsertVertex(int p);
    void PreSplit();
    double distToNearestPt(int v);
    virtual double findT(int p) {double t; return t;};
    virtual Point findP(double t) {Point p; return p;};
    virtual void PreSplit(double tv) {};

};

#endif
