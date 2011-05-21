#ifndef __BEZIER2__
#define __BEZIER2__

#include "bezier.h" 
#include "psc.h"

// see curve.h for inherited class
class Bezier2 : public Bezier {
//Quadratic Bezier curve (1 control point)
public:
	
    Bezier2(PSC * c, int a, int b, double x, double y);
	
    void PreSplit(double tv);
	
    double findT(int p);
    Point findP(double t);
    Vector findTang(double t);
    void findAngle(int p, vector<Vector> & angs);
	
    Point p1;
	
};

#endif
