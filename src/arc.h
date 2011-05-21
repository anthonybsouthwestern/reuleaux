#ifndef __ARC__
#define __ARC__

#include "circle.h" 
#include "psc.h"

// see circle.h for inherited class
class Arc : public Circle {
//Circular arcs
public:
	
    Arc(PSC * c, Point & a, int le, int re);
    bool isOn(Point & p);
    void SplitAtDist(int l, double r, vector<Point> & pts);
};

#endif
