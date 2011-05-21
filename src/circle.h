#ifndef __CIRCLE__
#define __CIRCLE__

#include "curve.h" 
#include "psc.h"

// see curve.h for inherited class
class Circle : public Curve {
//Circles
public:
	
    Circle(PSC * c);
    Circle(PSC * c, Point & a, double b);
	
    //The point midway between l and right[l] on the circle
    Point MidPoint(int l);
	
    //Add points to circle to reduce the variation in orientation.
    void PreSplit();
    void PreSplit(double tv);
    void PreSplit(vector<int> &newPts);
    void PreSplit(vector<int> &newPts, double tv);
	
    //Find the points on the circle a distance r from the point l
    void SplitAtDist(int l, double r, vector<Point> & pts);
	
    //Arclength in radians between x any y
    double anglebtw(Point x, Point y);
    //Clockwise rotation distance ang in radians
    Point rotatePoint(int p, double ang);	
	
    //Add point p to circle, first must locate point's left and right
    void InsertVertex(int p);
    //Insert p to the right of l
    void InsertVertex(int p, int l);	
	
    double radius;
    Point center;
};

#endif
