
#include "allh.h"


/*
 Arguments below are:
 Point a: the center point of the circle that describes the arc.
 int le: the left endpoint of the arc
 int re: the right endpoint.
 */
Arc::Arc(PSC * c, Point & a, int le, int re):Circle(c) {
    Vector v = c->pts->V[le] - a;
    type = CURVE_ARC;
    this->c = c;
    center = a;
    radius = v.len();
	
    left[le] = -1;
    right[le] = re;
    left[re] = le;
    right[re] = -1;
	
    leftEndpoint = le;
    rightEndpoint = re;
}

/*
 Check if point on circle describing arc
 is on the section described by the arc.
 
 Arguments below are:
 p: the point on the arc
 */
bool Arc::isOn(Point & p) {
    return anglebtw(c->pts->V[leftEndpoint], p) <= anglebtw(c->pts->V[leftEndpoint], c->pts->V[rightEndpoint]);
}

/*
 This function is used for handling small angles. 
 
 Arguments below are:
 int l: location of the input vertex to locate points around
 double r: radius of the circle to make
 vector<Point> pts: vector returning 1 or 2 points
 */
void Arc::SplitAtDist(int l, double r, vector<Point> & pts) {
    int nb1 = left[l];
    int nb2 = right[l];
    double theta = 2 * asin((r / 2) / radius);
	
    if (nb1 != -1) {
        if (d2(c->pts->V[nb1],c->pts->V[l]) < r*r) {
            cout << "ERROR: SplitAtDist given an overly long length" << endl;
        }
		
        Point p1 = rotatePoint(l, -theta);
        pts.push_back(p1);
    }
	
    if (nb2 != -1) {
        if (d2(c->pts->V[nb2],c->pts->V[l]) < r*r) {
            cout << "ERROR: SplitAtDist given an overly long length" << endl;
        }
		
        Point p2 = rotatePoint(l, theta);
        pts.push_back(p2);
    }
}