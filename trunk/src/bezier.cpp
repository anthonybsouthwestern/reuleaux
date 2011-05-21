#include "allh.h"

Bezier::Bezier(PSC * c):Curve(c){}


/*
 Find midpoint between l and right of l
 */
Point Bezier::MidPoint(int l) {
	double tmid;
	
	double tl = findT(l);
	double tr = findT(right[l]);
	
	tmid = (tl + tr) / 2.0;
	
	return findP(tmid);
}

/*
 This function will be useful for handling small angles. 
 
 Arguments below are:
 int l: location of the input vertex to locate points around
 double r: radius of the circle to make
 vector<Point> pts: vector returning 1 to 2 points
 */
void Bezier::SplitAtDist(int l, double r, vector<Point> & pts) {
	double tvertex = findT(l);
	
	if (tvertex < (1 - c->config->epsilon)) { //Vertex is not the rightEndpoint
		double tupper = findT(right[l]);
		double tlower = tvertex;
		double tcurrent = (tupper + tlower)/2.0;
		Vector v = findP(tcurrent) - c->pts->V[l];
		
		while (fabs(v.len() - r) > c->config->epsilon) {
			if (v.len() > r) {
				tupper = tcurrent;
			} else {
				tlower = tcurrent;
			}
			tcurrent = (tupper + tlower)/2.0;
			v = findP(tcurrent) - c->pts->V[l];
		}
		pts.push_back(findP(tcurrent));
	}
	
	
	if (tvertex > (c->config->epsilon)) {	//vertex is not the leftEndpoint
		double tupper = findT(left[l]);
		double tlower = tvertex;
		double tcurrent = (tupper + tlower)/2.0;
		Vector v = findP(tcurrent) - c->pts->V[l];
		
		while (fabs(v.len() - r) > c->config->epsilon) { 
			if (v.len() > r) {
				tupper = tcurrent;
			} else {
				tlower = tcurrent;
			}
			tcurrent = (tupper + tlower)/2.0;
			v = findP(tcurrent) - c->pts->V[l];
		}
		pts.push_back(findP(tcurrent));
	}
}

/*
 Insert point p into bezier if location on bezier is not known.
 */
void Bezier::InsertVertex(int p) {
	int l = leftEndpoint;
	int r = right[l];
	
	if (r != rightEndpoint) {
		for (double t = findT(p); findT(r) < t; r = right[l]) { l = r; }
	}
	
	right[l] = p;
	right[p] = r;
	left[r] = p;
	left[p] = l;
	
	c->pts->l2s[p].insert(this);
}
/*
 Adds points to bezier curve so that total variation
 (angle between tangent vectors) of two neighboring
 points is no more than 45 degrees
 */
void Bezier::PreSplit() {
    PreSplit(PI/4.0);
    // presplit needs to be at 45 degrees for bezier curves
    // rather than the 90 for arcs/circles
    // generally true for ellipses as well -- anything non-circular
    // related to total variation
}

/*
 Returns the smallest distance from a given
 point to the nearest point on the bezier as
 determined using straight line approximations
 between points on the bezier.
 
 Arguments below are:
 int v: the vertex
 */
double Bezier::distToNearestPt(int v) {
    int l = leftEndpoint;
    double radius = INFINITY, theta;
    Vector v1, v2;
	
    for(int r = right[l]; l != rightEndpoint;r = right[l]) {
        v1 = c->pts->V[v] - c->pts->V[r];
        if (fabs(v1.len()) < radius) {
            radius = fabs(v1.len());
        }		
        v1 = c->pts->V[v] - c->pts->V[l];
        if (fabs(v1.len()) < radius) {
            radius = fabs(v1.len());
        }
        v2 = c->pts->V[r] - c->pts->V[l];
        theta = (PI/2) - min(c->anglebtw(v2, v1), c->anglebtw(v1, v2));
        if (theta > 0) {
            if (fabs(v1.len() * cos(theta)) < radius) {
                radius = fabs(v1.len() * cos(theta));
            }
        }
        l = r;
    }
    return radius;
}