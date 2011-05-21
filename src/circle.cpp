#include "allh.h"


Circle::Circle(PSC * c):Curve(c) {}

/*
 Arguments below are:
 Psc * c: the piecewise smooth complex being meshed.
 Point & a: the centerpoint of the circle.
 double b: the radius of the circle.
 */
Circle::Circle(PSC * c, Point & a, double b):Curve(c) {
	
    type = CURVE_CIRCLE;
    this->c = c;
    center = a;
    radius = b;
}


/*
 Computes the midpoint of l and right[l].
 */
Point Circle::MidPoint(int l) {
	
    double theta = anglebtw(c->pts->V[l], c->pts->V[right[l]]) / 2;
    return rotatePoint(l, theta);
}

/*
 The default presplit of circles is 90 degrees.
 Adds in points on the circe so that no two neighboring points
 are more than 90 degrees apart.
 */
void Circle::PreSplit() {
    PreSplit(PI/2);
}

/*
 Adds in points on the circe so that no
 two neighboring points are more than tv
 radians apart.
 */
void Circle::PreSplit(double tv) {
    int loc;
	
    if (left.size() == 0) {		//No points on the circle
        // compute the first point
        Vector v1(0.0,radius);
        Point p1;
        p1 = center + v1;
		
        //Adjust bounding box
        for (int i = 0; i < 2; i++) {
            if (p1.xy[i] > c->tri->upperbound[i]) {
                c->tri->upperbound[i] = p1.xy[i];
            }
            if (p1.xy[i] < c->tri->lowerbound[i]) {
                c->tri->lowerbound[i] = p1.xy[i];
            }
        }
		
        // add point to Points data structure;
        loc = c->pts->AddPoint(p1,POINT_PRESPLIT);
        right[loc] = loc;
        left[loc] = loc;
		
        // add to list to be inserted
        c->presplitVertices.push_back(loc);
		
        // set leftEndpoint and rightEndpoint
        leftEndpoint = loc;
        rightEndpoint = loc;
		
        // set l2s (location2segment)
        c->pts->l2s[loc].insert(this);
    }
    double theta;
	
    int current = leftEndpoint;
    do {
        int r = right[current];
        theta = anglebtw(c->pts->V[current], c->pts->V[r]); //The arc length between the current point and its right.
		
        int numSplits = ceil(theta / tv);
        theta/= numSplits;
        //Adds points in to split the circle segment into numSplits segments
        for (int i = 1; i < numSplits; i++) {
			
            Point p = rotatePoint(current, theta);//Compute next point
			
            //Adjust bounding box
            for (int ii = 0; ii < 2; ii++) {
                if (p.xy[ii] > c->tri->upperbound[ii])
                    c->tri->upperbound[ii] = p.xy[ii];
                if (p.xy[ii] < c->tri->lowerbound[ii])
                    c->tri->lowerbound[ii] = p.xy[ii];
            }
			
            //Add new point to points data structure
            loc = c->pts->AddPoint(p,POINT_PRESPLIT);
			
            right[current] = loc;
            left[loc] = current;
			
            // add to list to be inserted
            c->presplitVertices.push_back(loc);
            // set l2s (location2segment)
            c->pts->l2s[loc].insert(this);
			
            current = loc;
        }
        right[current] = r;
        left[r] = current;
        current = r;
		
    }while(current != rightEndpoint);
	
}

/*
 Auxiliary PreSplit used for adding small angle
 circles. Since these circles always have at
 least 2 points on them, it is not necessary to
 check if no points are already present. Returns
 the keys of the new points added.
 */
void Circle::PreSplit(vector<int> &newPts) {
    PreSplit(newPts, PI/2);
}
void Circle::PreSplit(vector<int> &newPts, double tv) {
    int loc;
    double theta;
    int current = leftEndpoint;
	
    do {
		
        int r = right[current];
        theta = anglebtw(c->pts->V[current], c->pts->V[r]);
		
        int numSplits = ceil(theta / tv);
		theta/= numSplits;
		
        for (int i = 1; i < numSplits; i++) {
            Point p = rotatePoint(current, theta);//Compute next point
			
            //Adjust bounding box
            for (int ii = 0; ii < 2; ii++) {
                if (p.xy[ii] > c->tri->upperbound[ii])
                    c->tri->upperbound[ii] = p.xy[ii];
                if (p.xy[ii] < c->tri->lowerbound[ii])
                    c->tri->lowerbound[ii] = p.xy[ii];
            }
			
            //Add new point to points data structure
            loc = c->pts->AddPoint(p,POINT_PRESPLIT);
            newPts.push_back(loc);
            current = loc;
        }
        current = r;
		
    }while(current != rightEndpoint);
}

/*
 Returns the arc length in radians between
 point x on the circle and point y on the circle.
 
 Arguments below are:
 Point x: the right point.
 Point y: the left point. (clockwise rotation)
 */
double Circle::anglebtw(Point x, Point y) {
	
    if(x == y)		//360 degree angle.
        return 2*PI;
	
    Vector l = x - center;
    Vector r = y - center;
    l.normalize();
    r.normalize();
	
    double mysin = r.xy[0] * l.xy[1] - r.xy[1] * l.xy[0];
    // After normalizing, the product can in fact be slightly larger than 1, just due to precision issues.
    // The check below ensures that the asin function will succeed. 
    // Similarly for arccos, acos.
	
    if (mysin > 1.0) mysin = 1.0;
    if (mysin < -1.0) mysin = -1.0;
    double sinangle = asin(mysin);  //Arcsin of magnitude of unit vector cross product,
	
    double mycos = l*r;
    if (mycos > 1.0) mycos = 1.0;
    if (mycos < -1.0) mycos = -1.0;
    //Arccos of the unit vector dot product, both of which compute the angle between the two vectors.
    double cosangle = acos(mycos);
	
    // 1st quaadrant
    if (sinangle > 0.0 && cosangle <PI/2.0) {
        return cosangle;
    }
    // 2nd quadrant
    if (sinangle > 0.0) {
        return cosangle;
    }
    // 3rd quadrant
    if (cosangle >= PI/2.0) {
        return PI - sinangle;
    }
    // 4th quadrant
    return 2*PI + sinangle;
}

/*
 Returns a point on the circle rotated clockwise 
 by the amount ang given in radians.
 
 Arguments below are:
 int p: the starting point.
 double ang: the amount in radians to rotate point on circle.
 */
Point Circle::rotatePoint(int p, double ang){
    double x = c->pts->V[p].xy[0]-center.xy[0]; // Rotate around center
    double y = c->pts->V[p].xy[1]-center.xy[1];
    double costheta = cos(ang);
    double sintheta = sin(-ang);
	
    //Get rid of rounding errors with zero
    if (fabs(costheta) < c->config->epsilon)
        costheta = 0;
    if (fabs(sintheta) < c->config->epsilon)
        sintheta = 0;
	
    Vector pt(x * costheta - y * sintheta, x * sintheta + y * costheta);
    if (fabs(pt.xy[0]) < c->config->epsilon)
        pt.xy[0] = 0;
    if (fabs(pt.xy[1]) < c->config->epsilon)
        pt.xy[1] = 0;
    return center+pt; // shift based on circle center
}

/*
 This function returns the locations of points on the circle
 a distance r from the given point l.
 
 Arguments below are:
 int l: location of the input vertex to locate points around
 double r: radius of the circle to make
 vector<Point> pts: vector returning 1 or 2 points
 */
void Circle::SplitAtDist(int l, double r, vector<Point> & pts) {
    int nb1 = left[l];
    int nb2 = right[l];
	
	
    if (d2(c->pts->V[nb1],c->pts->V[l]) < r*r || d2(c->pts->V[nb2],c->pts->V[l]) < r*r) {
        cout << "ERROR: SplitAtDist given an overly long length" << endl;
    }
    double theta = 2 * asin((r / 2) / radius);
	
    Point p1 = rotatePoint(l, theta);
    pts.push_back(p1);
	
    Point p2 = rotatePoint(l, -theta);
    pts.push_back(p2);
}

/*
 Adds point p to circle construct. First locates the 
 left and right of p and places it in the circle.
 
 Arguments below are:
 int p: the point to add to the circle.
 */
void Circle::InsertVertex(int p) {
	
    int found = -1;
	
    if (left.size() == 0) {	//No points already on the circle
        leftEndpoint = p;
        rightEndpoint = p;
        left[p] = p;
        right[p] = p;
    } else if (right[leftEndpoint] == rightEndpoint) {	//Either one point on the circle, or two points if it is an arc.
        left[rightEndpoint] = p;
        right[leftEndpoint] = p;
        left[p] = leftEndpoint;
        right[p] = rightEndpoint;
    } else {	//Multiple points already on circle, must locate the left of point to be inserted
        double theta = anglebtw(c->pts->V[leftEndpoint], c->pts->V[p]);	//Position of point to be inserted.
		
        int loc = leftEndpoint;
        do {	//Check for loc to be left of p, and right[loc] to the right
            if ( ( anglebtw(c->pts->V[leftEndpoint], c->pts->V[loc]) < theta ||
                  anglebtw(c->pts->V[leftEndpoint], c->pts->V[loc]) > 2*PI - c->config->epsilon) && 
                (anglebtw(c->pts->V[leftEndpoint], c->pts->V[right[loc]]) > theta ||
                 anglebtw(c->pts->V[leftEndpoint], c->pts->V[right[loc]]) < c->config->epsilon) ) {
					
                    //Add point to circle
                    int r = right[loc];
                    right[loc] = p;
                    right[p] = r;
                    left[r] = p;
                    left[p] = loc;
                    found = loc;
					
                    break;
                } 
            loc = right[loc];
        } while (loc != rightEndpoint);
    }
    //Tell points complex that p is on the circle
    c->pts->l2s[p].insert(this);
}

/*
 Insert point p to the right of l into circle.
 
 Arguments below are:
 int p: the vertex to be added to the circle.
 int l: the point to the left of p on the circle.
 */
void Circle::InsertVertex(int p, int l) {
	
    if (left.size() < 1) {
        InsertVertex(p);
    }
	
    int r = right[l];
    left[r] = p;
    right[l] = p;
    left[p] = l;
    right[p] = r;
	
    c->pts->l2s[p].insert(this);
}

