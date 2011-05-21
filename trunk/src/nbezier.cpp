#include "allh.h"

/*
 Arguments below are:
 PSC * c: the piecewise smooth complex being meshed.
 int a: the left endpoint of the nth-degree bezier curve.
 int b: the right endpoint.
 int degree: the degree of the bezier curve.
 */
NBezier::NBezier(PSC * c, int a, int b, int degree):Bezier(c) {
	
	type = CURVE_NBEZIER;
	this->c = c;
	left[b] = a;
	right[b] = -1;
	left[a] = -1;
	right[a] = b;
	this->degree = degree;
	
	leftEndpoint = a;
	pt2t[a] = 0.0;
	rightEndpoint = b;
	pt2t[b] = 1.0;
	
	controlPts.resize(degree + 1);
	controlPts[0] = c->pts->V[a];
	controlPts[degree] = c->pts->V[b];
}

/*
 Add control point to bezier.
 
 Arguments below are:
 Point pt: the control point being added to the bezier curve.
 int pos: the number of the of the control point.
 controlPts[0] is the left endpoint and
 controlPts[degree] is the right endpoint.
 */
void NBezier::addPt(Point & pt, int pos) {
    controlPts[pos] = pt;
	
    for (int ii=0; ii<2; ii++) {
        if (pt.xy[ii] > c->tri->upperbound[ii])
            c->tri->upperbound[ii] = pt.xy[ii];
        if (pt.xy[ii] < c->tri->lowerbound[ii])
            c->tri->lowerbound[ii] = pt.xy[ii];
    }
}

/*
 Checks the angle between the tangent vectors of
 points on the beier, adding points as necessary
 to keep the angle below tv.
 
 Arguments below are:
 double tv: the total variation in orientation in radians.
 */
void NBezier::PreSplit(double tv) {
    int p, newloc, leftloc = leftEndpoint;
    double tl=0.0, tr = 0.0, tmid = 0.0;
    double angle;
    Vector leftTan, rightTan;
	
    // handles update inside the loop
    for (p = right[leftEndpoint]; p != -1; 1) {
        tl = findT(leftloc);
        tr = findT(p);
        tmid = (tr-tl)/(degree - 1);
		
        angle = 0;
        leftTan = findTang(tl);
		
        for (int i = 1; i < (degree - 1); i++) {
            rightTan = findTang(tl + tmid * i);
            angle += c->anglebtw1(leftTan, rightTan);
            leftTan = rightTan;
        }
        rightTan = findTang(tr);
        angle += c->anglebtw1(leftTan, rightTan);
		
		
        //Check some of angle between and the interior points to ensure total rotation of vector is obtained.
        if (angle > tv) {
            tmid = (tr+tl)/2.0;
            // add point to Points data structure;
            newloc = c->pts->AddPoint(findP(tmid),POINT_PRESPLIT);
            pt2t[newloc] = tmid;
            left[newloc] = leftloc;
            right[leftloc] = newloc;
			
            // add to list to be inserted
            c->presplitVertices.push_back(newloc);
			
            // set l2s (location2segment)
            c->pts->l2s[newloc].insert(this);
			
            right[newloc] = p;
            left[p] = newloc;
			
            //Adjust bounding box
            for (int ii = 0; ii < 2; ii++) {
                if (c->pts->V[newloc].xy[ii] > c->tri->upperbound[ii])
                    c->tri->upperbound[ii] = c->pts->V[newloc].xy[ii];
                if (c->pts->V[newloc].xy[ii] < c->tri->lowerbound[ii])
                    c->tri->lowerbound[ii] = c->pts->V[newloc].xy[ii];
            }
			
            p = newloc;
        }
        else {
            leftloc = p;
            tl = tr;
            p = right[leftloc];
        }
    }
}


/*
 Finds the parametric t-value for a point on
 the Bezier curve.
 
 Arguments below are:
 int p: the point on the bezier curve.
 */
double NBezier::findT(int p) {
	
    if (pt2t.find(p) != pt2t.end()) { //If the tvalue is already known
        return pt2t[p];
    }
    for (vector<double>::iterator it = splitAtDistTvals.begin();
         it != splitAtDistTvals.end();it++) {
        if (findP(*it) == c->pts->V[p]) {
            pt2t[p] = *it;
            splitAtDistTvals.erase(it);
            return pt2t[p];
        }
    }
    //Binary search to locate the t-value.
    double tl = 0.0, tr, tin, ldist, rdist, tpoint;
	
    tin = 1.0/((double)degree);
    ldist = d(c->pts->V[leftEndpoint], c->pts->V[p]);
    rdist = d(c->pts->V[rightEndpoint], c->pts->V[p]);
	
    for (int i = 1; i <= degree; i++) {
    tr = tin * (double) i;
        if (fabs(tr - 1) < c->config->epsilon)
            tr = 1.0;
        tpoint = findT(p, tl, tr);
        if(tpoint > 0) {
            pt2t[p] = tpoint;
            return tpoint;
        }
        tl = tin * (double) i;
    }
	
    cout << "Error in findT at point " << p << " " << c->pts->V[p] << endl;
    return 0;
}

/*
 Recursive findT
 
 Arguments below are:
 int p: the point on the bezier curve.
 double tlower: the lower bound for the binary search.
 double tupper: the upper bound for the binary search.
 */
double NBezier::findT(int p, double tlower, double tupper) {
    if (tupper - tlower < c->config->epsilon) {
        return -1;
    }
    double tmid = (tlower + tupper)/2.0, dmid = d(findP(tmid), c->pts->V[p]);
	
    if(fabs(dmid) < c->config->epsilon)
        return tmid;
	
    double dlow = d(findP(tlower), c->pts->V[p]), dup = d(findP(tupper), c->pts->V[p]);
    if (dlow < dup) {
        return findT(p, tlower, tmid);
    } else {
        return findT(p, tmid, tupper);
    }
}

/*
 Find midpoint between l and right of l
 */
Point NBezier::MidPoint(int l) {
    double tmid;
	
    double tl = findT(l);
    double tr = findT(right[l]);
	
    tmid = (tl + tr) / 2.0;
    splitAtDistTvals.push_back(tmid);
	
    return findP(tmid);
}

/*
 Returns point on bezier curve at t.
 
 Arguments below are:
 double t: the parametric t-value for the point returned.
 */
Point NBezier::findP(double t) {
    Point p(0,0);
    for (int i = 0; i <= degree; i++) {
        p = p + Binom(degree, i) * pow(1 - t, degree - i) * pow(t, i) * controlPts[i];
    }
    return p;
}

/*
 Calculates the binomial coefficient using
 dynamic programming. i.e. n choose m
 
 Arguments below are:
 int n: the total number of choices
 int m: the number being chosen.
 */
int NBezier::Binom(int n, int m)
{
    vector<int> b (n + 1);
    b [0] = 1;
    for (int i = 1; i <= n; ++i)
    {
        b [i] = 1;
    for (int j = i - 1U; j > 0; --j)
        b [j] += b [j - 1U];
    }
    return b [m];
}

/*
 Returns the tangent vector at t.
 
 Arguments below are:
 double t: the parametric t-value at the point of the tangent vector being found.
 */
Vector NBezier::findTang(double t) {
    Vector tangent(0, 0);
    Point p(0,0);
    if (t == 0) {
        tangent = degree*(controlPts[1] - controlPts[0]);
        return tangent;
    }
    if (t == 1) {
        tangent = degree*(controlPts[degree] - controlPts[degree-1]);
        return tangent;
    }
    for (int i = 0; i <= degree; i++) {
        if (i == degree) {
            tangent = tangent + p - Binom(degree, i) * controlPts[i] * (- i * pow(t, i - 1));
        } else {
            tangent = tangent + p - Binom(degree, i) * controlPts[i] * 
            ((degree - i) * pow(1 - t, degree - i - 1) * pow(t, i) - i * pow(t, i - 1) * pow(1 - t, degree - i));
        }
    }
    return tangent;
}

/*
 Find the tangent vector(s) in both directions of
 the bezier curve, unless on an endpoint.
 
 Arguments below are:
 int p: the point at which the tangent vectors are being found.
 vector<Vector> angs: where the resulting tangent vectors are stored and returned.
 */
void NBezier::findAngle(int p, vector<Vector> & angs) {
    if (p == leftEndpoint) {
        angs.push_back(findTang(0.0));
    } else if (p == rightEndpoint) {
        // tangent vector points towards p1 in both cases...
        angs.push_back(-1.0 * findTang(1.0));
    } else {
        double t = findT(p);
        angs.push_back(findTang(t));
        angs.push_back(- angs[0]);
    }
}

/*
 Used for the SVG output of the nth degree bezier curve. Returns the
 endpoints of straight line segments to be used to approximate the
 curve.
 
 Arguments below are:
 int n: the number of linesegments to use in the approximation.
 vector<Point> pts: the vector of points for the line segment approximation.
    The first point is the left endpoint of the first linesegment, the second
    point is that linesegments right endpoint as well as the next linesegment's
    leftendpoint.
 */
void NBezier::getOutput(int n, vector<Point> & pts) {
    double tin = 1;
    tin/=(double)n;
	
    pts.push_back(controlPts[0]);
    for (double i = tin; i < 1; i+=tin) {
        pts.push_back(findP(i));
    }
    pts.push_back(controlPts[degree]);
}

/*
 Used for the SVG output of the nth degree bezier curve. Creates quadratic bezier
 curves as approximations of the nth degree bezier curve between points on the curve.
 Finds control points from the intersection of tangent lines between endpoints of the
 quadratic bezier approximations.
 
 Arguments below are:
 vector<Point> pts: the points and control points of the quadratic beziers created.
 The first point is the left endpoint, then the first quadratic bezier's control
 point, then the right endpoint which is also the left endpoint of the next bezier
 curve.
 */
void NBezier::get2Beziers(vector<Point> & pts) {
    int lp = leftEndpoint, rp;
    Vector rtang, ltang = findTang(findT(lp));
    double paraVar; //Will store the t parametric value for the intersection of tangent lines.
    Point lcoord, rcoord, p1;
	
    lcoord = c->pts->V[lp];
    pts.push_back(lcoord);
	
    for (rp = right[leftEndpoint]; rp != -1; rp = right[lp]) {
        rcoord = c->pts->V[rp];
        rtang = findTang(findT(rp));
        paraVar = (rtang.xy[0] * (lcoord.xy[1] - rcoord.xy[1]) + rtang.xy[1] * (rcoord.xy[0] - lcoord.xy[0]))
        / (rtang.xy[1] * ltang.xy[0] - ltang.xy[1] * rtang.xy[0]);
		
        if (fabs(c->anglebtw(rtang, ltang)) < PI/8) {
            p1 = (lcoord + rcoord) / 2;
        } else {
            p1 = lcoord + ltang * paraVar;
        }
        pts.push_back(p1);
        pts.push_back(rcoord);
		
        lp = rp;
        ltang = rtang;
        lcoord = rcoord;
    }
}

/*
 This function will be useful for handling small angles. 
 
 Arguments below are:
 int l: location of the input vertex to locate points around
 double r: radius of the circle to make
 vector<Point> pts: vector returning 1 to 2 points
 */
void NBezier::SplitAtDist(int l, double r, vector<Point> & pts) {
    double tvertex = findT(l);
    //If vertex is not rightEndpoint
    if (tvertex < (1 - c->config->epsilon)) {
        double tupper = 1.0;
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
        splitAtDistTvals.push_back(tcurrent);
    }
	
    //If vertex is not leftEndpoint
    if (tvertex > c->config->epsilon) {
        double tupper = 0.0;
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
        splitAtDistTvals.push_back(tcurrent);
    }
}
