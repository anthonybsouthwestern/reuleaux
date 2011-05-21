#include "allh.h"

/*
 Uses the standard quadratic equation to
 locate roots of a cubic equation
 (pt^3+qt^2+rt+s=0) with one root known.
 
 Arguments below are:
 double p: the coeffiecent of t^3 in the equation.
 double q: the coeffiecent of t^2 in the equation.
 double r: the coeffiecent of t in the equation.
 double s: the constant value in the equation.
 double t: the known root to the equation.
 vector< double> & roots: the container for the roots to be located.
 */
bool Bezier3::QuadraticEquation(double p,  double q,  double r,  double s,  double t, vector< double> & roots)
{	
    double a,b,cx,i;
    double xreal,ximg;
	
    if (fabs(t) > c->config->epsilon) {
        a=p;
        b=q+(t*p);
        cx=(-s)/t;
    } else {
        a = p;
        b = q;
        cx = r;
    }
	
    i=((b*b)-(4*a*cx));
	
    if (i>=0 && fabs(a) > c->config->epsilon)
    {
        xreal =((-b)/(2*a));     
        ximg =(sqrt(i))/(2*a);
		
        roots.push_back(xreal+ximg);
        roots.push_back(xreal-ximg);
        return true;
    } else if (fabs(b) > c->config->epsilon) {
        roots.push_back(-b/cx);
    }
    return false;
}

/*
 Uses Newton's Method to approximate a root
 to the equation (at^3+bt^2+cx*t+d=0) and
 then uses QuadraticEquation to locate the
 rest of the roots.
 
 Arguments below are:
 double a: the coeffiecent of t^3 in the equation.
 double b: the coeffiecent of t^2 in the equation.
 double cx: the coeffiecent of t in the equation.
 double d: the constant value in the equation.
 vector <double> & roots: the container for the roots to be located.
 */
void Bezier3::CubicEquation(double a, double b, double cx, double d, vector <double> & roots) {
    int y;
	
    double x1,x2;
	
	
    if (fabs(d) > c->config->epsilon) {     
        // find a good starting point...
        x1 = 0.0;
        if (fabs((3*a*x1*x1)+(2*b*x1)+cx) < c->config->epsilon) {
            x1=1.0;
        } else if (fabs((3*a*x1*x1)+(2*b*x1)+cx) > c->config->epsilon){
            x1 = 0.5;
        }
        for(y=0;y<=100;y++)  {
            x2=x1-(((a*x1*x1*x1)+(b*x1*x1)+(cx*x1)+d)/((3*a*x1*x1)+(2*b*x1)+cx));
            x1=x2;
        }
    }
    else {
        x1 = 0.0;
    }
	
    roots.push_back(x1);
	
    QuadraticEquation(a,b,cx,d,x1,roots);
}

/*
 Arguments below are:
 PSC * c: the piecewise curved complex being meshed.
 int a: the left endpoint of the bezier curve.
 int b: the right endpoint of the bezier curve.
 Point p1: the first control point of the cubic bezier.
 Point p2: the second control point of the cubic bezier.
 */
Bezier3::Bezier3(PSC * c, int a, int b, Point p1, Point p2):Bezier(c) {
	
    type = CURVE_CUBBEZ;
    this->c = c;
    left[b] = a;
    right[b] = -1;
    left[a] = -1;
    right[a] = b;
	
    leftEndpoint = a;
    rightEndpoint = b;
    controlPts.push_back(p1);
    controlPts.push_back(p2);
	
    for (int jj=0; jj<2; jj++) {
        for (int ii=0; ii<2; ii++) {
            if (controlPts[jj].xy[ii] > c->tri->upperbound[ii])
                c->tri->upperbound[ii] = controlPts[jj].xy[ii];
            if (controlPts[jj].xy[ii] < c->tri->lowerbound[ii])
                c->tri->lowerbound[ii] = controlPts[jj].xy[ii];
        }
    }
}

/*
 Limits the total variation in orientation
 of the curve.
 Checks the angle between the tangent vectors
 of points on the bezier, adding points as
 necessary to keep the angle below tv.
 
 Arguments below are:
 double tv: the limit on total variation in orientation.
 */
void Bezier3::PreSplit(double tv) {
    Vector v1, v2, vmid;
    int p, newloc, leftloc = leftEndpoint;
    double tl=0.0, tr = 0.0, tmid = 0.0;
    double angle1 = 0.0, angle2 = 0.0;
	
    // handles update inside the loop
    for (p = right[leftEndpoint]; p != -1; 1) {
        tl = findT(leftloc);
        tr = findT(p);
        tmid = (tl+tr)/2.0;
		
        v1 = findTang(tl);
        v2 = findTang(tr);
        vmid = findTang(tmid); 
		
		
        angle1 = c->anglebtw1(v1, vmid);
        angle2 = c->anglebtw1(vmid, v2);
		
        //Check sum of angle between and the two points
        //to ensure total rotation of vector is obtained.
        if (angle1 + angle2 > tv) {
			
            // add point to Points data structure;
            newloc = c->pts->AddPoint(findP(tmid),POINT_PRESPLIT);
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
 Finds the t value for the parameteric equation
 of the bezier curve using CubicEquation.
 
 Arguments below are:
 int p: a point on the bezier curve.
 */
double Bezier3::findT(int p) {
	
    double a1,b1,c1,d1;
    // to solve:
    // a1*t^3 + b1*t^2 + c1*t + d1 = 0 
	
    //Use x-values to find a solution
    a1 = -1.0*c->pts->V[leftEndpoint].xy[0] + 3.0*controlPts[0].xy[0] 
    - 3.0*controlPts[1].xy[0] + c->pts->V[rightEndpoint].xy[0];
    b1 = 3.0*c->pts->V[leftEndpoint].xy[0] - 6.0*controlPts[0].xy[0] 
    + 3.0*controlPts[1].xy[0];
    c1 = -3.0*c->pts->V[leftEndpoint].xy[0] + 3.0*controlPts[0].xy[0];
    d1 = 1.0*c->pts->V[leftEndpoint].xy[0] - c->pts->V[p].xy[0];
	
    vector<double> roots;
    CubicEquation(a1,b1,c1,d1,roots);
	
    //Check x-values solutions with the y-values
    for (unsigned int i=0; i<roots.size(); i++) {
		
        a1 = -1.0*c->pts->V[leftEndpoint].xy[1] + 3.0*controlPts[0].xy[1] 
        - 3.0*controlPts[1].xy[1] + c->pts->V[rightEndpoint].xy[1];
        b1 = 3.0*c->pts->V[leftEndpoint].xy[1] - 6.0*controlPts[0].xy[1] 
        + 3.0*controlPts[1].xy[1];
        c1 = -3.0*c->pts->V[leftEndpoint].xy[1] + 3.0*controlPts[0].xy[1];
        d1 = 1.0*c->pts->V[leftEndpoint].xy[1] - c->pts->V[p].xy[1];
		
        double value = a1*roots[i]*roots[i]*roots[i] + b1*roots[i]*roots[i] +
        c1*roots[i] + d1;
        if (fabs(value)< 0.0000001) {
			return roots[i];
        }
    }
	
    //Use y-values to find a solution
    a1 = -1.0*c->pts->V[leftEndpoint].xy[1] + 3.0*controlPts[0].xy[1] 
    - 3.0*controlPts[1].xy[1] + c->pts->V[rightEndpoint].xy[1];
    b1 = 3.0*c->pts->V[leftEndpoint].xy[1] - 6.0*controlPts[0].xy[1] 
    + 3.0*controlPts[1].xy[1];
    c1 = -3.0*c->pts->V[leftEndpoint].xy[1] + 3.0*controlPts[0].xy[1];
    d1 = 1.0*c->pts->V[leftEndpoint].xy[1] - c->pts->V[p].xy[1];
	
    roots.clear();
    CubicEquation(a1,b1,c1,d1,roots);
	
    //Check y-values solutions with the y-values
    for (unsigned int i=0; i<roots.size(); i++) {
		
        a1 = -1.0*c->pts->V[leftEndpoint].xy[0] + 3.0*controlPts[0].xy[0] 
        - 3.0*controlPts[1].xy[0] + c->pts->V[rightEndpoint].xy[0];
        b1 = 3.0*c->pts->V[leftEndpoint].xy[0] - 6.0*controlPts[0].xy[0] 
        + 3.0*controlPts[1].xy[0];
        c1 = -3.0*c->pts->V[leftEndpoint].xy[0] + 3.0*controlPts[0].xy[0];
        d1 = 1.0*c->pts->V[leftEndpoint].xy[0] - c->pts->V[p].xy[0];
		
        double value = a1*roots[i]*roots[i]*roots[i] + b1*roots[i]*roots[i] +
        c1*roots[i] + d1;
        if (fabs(value)< 0.0000001) {
            return roots[i];
        }
    }
	
	
    cout << "WARNING: findT is having problems with " << name << endl;
    cout << roots.size() << endl;
    for (unsigned int i=0; i<roots.size(); i++) {
        cout << roots[i] << "," << a1*roots[i]*roots[i]*roots[i] + b1*roots[i]*roots[i] +
        c1*roots[i] + d1 << " :: ";
    } cout << endl;
    return -100.0;
}

/*
 Returns Point on bezier curve at t.
 
 double t: parametric value for the location of a point on the bezier curve.
 */
Point Bezier3::findP(double t) {
    return (1.0 - t) * (1.0 - t) * (1.0 - t) * c->pts->V[leftEndpoint]
    + 3.0 * t * (1.0 - t) * (1.0 - t) * controlPts[0]
    + 3 * t * t * (1.0 - t) * controlPts[1]
    + t * t * t * c->pts->V[rightEndpoint];
}

/*
 Returns the tangent vector at t.
 
 Arguments below are:
 double t: parametric value for the location of a point on the bezier curve
 */
Vector Bezier3::findTang(double t) {
    return 3.0 * (1.0 - t) * (1.0 - t) * (controlPts[0] - c->pts->V[leftEndpoint])
    + 6.0 * t * (1.0 - t) * (controlPts[1] - controlPts[0])
    + 3 * t * t * (c->pts->V[rightEndpoint] - controlPts[1]);
}

/*
 Find the tangent vector(s) in both directions
 of bezier curve, unless on an endpoint.
 
 Arguments below are:
 int p: the point on the bezier curve to locate tangent Vectors at.
 vector<Vector> angs: vector returned with the Vectors tangent at p.
 */
void Bezier3::findAngle(int p, vector<Vector> & angs) {
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