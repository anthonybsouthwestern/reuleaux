#include "allh.h"

/*
 Arguments below are:
 PSC * c: the piecewise smooth complex being meshed.
 int a: the left endpoint of the quadratic bezier curve.
 int b: the right endpoint.
 double x: the x-coordinate of the control point.
 double y: the y-coordinate of the control point.
 */
Bezier2::Bezier2(PSC * c, int a, int b, double x, double y):Bezier(c) {
	
    type = CURVE_BEZIER;
    this->c = c;
    left[b] = a;
    right[b] = -1;
    left[a] = -1;
    right[a] = b;
	
    leftEndpoint = a;
    rightEndpoint = b;
    p1 = Point(x, y);
}

/*
 Checks the angle between the tangent vectors
 of points on the bezier, adding points as
 necessary to keep the angle below tv.
 
 Arguments below are:
 double tv: the total variation in orientation
 in radians.
*/
void Bezier2::PreSplit(double tv) {
    Vector v1, v2, vmid;
    int p, newloc, leftloc = leftEndpoint;
    double tl=0.0, tr = 0.0, tmid = 0.0;
    double angle1 = 0.0, angle2 = 0.0;
	
    p = right[leftEndpoint];
    while (p != -1) {
        tl = findT(leftloc);
        tr = findT(p);
        tmid = (tl+tr)/2.0;
		
        v1 = findTang(tl);
        v2 = findTang(tr);
        vmid = findTang(tmid); 
		
        //Finds middle angle to ensure that the whole rotation of the point is obtained.
        angle1 = c->anglebtw1(v1, vmid);
        angle2 = c->anglebtw1(vmid, v2);
		
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
 Find the t value of point p on bezier
 curve using the quadratic equation.
 
 Arguments below are:
 int p: the point on the bezier curve.
 */
double Bezier2::findT(int p) {
    Point aa, bb, cc;
    double tx, tx1;
	
    aa = c->pts->V[leftEndpoint] - 2 * p1 + c->pts->V[rightEndpoint];	//Coefficient of t^2
    bb = 2 * (p1 - c->pts->V[leftEndpoint]);		//Coefficient of t
    cc = c->pts->V[leftEndpoint];		//Coefficient of t^0
	
    if(aa.xy[0] == 0) {	  
        tx = (c->pts->V[p].xy[0] - cc.xy[0]) / bb.xy[0];
        tx1 = tx;	  
    } else {		
        tx = (- bb.xy[0] + sqrt(bb.xy[0] * bb.xy[0] - 4 * aa.xy[0] * (cc.xy[0]-c->pts->V[p].xy[0]))) / (2 * aa.xy[0]);
        tx1 = (- bb.xy[0] - sqrt(bb.xy[0] * bb.xy[0] - 4 * aa.xy[0] * (cc.xy[0]-c->pts->V[p].xy[0]))) / (2 * aa.xy[0]);
    }
	
    //Check if the value of t found using x values of point work for the y values
    if (fabs(aa.xy[1]*tx*tx + bb.xy[1]*tx + cc.xy[1] - c->pts->V[p].xy[1]) < 0.000001) {//FIXME: Arbitrary value	
        return tx;
    }
    if (fabs(aa.xy[1]*tx1*tx1 + bb.xy[1]*tx1 + cc.xy[1] - c->pts->V[p].xy[1]) < 0.000001) {//FIXME: Arbitrary value	
        return tx1;
    }
	
    //If those values do not work, try solving using y values.
    if(aa.xy[1] == 0) {	  
        tx = (c->pts->V[p].xy[1] - cc.xy[1]) / bb.xy[1];
        tx1 = tx;	  
    } else {		
        tx = (- bb.xy[1] + sqrt(bb.xy[1] * bb.xy[1] - 4 * aa.xy[1] * (cc.xy[1]-c->pts->V[p].xy[1]))) / (2 * aa.xy[1]);
        tx1 = (- bb.xy[1] - sqrt(bb.xy[1] * bb.xy[1] - 4 * aa.xy[1] * (cc.xy[1]-c->pts->V[p].xy[1]))) / (2 * aa.xy[1]);
    }
	
    //Check if the value of t found using y values of point work for the x values
    if (fabs(aa.xy[0]*tx*tx + bb.xy[0]*tx + cc.xy[0] - c->pts->V[p].xy[0])
            < 0.000001) {	
        return tx;
    }
    if (fabs(aa.xy[0]*tx1*tx1 + bb.xy[0]*tx1 + cc.xy[0] - c->pts->V[0].xy[0])
            < 0.000001) {
        return tx1;
    }
	
    cout << "Error: Finding t value on bezier, t:" << tx << " " << tx1 << " ";
    cout << "Point " << c->pts->V[p] << " with id " << c->pts->ids[p] << " on " << name << endl;
    return tx;
}

/*
 Return point determined by parametric t value.
 
 Arguments below are:
 double t: the t-value of the Point returned.
 */
Point Bezier2::findP(double t) {
    return (1.0 - t) * (1.0 - t) * c->pts->V[leftEndpoint]
        + 2.0 * t * (1.0 - t) * p1 + t * t * c->pts->V[rightEndpoint];
}

/*
 Return the tangent vector at t-value.
 
 Arguments below are:
 double t: the t-value of the Vector returned.
 */
Vector Bezier2::findTang(double t) {
    return 2.0 * (1.0 - t) * (p1 - c->pts->V[leftEndpoint])
        + 2.0 * t * (c->pts->V[rightEndpoint] - p1);
}

/*
 Return 1 or 2 vectors off of point p in the
 direction of the tangent line at that point.
 If p is an endpoint, only the Vector pointing
 into the Bezier is returned.
 
 Arguments below are:
 int p: the point on the bezier curve
 vector<Vector> angs: the resulting tangent Vectors
 */
void Bezier2::findAngle(int p, vector<Vector> & angs) {
    if (p == leftEndpoint) {
        angs.push_back(2.0 * (p1 - c->pts->V[leftEndpoint]));
    } else if (p == rightEndpoint) {
        // tangent vector points towards p1 in both cases...
    angs.push_back(2.0 * (p1 - c->pts->V[rightEndpoint]));
    } else {
        double t = findT(p);
        angs.push_back(2.0 * (1.0 - t) * (p1 - c->pts->V[leftEndpoint])
            + 2.0 * t * (c->pts->V[rightEndpoint] - p1));
        angs.push_back(- angs[0]);
    }
}