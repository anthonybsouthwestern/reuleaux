#include "allh.h"

/*
 Arguments below are:
 
 PSC * c: the piecewise smooth complex being meshed.
 int a: the left endpoint of the linesegment.
 int b: the right endpoint.
 */
Segment::Segment(PSC * c, int a, int b):Curve(c) {

  type = CURVE_SEGMENT;
  this->c = c;
  left[b] = a;
  right[b] = -1;
  left[a] = -1;
  right[a] = b;

  leftEndpoint = a;
  rightEndpoint = b;
}

/*
 Find the midpoint between l and right[l]
 */
Point Segment::MidPoint(int l) {

  return (c->pts->V[l] + c->pts->V[right[l]])/2;
}

void Segment::PreSplit() {
  // do nothing... total variation in orientation is 0
}

void Segment::PreSplit(double tv) {
  // do nothing... total variation in orientation is 0
}


/*
This function will be useful for handling small angles. 

Arguments below are:
int l: location of the input vertex to locate points around
double r: radius of the circle to make
vector<Point> pts: vector returning 1 or 2 points
*/
void Segment::SplitAtDist(int l, double r, vector<Point> & pts) {
  
  
  int nb1 = left[l];
  int nb2 = right[l];

  if (nb1 != -1) {
    if (d2(c->pts->V[nb1],c->pts->V[l]) < r*r) {
      cout << "ERROR: SplitAtDist given an overly long length" << endl;
    }
    Vector v =c->pts->V[nb1]- c->pts->V[l];
    v.normalize();

    Point p = c->pts->V[l] + r*v;
    pts.push_back(p);
  }

  if (nb2 != -1) {
    if (d2(c->pts->V[nb2],c->pts->V[l]) < r*r) {
      cout << "ERROR: SplitAtDist given an overly long length" << endl;
    }
    Vector v =c->pts->V[nb2]- c->pts->V[l];
    v.normalize();

    Point p = c->pts->V[l] + r*v;
    pts.push_back(p);
  }

}
