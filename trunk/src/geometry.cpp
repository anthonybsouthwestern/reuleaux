#include <algorithm>
#include <cmath>

#include "allh.h"

extern "C" {
  double incircle(double * pa, double * pb, double * pc, double * pd);
  double orient2d(double * pa, double * pb, double * pc);
}

Geometry::Geometry(PSC * c1)
{
   c = c1;
   pts = c->pts;
}

double Geometry::RadiusEdgeRatio(Triangle * s)
{
  double cr,e1,e2,e3,se;
  
  cr = d(circumcenter(s), pts->V[s->v[0]]);
  
  e1 = d(pts->V[s->v[0]], pts->V[s->v[1]]);
  e2 = d(pts->V[s->v[0]], pts->V[s->v[2]]);
  e3 = d(pts->V[s->v[1]], pts->V[s->v[2]]);
   
  if (e1 < e2 && e1 < e3) {
    se = e1;
  } else if (e2 < e3) {
    se = e2;
  } else { 
    se = e3;
  }
  
  return cr / se;
}

bool Geometry::inSphereReallyStrict(int p, Triangle * s)
{
   if (d2(s->cc, pts->V[p]) < s->rsq - c->config->epsilon) {
      return true;
   }
   return false;

}

bool Geometry::inSphereReallyStrict(Point & p, Triangle * s)
{
   if (d2(s->cc, p) < s->rsq - c->config->epsilon) {
      return true;
   }
   return false;
}


bool Geometry::inSphereStrict(int p, Triangle * s)
{
   return inSphereStrict(pts->V[p], s);
}


bool Geometry::inSphereStrict(Point & p, Triangle * s)
{
  
  
  if (incircle(pts->V[s->v[0]].xy,pts->V[s->v[1]].xy,pts->V[s->v[2]].xy,p.xy)*
      orient2d(pts->V[s->v[0]].xy,pts->V[s->v[1]].xy,pts->V[s->v[2]].xy) > 0.0) {
    return true;
  }
  return false;

  return inSphereReallyStrict(p,s);
   
  // FIXME: include exact math
  //  by using exact predicates in predicates.c
  //if (c->config->exactinsphere) {
  //}
  return false;
}


bool Geometry::InDiametralBall(int i, int l, int r) {
  if (i == l || i == r) {
    return false;
  }
  return InDiametralBall(c->pts->V[i], l, r);
}

// FIXME: code this!
bool Geometry::InDiametralBall(Point p, int l, int r) {

  Point q = (c->pts->V[l] + c->pts->V[r])/2.0;

  if (d2(p,q) < d2(q,c->pts->V[l])*(1.0-c->config->epsilon) ) {
    return true;
  }  
  return false;

}


double Geometry::circumradius(Triangle * s)
{
   return d(pts->V[s->v[0]], circumcenter(s));
}

// faster version of circumcenter
// this is more robust for slivers than before...
Point Geometry::circumcenter(Triangle * t)
{
  double ax,bx,cx,ay,by,cy,d,x,y;

  ax = c->pts->V[t->v[0]].xy[0];
  ay = c->pts->V[t->v[0]].xy[1];
  bx = c->pts->V[t->v[1]].xy[0];
  by = c->pts->V[t->v[1]].xy[1];
  cx = c->pts->V[t->v[2]].xy[0];
  cy = c->pts->V[t->v[2]].xy[1];

  d = 2.0*(ay*cx + by*ax - by*cx - ay*bx - cy*ax + cy*bx);

  x = (by*ax*ax - cy*ax*ax - by*by*ay + cy*cy*ay + bx*bx*cy + ay*ay*by 
       + cx*cx*ay - cy*cy*by - cx*cx*by - bx*bx*ay + by*by*cy - ay*ay*cy)/d;
  y = (ax*ax*cx + ay*ay*cx + bx*bx*ax - bx*bx*cx + by*by*ax - by*by*cx 
       - ax*ax*bx - ay*ay*bx - cx*cx*ax + cx*cx*bx - cy*cy*ax + cy*cy*bx)/d;

  Point p(x,y);

  return p;
}

double Geometry::Volume(Triangle * s) {

  //FIXME: return the area of the triangle...


   return 0.0;
}

void Geometry::Angles(Triangle * s, vector < double >&angles) {
  Point a,b,c;

  a = pts->V[s->v[0]];
  b = pts->V[s->v[1]];
  c = pts->V[s->v[2]];

  double val;

  Vector v1, v2;

  v1 = b-a;
  v2 = c-a;
  v1.normalize();
  v2.normalize();
  val = v1*v2;
  if (val > 1.0 || val < -1.0) {
    angles.push_back(0.0);
  } else {
    angles.push_back(acos(val));
  }
  
  v1 = a-b;
  v2 = c-b;
  v1.normalize();
  v2.normalize();
  val = v1*v2;
  if (val > 1.0 || val < -1.0) {
    angles.push_back(0.0);
  } else {
    angles.push_back(acos(val));
  }

  v1 = b-c;
  v2 = a-c;
  v1.normalize();
  v2.normalize();
  val = v1*v2;
  if (val > 1.0 || val < -1.0) {
    angles.push_back(0.0);
  } else {
    angles.push_back(acos(val));
  }

  for(int i=0; i<3; i++) {
    angles[i] = angles[i]*180.0/PI;
  }
  
  // sort
  double tmp;
  if (angles[0] <= angles[1]) {
    if (angles[2] < angles[0]) {
      tmp = angles[0];
      angles[0] = angles[2];
      angles[2] = angles[1];
      angles[1] = tmp;
    } else if (angles[2] < angles[1]) {
      tmp = angles[1];
      angles[1] = angles[2];
      angles[2] = tmp;
    }
  } else if (angles[1] <= angles[2]) {
    if (angles[0] <= angles[2]) {
      tmp = angles[0];
      angles[0] = angles[1];
      angles[1] = tmp;
    } else {
      tmp = angles[0];
      angles[0] = angles[1];
      angles[1] = angles[2];
      angles[2] = tmp;
    }
  } else {
    tmp = angles[0];
    angles[0] = angles[2];
    angles[2] = tmp;
  }  
}
