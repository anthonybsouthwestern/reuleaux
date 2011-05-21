#ifndef __GEOMETRY__
#define __GEOMETRY__

#include <vector>
#include <set>
#include <map>

#define EPSILON .000000000001

//
//  Geometry class
//
//  This is a collection of geometric utilities.
//  Includes functions for computing circumcenters,
//  circumradii, and angles of triangles, among other
//  things.
//

class Geometry {

public:

   Points * pts;
   PSC *c;

   Geometry(PSC * c1);
   
   bool inSphereStrict(int p, Triangle * s);
   bool inSphereStrict(Point & p, Triangle * s);
   bool inSphereReallyStrict(int p, Triangle * s);
   bool inSphereReallyStrict(Point & p, Triangle * s);

   bool InDiametralBall(int i, int l, int r);
   bool InDiametralBall(Point p, int l, int r);

   Point circumcenter(Triangle * s);

   double circumradius(Triangle * s);

   double RadiusEdgeRatio(Triangle * s);

// never used!
   double Orient(Triangle * s);

   double Volume(Triangle * s);
   void Angles(Triangle * s, vector < double >&angles);

private:

};


#endif
