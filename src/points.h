#ifndef __POINTS__
#define __POINTS__

#include <map>
#include <set>
#include <vector>

#include "common.h"


enum PointType { POINT_BBOX,POINT_INPUT,POINT_SEGMENT,POINT_2D, POINT_PRESPLIT};

class Points {
public:

   Points();

   vector < Point > V;          // points
   vector < PointType > type;   // POINT_BBOX, POINT_3D, etc

   map < int, int >ids;         // map from point's index to its id in the datafile...
   map < int, int >id2loc;      // id to loc map from input...

   map < int, set < Curve * > >l2s;      // list of containing segments...
   vector < int >parent;        // this is not really used at this point...

   vector<int> found;

   int AddPoint(Point p, PointType t);
   int AddPoint(Point p, int hint, PointType t);
   int AddPoint(int id, Point p, PointType t);

   bool IsIn(int p, Curve * s);
};

ostream & operator<<(ostream &, PointType);

#endif                          // !defined(__POINTS__)
