#include <stdlib.h>

#include "allh.h"

Points::Points()
{
   // nothing to initialize...
}


// note: the ordering of the int and point arguments is important below
// if the int is first, it is treated as an ID for the point...
// otherwise it is treated as a hint or parent vertex...
int Points::AddPoint(int id, Point p, PointType t)
{
   id2loc[id] = V.size();
   ids[V.size()] = id;

   V.push_back(p);
   type.push_back(t);
   parent.push_back(-1);
   found.push_back(0);
   return V.size() - 1;
}

int Points::AddPoint(Point p, PointType t)
{
   V.push_back(p);
   parent.push_back(-1);
   type.push_back(t);
   found.push_back(0);
   return V.size() - 1;
}

int Points::AddPoint(Point p, int hint, PointType t)
{
   V.push_back(p);
   parent.push_back(hint);
   type.push_back(t);
   found.push_back(0);
   return V.size() - 1;
}

bool Points::IsIn(int p, Curve * s)
{
   if (l2s[p].find(s) != l2s[p].end()) {
      return true;
   }
   return false;
}

ostream & operator<<(ostream & output, PointType t)
{

   if (t == POINT_BBOX) {
      output << "POINT_BBOX";
   }
   if (t == POINT_INPUT) {
      output << "POINT_INPUT";
   }
   if (t == POINT_SEGMENT) {
      output << "POINT_SEGMENT";
   }
   if (t == POINT_2D) {
      output << "POINT_2D";
   }
   if (t == POINT_PRESPLIT) {
      output << "POINT_PRESPLIT";
   }

   return output;
}
