#ifndef __Common_H__
#define __Common_H__


#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <vector>
#include <time.h>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <deque>
#include <list>
#include <algorithm>

#include "point.h"

#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH

class Points;
class Geometry;
class PSC;
class Curve;
class Segment;
class Circle;
class Subsegment;
class Triangulation;
class Triangle;
class TrianglePool;
class Configuration;
class Stats;


struct Face {
   Triangle *s;
   int v;
};

#endif                          //__Common_H__
