#include "common.h"
#include "triangle.h"

Triangle::Triangle(Triangle * s)
{
   v[0] = s->v[0];
   v[1] = s->v[1];
   v[2] = s->v[2];
   n[0] = s->n[0];
   n[1] = s->n[1];
   n[2] = s->n[2];
   rsq = s->rsq;
   cc = s->cc;
   visited = s->visited;
}

Triangle::Triangle()
{
   v[0] = -1;
   v[1] = -1;
   v[2] = -1;
   n[0] = NULL;
   n[1] = NULL;
   n[2] = NULL;
   rsq = -1.0;
   visited = 0;
}

Triangle::Triangle(int a, int b, int c)
{
   v[0] = a;
   v[1] = b;
   v[2] = c;
   n[0] = NULL;
   n[1] = NULL;
   n[2] = NULL;
   rsq = -1.0;
   visited = 0;
}


// set pointer to adjacent triangles
void Triangle::link(Triangle * s1, Triangle * s2, Triangle * s3)
{
   n[0] = s1;
   n[1] = s2;
   n[2] = s3;
}

void Triangle::SetVertex(int i, int myv)
{
   v[i] = myv;
}

void Triangle::SetNeighbor(int i, Triangle * myn)
{
   n[i] = myn;
}

// does triangle contain vertex i?
bool Triangle::contains(int i)
{
   if (v[0] == i || v[1] == i || v[2] == i) {
      return true;
   }
   return false;
}

// is triangle i a neighbor to the triangle?
bool Triangle::ncontains(Triangle * i)
{
   if (n[0] == i || n[1] == i || n[2] == i) {
      return true;
   }
   return false;
}

// do my triangle and triangle s share two vertices so 
// they should be linked together?
int Triangle::isNeighbor(Triangle * s)
{
   int a, b, c;
   a = 0;
   b = 0;
   c = 0;

   if (s->contains(v[0])) {
      a = 1;
   }
   if (s->contains(v[1])) {
      b = 1;
   }
   if (s->contains(v[2])) {
      c = 1;
   }

   if (a + b + c < 3) {
      return -1;
   } else {
      if (a == 0) {
         return 0;
      }
      if (b == 0) {
         return 1;
      }
      if (c == 0) {
         return 2;
      }
   }
   return -2;
}


bool operator<(Triangle & x, Triangle & y)
{

  for (int i = 0; i < 3; i++) {
    if (x.v[i] < y.v[i]) {
      return true;
    }
    if (x.v[i] > y.v[i]) {
      return false;
    }
  }

  return false;
}

// ordering on faces... this is important for
// declaring sets of faces
bool operator<(const Face & x, const Face & y)
{

   if (x.s < y.s) {
      return true;
   } else if (x.s == y.s && x.v < y.v) {
      return true;
   }
   return false;
}

// usefule output routine
ostream & operator<<(ostream & output, Triangle * s)
{
   if (s == NULL) {
      output << "[]" << endl;
      return output;
   }

   output << "[" << s->v[0] << ", " << s->v[1] << ", " << s->v[2] << "]";

   return output;
}
