#include "common.h"
#include "trianglepool.h"

using namespace std;


TrianglePool::TrianglePool()
{

   reuse = NULL;
   newtriangles = new Triangle[POOL_SIZE];
   nexttriangle = 0;
   size = POOL_SIZE;
}

Triangle *TrianglePool::NewTriangle()
{
  Triangle *newtriangle = NULL;
  if (reuse != NULL) {
    // try to reuse first...
    newtriangle = reuse;
    reuse = newtriangle->n[0];
    
  } else if (nexttriangle < POOL_SIZE) {
    // get the next in the pool...
    newtriangle = &newtriangles[nexttriangle];
    nexttriangle++;
  } else {
    size += POOL_SIZE;
    
    //cout << "Allocating more space " << size  << endl;
    // allocate more space...
    newtriangles = new Triangle[POOL_SIZE];
    nexttriangle = 0;
    newtriangle = &newtriangles[nexttriangle];
    nexttriangle++;
  }
  
  newtriangle->v[0] = -1;
  newtriangle->v[1] = -1;
  newtriangle->v[2] = -1;
  newtriangle->n[0] = NULL;
  newtriangle->n[1] = NULL;
  newtriangle->n[2] = NULL;
  newtriangle->rsq = -1.0;
  
  return newtriangle;
}

Triangle *TrianglePool::NewTriangle(int a, int b, int c)
{
   Triangle *newtriangle = NewTriangle();
   newtriangle->v[0] = a;
   newtriangle->v[1] = b;
   newtriangle->v[2] = c;

   return newtriangle;

}

Triangle *TrianglePool::NewTriangle(Triangle * s)
{
   Triangle *newtriangle = NewTriangle();
   newtriangle->v[0] = s->v[0];
   newtriangle->v[1] = s->v[1];
   newtriangle->v[2] = s->v[2];
   newtriangle->n[0] = s->n[0];
   newtriangle->n[1] = s->n[1];
   newtriangle->n[2] = s->n[2];
   newtriangle->rsq = s->rsq;
   newtriangle->cc = s->cc;

   return newtriangle;
}

void TrianglePool::DeleteTriangle(Triangle * s)
{

    s->n[0] = reuse;
    reuse = s;

   s->rsq = -1.0;
   s->v[0] = -1;
   s->v[1] = -1;
   s->v[2] = -1;
   // n[0] points to the list of triangles to reuse...
   s->n[1] = NULL;
   s->n[2] = NULL;
}
