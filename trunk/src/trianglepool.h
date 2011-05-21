#ifndef __TRIANGLEPOOL__
#define __TRIANGLEPOOL__

//
//  TrianglePool data structure
//
// TrianglePool is used toallocate triangles in 
// batches and reuse the memory when triangles 
// deleted.
// 
#define POOL_SIZE 5000

using namespace std;

#include "triangle.h"


class TrianglePool {
public:

   TrianglePool();

   Triangle *NewTriangle();
   Triangle *NewTriangle(Triangle * s);

   Triangle *NewTriangle(int a);
   Triangle *NewTriangle(int a, int b);
   Triangle *NewTriangle(int a, int b, int c);

   void DeleteTriangle(Triangle * s);
   int size;

private:

   Triangle * newtriangles;
   int nexttriangle;
   Triangle *reuse;


};

#endif
