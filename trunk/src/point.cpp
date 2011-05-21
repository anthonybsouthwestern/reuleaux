#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>

#include "point.h"


istream & operator>>(istream & input, Point & P)
{
   input >> P.xy[0];
   input >> P.xy[1];
   return input;
}

ostream & operator<<(ostream & output, Point P)
{

   output << "(" << P.xy[0] << ", " << P.xy[1] << ")";
   return output;
}

ostream & operator<<(ostream & output, Vector P)
{

   output << "<" << P.xy[0] << ", " << P.xy[1] << ">";
   return output;
}


Point & Point::operator=(const Point & p)
{
   xy[0] = p.xy[0];
   xy[1] = p.xy[1];
   return *this;
}

int Point::operator==(Point Q)
{
   return (xy[0] == Q.xy[0] && xy[1] == Q.xy[1]);
}

int Point::operator!=(Point Q)
{
   return (xy[0] != Q.xy[0] || xy[1] != Q.xy[1]);
}


Vector Point::operator-(Point Q)
{
   Vector v;
   v.xy[0] = xy[0] - Q.xy[0];
   v.xy[1] = xy[1] - Q.xy[1];
   return v;
}

Point Point::operator+(Vector v)
{
   Point P;
   P.xy[0] = xy[0] + v.xy[0];
   P.xy[1] = xy[1] + v.xy[1];
   return P;
}

Point Point::operator-(Vector v)
{
   Point P;
   P.xy[0] = xy[0] - v.xy[0];
   P.xy[1] = xy[1] - v.xy[1];
   return P;
}

Point & Point::operator+=(Vector v)
{
   xy[0] += v.xy[0];
   xy[1] += v.xy[1];
   return *this;
}

Point & Point::operator-=(Vector v)
{
   xy[0] -= v.xy[0];
   xy[1] -= v.xy[1];
   return *this;
}


Point operator*(int c, Point Q)
{
   Point P;
   P.xy[0] = c * Q.xy[0];
   P.xy[1] = c * Q.xy[1];
   return P;
}

Point operator*(double c, Point Q)
{
   Point P;
   P.xy[0] = c * Q.xy[0];
   P.xy[1] = c * Q.xy[1];
   return P;
}

Point operator*(Point Q, int c)
{
   Point P;
   P.xy[0] = c * Q.xy[0];
   P.xy[1] = c * Q.xy[1];
   return P;
}

Point operator*(Point Q, double c)
{
   Point P;
   P.xy[0] = c * Q.xy[0];
   P.xy[1] = c * Q.xy[1];
   return P;
}

Point operator/(Point Q, int c)
{
   Point P;
   P.xy[0] = Q.xy[0] / c;
   P.xy[1] = Q.xy[1] / c;
   return P;
}

Point operator/(Point Q, double c)
{
   Point P;
   P.xy[0] = Q.xy[0] / c;
   P.xy[1] = Q.xy[1] / c;
   return P;
}

Point operator+(Point Q, Point R)
{
   Point P;
   P.xy[0] = Q.xy[0] + R.xy[0];
   P.xy[1] = Q.xy[1] + R.xy[1];
   return P;
}


double d(Point P, Point Q)
{                               // Euclidean distance
   double dx = P.xy[0] - Q.xy[0];
   double dy = P.xy[1] - Q.xy[1];
   return sqrt(dx * dx + dy * dy);
}

double d2(Point P, Point Q)
{                               // squared distance
   double dx = P.xy[0] - Q.xy[0];
   double dy = P.xy[1] - Q.xy[1];
   return (dx * dx + dy * dy);
}



Vector Vector::operator-()
{
   Vector v;
   v.xy[0] = -xy[0];
   v.xy[1] = -xy[1];
   return v;
}

//------------------------------------------------------------------
//  Scalar Ops
//------------------------------------------------------------------

// Scalar multiplication
Vector operator*(int c, Vector w)
{
   Vector v;
   v.xy[0] = c * w.xy[0];
   v.xy[1] = c * w.xy[1];
   return v;
}

Vector operator*(double c, Vector w)
{
   Vector v;
   v.xy[0] = c * w.xy[0];
   v.xy[1] = c * w.xy[1];
   return v;
}

Vector operator*(Vector w, int c)
{
   Vector v;
   v.xy[0] = c * w.xy[0];
   v.xy[1] = c * w.xy[1];
   return v;
}

Vector operator*(Vector w, double c)
{
   Vector v;
   v.xy[0] = c * w.xy[0];
   v.xy[1] = c * w.xy[1];
   return v;
}

Vector operator/(Vector w, int c)
{
   Vector v;
   v.xy[0] = w.xy[0] / c;
   v.xy[1] = w.xy[1] / c;
   return v;
}

Vector operator/(Vector w, double c)
{
   Vector v;
   v.xy[0] = w.xy[0] / c;
   v.xy[1] = w.xy[1] / c;
   return v;
}


Vector Vector::operator+(Vector w)
{
   Vector v;
   v.xy[0] = xy[0] + w.xy[0];
   v.xy[1] = xy[1] + w.xy[1];
   return v;
}

Vector Vector::operator-(Vector w)
{
   Vector v;
   v.xy[0] = xy[0] - w.xy[0];
   v.xy[1] = xy[1] - w.xy[1];
   return v;
}

double Vector::operator*(Vector w)
{
   return (xy[0] *w.xy[0] + xy[1] *w.xy[1]);
}

double Vector::operator*(Point w)
{
   return (xy[0] *w.xy[0] + xy[1] *w.xy[1]);
}

Vector Vector::operator^(Vector w)
{
   Vector v;
   v.xy[0] = xy[1] *w.xy[2] - xy[2] *w.xy[1];
   v.xy[1] = xy[2] *w.xy[0] - xy[0] *w.xy[2];
   return v;
}

Vector & Vector::operator*=(double c)
{
   xy[0] *= c;
   xy[1] *= c;
   return *this;
}

Vector & Vector::operator/=(double c)
{
   xy[0] /= c;
   xy[1] /= c;
   return *this;
}

Vector & Vector::operator+=(Vector w)
{
   xy[0] += w.xy[0];
   xy[1] += w.xy[1];
   return *this;
}

Vector & Vector::operator-=(Vector w)
{
   xy[0] -= w.xy[0];
   xy[1] -= w.xy[1];
   return *this;
}


void Vector::normalize()
{
   double ln = sqrt(xy[0] *xy[0] + xy[1] *xy[1]);
   if (ln == 0) {
      return;                   // do nothing to zero vector...
   }

   xy[0] /= ln;
   xy[1] /= ln;
}
