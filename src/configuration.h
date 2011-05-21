#ifndef __CONFIGURATION_H__
#define __CONFIGURATION_H__

//
// Configuration class
//
// Read user parameters from the .cfg file and 
// set default parameter values.
//


#include <iostream>
#include <string>
#include <vector>
#include <set>

#include "common.h"

#include "iniparser.h"

using namespace std;


class Configuration {
public:

   double epsilon;
   double infinity;

   double maxradiusedgeratio;
   double boundingboxratio;
   double maxoutputradius;
   double hsq;
   double presplitvariation;
   
   int maxiterations;

   bool exactinsphere;

   bool printtplc;

   bool removeboundingbox;
   bool printNBezierWSeg;

   bool verbose;

   PSC *C;

   Configuration(PSC * C1);
   void ReadConfiguration(string fname);
   void PrintConfiguration();
   void LogConfiguration();

   dictionary * ini;

private:
};

#endif              // !defined(__CONFIGURATION_H__)
