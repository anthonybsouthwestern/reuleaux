#include "allh.h"
#include <iostream>
#include <fstream>

#include "iniparser.h"

Configuration::Configuration(PSC * C1)
{

   C = C1;

   ini = NULL;

   //epsilon = .000000000001;
   //infinity = 1000000000000.0;

   epsilon = .00000001;
   infinity = 100000000.0;

   boundingboxratio = 6.0;
   maxradiusedgeratio = 1.5;

   maxoutputradius = -1.0;
   hsq = -1.0;
   presplitvariation = -1.0;

   maxiterations = -1;

   exactinsphere = true;

   printtplc = true;

   verbose = false;

   removeboundingbox = false;
	
   printNBezierWSeg = false;

}

void Configuration::ReadConfiguration(string fname)
{

   ini = iniparser_load(fname.c_str());

   if (ini == NULL) {
      if (C->logfileisopen) {
         (*(C->LOGFILE)) << "Config file " << fname << " not found.   Defaults used."
                         << endl;
      }

   } else {
      (*(C->LOGFILE)) << "Config file: " << fname
                      <<" loaded successfully."  << endl;

   // double options
   epsilon = iniparser_getdouble(ini,"reuleaux:epsilon",epsilon);
   infinity = iniparser_getdouble(ini,"reuleaux:infinity",infinity);
   maxradiusedgeratio = iniparser_getdouble(ini,"reuleaux:maxradiusedgeratio", maxradiusedgeratio);
   maxoutputradius = iniparser_getdouble(ini,"reuleaux:maxoutputradius", maxoutputradius);
   boundingboxratio = iniparser_getdouble(ini,"reuleaux:boundingboxratio",boundingboxratio);
   presplitvariation = iniparser_getdouble(ini, "reuleaux:presplitvariation", presplitvariation);
	
   if (maxoutputradius >0.0)
     hsq = maxoutputradius*maxoutputradius;

   // integer options
   maxiterations = iniparser_getint(ini,"reuleaux:maxiterations",maxiterations);

   // boolean options
   exactinsphere = iniparser_getboolean(ini,"reuleaux:exactinsphere",exactinsphere);


   printtplc = iniparser_getboolean(ini,"reuleaux:printtplc",printtplc);

   verbose = iniparser_getboolean(ini,"reuleaux:verbose",verbose);

   removeboundingbox = iniparser_getboolean(ini,"reuleaux:removeboundingbox",removeboundingbox);
	   
   printNBezierWSeg = iniparser_getboolean(ini, "reuleaux:printNBezierWSeg", printNBezierWSeg);
   

   }

   LogConfiguration();
}

void Configuration::LogConfiguration()
{
   if (C->logfileisopen) {
     (*(C->LOGFILE)) << setfill('_') << setw(25) << left << "Epsilon" 
		     << setfill('_') << setw(20) << right << epsilon << endl;
     (*(C->LOGFILE)) << setfill('_') << setw(25) << left << "Infinity" 
		     << setfill('_') << setw(20) << right << infinity << endl;
     (*(C->LOGFILE)) << setfill('_') << setw(25) << left << "BoundingBoxRatio" 
		     << setfill('_') << setw(20) << right << boundingboxratio << endl;
     (*(C->LOGFILE)) << setfill('_') << setw(25) << left << "MaxRadiusEdgeRatio" 
		     << setfill('_') << setw(20) << right << maxradiusedgeratio << endl;
     (*(C->LOGFILE)) << setfill('_') << setw(25) << left << "ExactInsphere" 
		     << setfill('_') << setw(20) << right << exactinsphere << endl;
     (*(C->LOGFILE)) << setfill('_') << setw(25) << left << "RemoveBoundingBox" 
		     << setfill('_') << setw(20) << right << removeboundingbox << endl;
   }
}
