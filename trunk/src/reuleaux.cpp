
#include <iostream>
#include <string>

#include "allh.h"

using namespace std;

int main(int argc, char *argv[])
{
	
    string fn;
    bool isPsc = false;
	
    if (argc == 1) {
        cout << "Usage: reuleaux filename" << endl;
        return 0;
    } else {
        
        string tmp(argv[1]);
        fn = tmp;
    }

    isPsc = true;
	
    // strip .plc from the file name
    if (fn.size() > 4 && 
        fn.compare(fn.size() - 4, fn.size() - 1, ".plc") == 0) {
        fn = fn.substr(0, fn.size() - 4);
	isPsc = false;
    }
    // strip ./psc from filename
    if (fn.size() > 4 && 
        fn.compare(fn.size() - 4, fn.size() - 1, ".psc") == 0) {
        fn = fn.substr(0, fn.size() - 4);

    }
	
	
    cout << endl;
	
    PSC * complex = new PSC();
    complex->OpenLogFile(fn + ".log");
    complex->ReadCfg(fn + ".cfg");
    if (isPsc) {
        complex->ReadPsc(fn+".psc");
    }
    else {
        complex->ReadPlc(fn+".plc");
    }
    complex->PrintSVG(fn+".inp.svg",false,true,true); 
	
    // PreSplit before FixSmallAngles to pick the correct radius
    // PreSplit after FixSmallAngles to correct small angle circles
    complex->PreSplit();
    complex->PrintSVG(fn+".pre.svg",false,true,true);
    complex->FixSmallAngles();
    complex->InitBoundingBox();
    complex->PrintSVG(fn+".fsa.svg",false,true,true);
    complex->TriangulateInputVertices();
    complex->PrintSVG(fn+".ini.svg",true,true,false);
    complex->RuppertsAlgorithm();
    complex->PrintSVG(fn+".out.svg",true,true,false);

    //Point mycenter(0.5,0.75);
    //complex->PrintSVG(fn+".zoom.svg",true,true, false,mycenter, 0.25);

    complex->SetUpExteriorTriangles();
    complex->PrintSVG(fn+".int.svg",true,true,false);

    // example of how to print a zoomed image
    // this works with tangentCirclesExtreme.psc

    // good for paperexamples/tangents
    Point mycenter(1.0,0.33333);
    complex->PrintSVG(fn+".zoom.svg",true,true, false, mycenter, .2);

    // good for paperexamples/bezier
    //Point mycenter(0.0,0.0);
    //complex->PrintSVG(fn+".zoom.svg",true,true, false, mycenter, 5.77777778);

    // good for shapes/texas.psc
    //Point mycenter(4250,-100);
    //complex->PrintSVG(fn+".zoom.svg",true,true, false, mycenter, 200);

    complex->SetUpSmallAngleTriangles();
    complex->PrintSVG(fn+".rem.svg",true,true,false);

    complex->PrintNodeEle(fn+".out");

    complex->LogResults();

    cout << endl;
	
    return 0;
}
