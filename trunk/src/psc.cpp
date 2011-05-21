#include "allh.h"


PSC::PSC() {
    pts = new Points();
    geom = new Geometry(this);
    config = new Configuration(this);
    stats = new Stats();
    triPool = new TrianglePool();
}


/*
 Opens a file for log information of the mesh for
 the input PSC.
 
 Arguments below are:
 string filename: the name of the logfile
 */
void PSC::OpenLogFile(string filename)
{
    LOGFILE = new ofstream(filename.c_str());
    logfileisopen = true;
}

/*
 Reads the contents of the configuration file
 associated with the PSC file.
 
 Arguments below are:
 sting filename: the name of the configuration file
 */
void PSC::ReadCfg(string filename) {
    PrintBegin("ReadCfg(): Reading Configuration File.");
    config->ReadConfiguration(filename);
    PrintEnd("ReadCfg(): Reading Configuration File.");
}

/*
 This routine currently reads the PSC file type
 containing linesegments, circles, circular arcs,
 and Bezier curves (quadratic, cubic, and nth degree).
 
 Arguments below are:
 string fname: the name of the PSC file being read
 */
void PSC::ReadPsc(string fname)	
{	
    PrintBegin("ReadPsc(): Reading Data From " + fname);
	
    map <int, int> e2p1;
    map <int, int> e2p2;
	
    tri = new Triangulation(this);
	
    // read an input file...
    ifstream in((fname).c_str());
    if (!in) {
        cerr << "ERROR: can't open file:" << fname << endl;
        exit(0);
    }
	
    string junk;
	
    int nVert, dim;
	
    in >> junk;
	
    //Ensure it is a PSC file
    if (junk.compare("PSC") != 0) {
        if (logfileisopen) {
            (*LOGFILE) << "Specified file is not in correct format." 
         << endl;
        }
    } else {
        if (logfileisopen) {
            (*LOGFILE) << "PSC file " << fname << " opened successfully."
            << endl;
        }
    }
	
    in >> junk;
    //Get the dimensions of the PSC (currently code only works with 2 dimensions)
    if (junk.compare("DIMENSION") != 0) {
        cout << "Error reading dimension section of data file." << endl;
        exit(0);
    }
	
    in >> dim;
	
	
    if (dim != 2) {
        cout << "ERROR: mesh dimension is " << dim << "." <<
        endl;
        exit(0);
    }
	
    in >> junk >> nVert;
	
    if (junk.compare("POINTS") != 0) {
        cout << "Error reading points section of data file." << endl;
        exit(0);
    }
	
    if (logfileisopen) {
        (*LOGFILE) << "Reading " << nVert << " Points" << endl;
    }
	
    // read the nodes...
    int id;
	
    for (int ii=0; ii<2; ii++) {
        tri->upperbound[ii] = -1.0*config->infinity;
        tri->lowerbound[ii] = config->infinity;
    }
    double *p = new double[2]; //2 dimensional point
    p[0] = 0.0;
    p[1] = 0.0;
	
    vector < int >myids;
    vector < Point > mypoints;
    map < int, int >id2point;
	
    //Add points to the triangulation, while updating the bounds of the triangulation
    for (int i = 0; i < nVert; i++) {
        in >> id;
        for (int j = 0; j < dim; j++) {
            double tmp;
            in >> tmp;
            p[j] = tmp;
            if (tmp > tri->upperbound[j]) {
                tri->upperbound[j] = tmp;
            }
            if (tmp < tri->lowerbound[j]) {
                tri->lowerbound[j] = tmp;
            }
        }
        Point p1(p[0], p[1]);
        mypoints.push_back(p1);
        myids.push_back(id);
        id2point[id] = i;
    
        //Add loc to points data structure
        int loc = pts->AddPoint(id, p1, POINT_INPUT);
        inputVertices.push_back(loc);
		
    }
	
    delete[] p;
	
    // This temporarily sets the bbox for the .in.svg output file
    Point center((tri->upperbound[0]+tri->lowerbound[0])/2.0,
                (tri->upperbound[1]+tri->lowerbound[1])/2.0);
	
    tri->center = center;
	
    double bboxsize = max(tri->upperbound[0]-tri->lowerbound[0],
                         tri->upperbound[1]-tri->lowerbound[1])/2.0;
	
    tri->bbox = bboxsize * config->boundingboxratio;
	
    // add curves...
    int nCurves;
	
    in >> junk >> nCurves;
	
    if (junk.compare("CURVES") != 0) {
    cout << "Error reading lines section of data file." << endl;
        exit(0);
    }
	
	
    if (logfileisopen) {
        (*LOGFILE) << "Reading " << nCurves << " Curves" << endl;
    }
	
    for (int i = 0; i < nCurves; i++) {
        int edgeid, npoints, ida, idb, degree;
        double loca, locb, radius, locc, locd;
		
        in >> edgeid >> junk;
		
        if(junk.compare("linesegment")==0){
			
            in >> npoints;
            in >> idb;
            //Creates a seperate linesegment between each set of points.
            for (int j = 0;j < npoints; j++){
                ida = idb;
                in >> idb;
                loca = pts->id2loc[ida];
                locb = pts->id2loc[idb];
				
                Segment * newSeg = new Segment(this,loca,locb);
				
                ostringstream oss1;
                oss1 << "Line-" << edgeid*10 + j;
                newSeg->name = oss1.str();
                newSeg->id = edgeid;
                //Add curve to points data structure
                pts->l2s[loca].insert(newSeg);
                pts->l2s[locb].insert(newSeg);
				
                curves.push_back(newSeg);
            }
        } else if(junk.compare("circle")==0) {
			
            in >> loca >> locb >> radius;
			
            Point cp(loca, locb);
			
			
            Circle * newCir = new Circle(this, cp, radius);
			
            in >> npoints;
            //Add intersection points to circle
            //Note: can add other points as well, algorithm is most efficient with only these points
            for(int j = 0; j < npoints; j++){ 
                in >> id;
                ida = pts->id2loc[id];
				
                if (j == 0) {
                    newCir->InsertVertex(ida);
                } else {
                    newCir->InsertVertex(ida, idb);
                }
                pts->l2s[ida].insert(newCir);
                idb = ida;
            }
			
			
            ostringstream oss1;
            oss1 << "Circle-" << edgeid;
            newCir->name = oss1.str();
            newCir->id = edgeid;
			
            curves.push_back(newCir);
			
        } else if(junk.compare("arc")==0) {
            in >> ida >> idb >> loca >> locb;
			
            Point cp(loca, locb);
			
            Arc * newArc = new Arc(this, cp, pts->id2loc[ida], pts->id2loc[idb]);
			
            pts->l2s[pts->id2loc[ida]].insert(newArc);
            pts->l2s[pts->id2loc[idb]].insert(newArc);
			
            in >> npoints;
            //Add intersection points to circular arc
            for(int j = 0; j < npoints; j++){ 
                in >> id;
                loca = pts->id2loc[id];
				
                newArc->InsertVertex(loca);
				
            }
			
			
            ostringstream oss1;
            oss1 << "Arc-" << edgeid;
            newArc->name = oss1.str();
            newArc->id = edgeid;
			
            curves.push_back(newArc);
        } else if(junk.compare("bezier")==0) {	//quadratic Bezier curve
            in >> ida >> idb >> loca >> locb;
			
            Bezier2 * newBez = new Bezier2(this, pts->id2loc[ida], pts->id2loc[idb], loca, locb);
			
            pts->l2s[pts->id2loc[ida]].insert(newBez);
            pts->l2s[pts->id2loc[idb]].insert(newBez);
			
            ostringstream oss1;
            oss1 << "Bezier-" << edgeid;
            newBez->name = oss1.str();
            newBez->id = edgeid;
			
            in >> npoints;
            //Add intersection points to quadratic bezier
            for (int j = 0; j < npoints; j++) {
                in >> id;
                newBez->InsertVertex(pts->id2loc[id]);
            }
            curves.push_back(newBez); 
			
        } else if (junk.compare("bezier3")==0) {	//cubic Bezier curve
			
            in >> ida >> idb >> loca >> locb >> locc >> locd;
			
            Point p1(loca, locb);
            Point p2(locc, locd);
			
            Bezier3 * newCubBez = new Bezier3(this, pts->id2loc[ida], pts->id2loc[idb], p1, p2);
			
            pts->l2s[pts->id2loc[ida]].insert(newCubBez);
            pts->l2s[pts->id2loc[idb]].insert(newCubBez);
			
            ostringstream oss1;
            oss1 << "CubicBezier-" << edgeid;
            newCubBez->name = oss1.str();
            newCubBez->id = edgeid;
			
            in >> npoints;
            //Add intersection points to cubic bezier
            for (int j = 0; j < npoints; j++) {
                in >> id;
                newCubBez->InsertVertex(pts->id2loc[id]);
            }
            curves.push_back(newCubBez); 
        }else if (junk.compare("nbezier")==0) {	//nth degree Bezier curve
			
            in >> ida >> idb >> degree >> npoints;
			
            NBezier * newNthBez = new NBezier(this, pts->id2loc[ida], pts->id2loc[idb], degree);
			
            for (int j = 1; j < degree; j++) {
                in >> loca >> locb;
                Point p(loca, locb);
                newNthBez->addPt(p, j);
            }
			
            pts->l2s[pts->id2loc[ida]].insert(newNthBez);
            pts->l2s[pts->id2loc[idb]].insert(newNthBez);
			
            ostringstream oss1;
            oss1 << "BezierDegree" << degree << "-" << edgeid;
            newNthBez->name = oss1.str();
            newNthBez->id = edgeid;
			
            //Add intersection points to n-degree bezier
            ida = newNthBez->leftEndpoint;
            for (int j = 0; j < npoints; j++) {
                in >> idb;
                ((Curve *) newNthBez)->InsertVertex(pts->id2loc[idb], pts->id2loc[ida]);
                ida = idb;
            }
            curves.push_back(newNthBez); 
        }
    }
	
    PrintEnd("ReadPsc(): Reading Data");
}

/*
 This routine currently reads the PLC file type
 containing only segments and no smooth curves.
 
 Arguments below are:
 string fname: the name of the PLC file being read
 */
void PSC::ReadPlc(string fname)	
{
    PrintBegin("ReadPlc(): Reading Data From " + fname);
	
    map <int, int> e2p1;
    map <int, int> e2p2;
	
    tri = new Triangulation(this);
	
    // read an input file...
    ifstream in((fname).c_str());
    if (!in) {
        cerr << "ERROR: can't open file:" << fname << endl;
        exit(0);
    }
	
    string junk;
	
    int nVert, dim;
	
    in >> junk;
	
    if (junk.compare("PLC") != 0) {
        if (logfileisopen) {
            (*LOGFILE) << "Specified file is not in correct format." 
         << endl;
        }
    } else {
        if (logfileisopen) {
            (*LOGFILE) << "PLC file " << fname << " opened successfully."
            << endl;
        }
    }
	
    in >> junk;
    if (junk.compare("DIMENSION") != 0) {
        cout << "Error reading dimension section of data file." << endl;
        exit(0);
    }
	
    in >> dim;
	
	
    if (dim != 2) {
        cout << "ERROR: mesh dimension is " << dim << "." <<
        endl;
        exit(0);
    }
	
    in >> junk >> nVert;
	
    if (junk.compare("POINTS") != 0) {
        cout << "Error reading points section of data file." << endl;
        exit(0);
    }
	
    if (logfileisopen) {
        (*LOGFILE) << "Reading " << nVert << " Points" << endl;
    }
	
    // read the nodes...
    int id;
	
    for (int ii=0; ii<2; ii++) {
        tri->upperbound[ii] = -1.0*config->infinity;
        tri->lowerbound[ii] = config->infinity;
    }
    double *p = new double[2];
    p[0] = 0.0;
    p[1] = 0.0;
	
    vector < int >myids;
    vector < Point > mypoints;
    map < int, int >id2point;
	
    for (int i = 0; i < nVert; i++) {
        in >> id;
        for (int j = 0; j < dim; j++) {
            double tmp;
            in >> tmp;
            p[j] = tmp;
            if (tmp > tri->upperbound[j]) {
                tri->upperbound[j] = tmp;
            }
            if (tmp < tri->lowerbound[j]) {
                tri->lowerbound[j] = tmp;
            }
        }
        Point p1(p[0], p[1]);
		
        mypoints.push_back(p1);
        myids.push_back(id);
        id2point[id] = i;
		
        int loc = pts->AddPoint(id, p1, POINT_INPUT);
        inputVertices.push_back(loc);
		
    }
	
    delete[] p;
	
    Point center((tri->upperbound[0]+tri->lowerbound[0])/2.0,
                (tri->upperbound[1]+tri->lowerbound[1])/2.0);
	
    tri->center = center;
	
    double bboxsize = max(tri->upperbound[0]-tri->lowerbound[0],
                         tri->upperbound[1]-tri->lowerbound[1])/2.0;
	
    tri->bbox = bboxsize * config->boundingboxratio;
	
    tri->Initialize(); // set up the initial square
	
    tri->bbox = tri->bbox + config->epsilon;  //create a little wiggle room... 
	
    if (logfileisopen) {
        (*LOGFILE) << "Bounding box center: " << tri->center << endl;
        (*LOGFILE) << "Bounding box size:   " << tri->bbox << endl;
    }
	
    // add segments...
    int nSegments;
	
    in >> junk >> nSegments;

    if (junk.compare("SEGMENTS") != 0) {
        cout << "Error reading lines section of data file." << endl;
        exit(0);
    }
	
	
    if (logfileisopen) {
        (*LOGFILE) << "Reading " << nSegments << " Segments" << endl;
    }
	
    for (int i = 0; i < nSegments; i++) {
        int edgeid, ida,idb,loca,locb;
		
        in >> edgeid >> ida >> idb;
		
        loca = pts->id2loc[ida];
        locb = pts->id2loc[idb];
		
        Segment * newSeg = new Segment(this,loca,locb);
		
        ostringstream oss1;
        oss1 << "Line-" << edgeid;
        newSeg->name = oss1.str();
        newSeg->id = edgeid;
		
        pts->l2s[loca].insert(newSeg);
        pts->l2s[locb].insert(newSeg);
		
        curves.push_back(newSeg);
    }
	
    PrintEnd("ReadPlc(): Reading Data");
}

/*
 Reduce the variation in orientation of all input
 curves. This means the max allowed angle change
 of a tangent line as it traces a subsegment of a
 curve cannot be more than 45 degrees (90 degrees
 for circles and circular arcs).
 */
void PSC::PreSplit() {
    PrintBegin("PreSplit(): Reduce all curves to 90 degrees");
	
    // loop through curves and call presplit
    foreach (Curve * c, curves) {
        if (config->presplitvariation > 0.0) {
            c->PreSplit(config->presplitvariation);
        } else {
            c->PreSplit();
        }
    }
}

/*
 Creates a Delaunay triangulation from the initial input.
 */
void PSC::TriangulateInputVertices() {
    //Adds all input vertices to the triangulation construct.
    PrintBegin("TriangulateInputVertices(): Computing Delaunay Triangulation");
    foreach (int loc, inputVertices) {
        tri->InsertPoint(loc);
		
    }
    //Adds all of the additional presplit vertices to the triangulation.
    foreach (int loc, presplitVertices) {
        tri->InsertPoint(loc);
		
    }
	
    PrintEnd("TriangulateInputVertices(): Computing Delaunay Triangulation");
}

/*
 Performs the core of Ruppert's Algorithm on the 
 Delaunay triangulation.
 */
void PSC::RuppertsAlgorithm() {
    PrintBegin("RuppertsAlgorithm(): Refining Triangulation");	
	
	
	
    // queue encroached segments
    foreach (Curve * s, curves) {
        Curve::iterator it;
        it=s->begin();
        do {
            Subsegment ss;
            ss.leftEndpoint = *it;
            ss.inputCurve = s;
			
            if (IsEncroached(ss) && !IsProtected(ss)) {
                segmentQueue.insert(ss);      
            }
            it++;
        } while (it != s->end());
    }
	
    // queue bad triangles
    for (Triangulation::titerator it=tri->tbegin();
        it != tri->tend(); it++) {
        Triangle * t = *it;		
        if (geom->RadiusEdgeRatio(t) > config->maxradiusedgeratio && 
            !IsProtected(t) ){  
            triangleQueue.insert(t);
        }
    }
	
    if (logfileisopen) {
        (*LOGFILE) << "Initial Queue Size: " << segmentQueue.size()
      << " segments and " << triangleQueue.size() 
      << " triangles." << endl;
    }
	
	
    bool insertVertex = false;
    int vertexLocation = -1;
    Triangle * vertexContainedIn = NULL;
    while (segmentQueue.size() > 0 || triangleQueue.size() > 0) {
		
        insertVertex = false;
		
        if (segmentQueue.size() > 0) {
			
            set<Subsegment>::iterator it = segmentQueue.begin();
            Subsegment ss = *it;
            segmentQueue.erase(it);
			
            int left = ss.leftEndpoint;	//Used after checkSteinerPoint to see if subsegment ss was affected
            int right = ss.inputCurve->right[left];
			
            Point p = ss.inputCurve->MidPoint(ss.leftEndpoint);
			
            //Test to see if midpoint is within 2 * the
            //distance to a potential small angle circle radius.
            // If so, inserts that small angle circle into the mesh
            if(checkSteinerPoint(p)) {	//No small angle circle added
                vertexLocation = pts->AddPoint(p,left,POINT_SEGMENT);
				
                ss.inputCurve->InsertVertex(vertexLocation, left);
				
                vertexContainedIn = tri->Locate(vertexLocation,left);
				
                insertVertex = true;
            } else	//Check to see if segment queue was changed by the addition of the circle
                if (IsEncroached(ss) && !IsProtected(ss) && ss.inputCurve->right[left] == right) {
					
                    segmentQueue.insert(ss);
                }
			
			
        } else {	//Work on bad triangles
			
            set<Triangle *>::iterator it = triangleQueue.begin();
            vertexContainedIn = *it;
            triangleQueue.erase(it);
			
            Point p = geom->circumcenter(vertexContainedIn);
			
            if(checkSteinerPoint(p, vertexContainedIn)) {
                // should I yield?
				
                vector<Subsegment> segmentsEncroached;
                if (DoesEncroach(p,vertexContainedIn,segmentsEncroached)) {
                    foreach(Subsegment ss, segmentsEncroached) {
                        if (!IsProtected(ss)) {
                            segmentQueue.insert(ss);
                        }
                    }
                } else {
                    vertexLocation = pts->AddPoint(p,vertexContainedIn->v[0],POINT_2D);
                    insertVertex = true;
                }
				
            }//Check to see if triangle queue was affected by the addition of circle
            else if (geom->RadiusEdgeRatio(vertexContainedIn) > config->maxradiusedgeratio && !IsProtected(vertexContainedIn) && 
                  tri->Locate(p) == vertexContainedIn) {
				
                triangleQueue.insert(vertexContainedIn);
            }
        }
        //A steiner point was added to improve on the triangulation
        if(insertVertex) {
            vector < Triangle * > oldCavity;
            vector < Triangle * > newCavity;
			
            tri->InsertPoint(vertexLocation, vertexContainedIn, oldCavity, newCavity);
			
            foreach (Triangle * t1, oldCavity) {
                // remove from the queue
                set<Triangle *>::iterator it1 = triangleQueue.find(t1);
                if(it1 != triangleQueue.end()) {
                    triangleQueue.erase(t1);
                }
				
                // return the memory
                triPool->DeleteTriangle(t1);
            }
			
            if (pts->type[vertexLocation] == POINT_SEGMENT) {
                // check new subsegments for encroachment...
                Subsegment lefts;
                Subsegment rights;
				
                Curve * inputSeg = *(pts->l2s[vertexLocation].begin());
				
                lefts.leftEndpoint = vertexLocation;
                lefts.inputCurve = inputSeg;
				
                if (IsEncroached(lefts) && !IsProtected(lefts)) {
                    segmentQueue.insert(lefts);
                }
				
                rights.leftEndpoint = inputSeg->left[vertexLocation];
                rights.inputCurve = inputSeg;
				
                if (rights.leftEndpoint != 1 && IsEncroached(rights) && !IsProtected(rights) ) {
                    segmentQueue.insert(rights);
                }
            }
			
			
            // check encroachment...
            vector<Subsegment> segmentsEncroached;
            vector<int> nbs;
            tri->Neighbors(vertexLocation, nbs);
            if (DoesEncroach(pts->V[vertexLocation],nbs,segmentsEncroached)) {
                foreach(Subsegment ss, segmentsEncroached) {
                    if (!IsProtected(ss)) {
                        segmentQueue.insert(ss);
                    }
                }
            }
			
            foreach (Triangle * t1, newCavity) {
                // check radius edge ratio and queue
                if (geom->RadiusEdgeRatio(t1) > config->maxradiusedgeratio && !IsProtected(t1) ) {
                    triangleQueue.insert(t1);
                }
            }
        }
    }
    PrintEnd("RuppertsAlgorithm(): Refining Triangulation");
	
}


/*
 Check to see if the addition of a Steiner point
 to the triangulation will be too close to a small
 angle vertex. If so, the circle is then added to
 the triangulation, and Ruppert's algorithm
 continues as usual.
 
 Arguments below are:
 Point p: the added steiner point
 */
bool PSC::checkSteinerPoint(Point & p) {
    return checkSteinerPoint(p, tri->Locate(p));
}

/*
 Check to see if the addition of a Steiner point
 to the triangulation will be too close to a small
 angle vertex. If so, the circle is then added to
 the triangulation, and Ruppert's algorithm
 continues as usual.
 
 Arguments below are:
 Point p: the added steiner point
 Triangle t: the triangle containing the added steiner point
 */
bool PSC::checkSteinerPoint(Point &p, Triangle * t) {
    vector<int> nbs;
    tri->Cavity(p, t, nbs);

    for (int i = 0; i < (int)nbs.size(); i++) {
        if (vertex2Radius.find(nbs[i]) != vertex2Radius.end()) {
            if (d(p, pts->V[nbs[i]]) < 2 * vertex2Radius[nbs[i]]) {
				
                addSmallAngleCircle(nbs[i], vertex2Radius[nbs[i]]);
				
                vertex2Radius.erase(nbs[i]);
				
                return false;
            }
        }
    }
	
    return true;
}

/*
 Check to see if the addition of a Steiner point
 to the triangulation will be too close to a small
 angle vertex. If so, the circle is then added to
 the triangulation, and Ruppert's algorithm
 continues as usual.
 
 Arguments below are:
 Point p: the added steiner point
 vector<Triangle *> tris: a vector of triangles that contain the added steiner point
 */
bool PSC::checkSteinerPoint(Point & p, vector<Triangle *> tris) {
    for (int j = 0;j < tris.size();j++) {
        if(!checkSteinerPoint(p, tris[j]))
            return false;
    }
    return true;
}

/*
 Determines if a subsegment's diametral circle
 contains a point in the triangulation.
 
 Arguments below are:
 Subsegment s: the subsegment
 */
bool PSC::IsEncroached(Subsegment s) {
    int l = s.leftEndpoint;
    int r = s.inputCurve->right[l];
	
    if (r == -1) {
        return false;
    }
	
    vector<int> nbs;
    tri->Neighbors(l,nbs);
    tri->Neighbors(r,nbs);// adds the right neighbors to set
    foreach(int nb, nbs) {
		
        // does it encroach?
        if (nb != r && nb != l) {
            if (geom->InDiametralBall(nb,l,r)) {
                return true;
            }      
        }
    }
    nbs.clear();
	
    return false;
	
}

/*
 Does p encroach a subcurve with an endpoint in the
 neighbors of the triangle t?
 
 Arguments below are:
 Point p: the point.
 Triangle *t: the triangle near point p.
 set<Subsegment> s: is returned with the subsegments that were encroached.
 */
bool PSC::DoesEncroach(Point & p, Triangle * t, set<Subsegment> & s) {
    set<int> nbs;
    tri->Cavity(p,t, nbs);
	
    return DoesEncroach(p, nbs, s);
	
}

/*
 Does p encroach a subcurve with an endpoint in the
 neighbors of the triangle t?
 
 Arguments below are:
 Point p: the point.
 Triangle *t: the triangle near point p.
 vector<Subsegment> s: is returned with the subsegments that were encroached.
 */
bool PSC::DoesEncroach(Point & p, Triangle * t, vector<Subsegment> & s) {
    vector<int> nbs;
    tri->Cavity(p,t, nbs);
	
    return DoesEncroach(p, nbs, s);
	
}

/*
 Does p encroach a subcurve with an endpoint in nbs?
 
 Arguments below are:
 Point p: the point.
 set<int> nbs: neighboring points of p.
 set<Subsegment> s: is returned with the subsegments that were encroached.
 */
bool PSC::DoesEncroach(Point & p, set<int> nbs, set<Subsegment> & encroached) {
    foreach (int iNb, nbs) {
        foreach (Curve * seg, pts->l2s[iNb]) {
            // check segments for encroachment...
            int left, right;
			
            left = seg->left[iNb];
            right = seg->right[iNb];
			
            // check the two subsegments...
            if (left != -1) {
                if (geom->InDiametralBall(p,left,iNb)) {
                    Subsegment tmp;
                    tmp.leftEndpoint = left;
                    tmp.inputCurve = seg;
                    encroached.insert(tmp);
                }      
            }
			
            if (right != -1) {
                if (geom->InDiametralBall(p,iNb,right)) {
                    Subsegment tmp;
                    tmp.leftEndpoint = iNb;
                    tmp.inputCurve = seg;
                    encroached.insert(tmp);
                }      
            }
        } 
    }
    if (encroached.size() > 0) {
        return true;
    }
	
    return false;
	
}

/*
 Does p encroach a subcurve with an endpoint in nbs?
 
 Arguments below are:
 Point p: the point.
 vector<int> nbs: neighboring points of p.
 vector<Subsegment> s: is returned with the subsegments that were encroached.
 */
bool PSC::DoesEncroach(Point & p, vector<int> nbs, vector<Subsegment> & encroached) {
    foreach (int iNb, nbs) {
        pts->found[iNb] = 1;
        foreach (Curve * seg, pts->l2s[iNb]) {
            // check segments for encroachment...
            int left, right;
			
            left = seg->left[iNb];
            right = seg->right[iNb];
			
            // check the two subsegments...
            if (left != -1 && pts->found[left] == 0) {
                if (geom->InDiametralBall(p,left,iNb)) {
                    Subsegment tmp;
                    tmp.leftEndpoint = left;
                    tmp.inputCurve = seg;
                    encroached.push_back(tmp);
                }      
            }
			
            if (right != -1 && pts->found[right] == 0) {
                if (geom->InDiametralBall(p,iNb,right)) {
                    Subsegment tmp;
                    tmp.leftEndpoint = iNb;
                    tmp.inputCurve = seg;
                    encroached.push_back(tmp);
                }      
            }
        } 
    }
	
    foreach (int iNb, nbs) {
        pts->found[iNb] = 0;
    }
	
    if (encroached.size() > 0) {
        return true;
    }
	
    return false;
}

/*
 Print the mesh in Triangle format.
 
 Arguments below are:
 string filename: the name of the file the mesh output is printed to.
 */
void PSC::PrintMesh(string filename)
{
    ofstream psout;
    psout.open(filename.c_str());
	
    double bbox = tri->bbox;
	
    double xmin = tri->center.xy[0]-bbox/config->boundingboxratio,
    ymin = tri->center.xy[1]-bbox/config->boundingboxratio,
    xmax = tri->center.xy[0]+bbox/config->boundingboxratio,
    ymax = tri->center.xy[1]+bbox/config->boundingboxratio;
	
    // Scale and translate so that
    // lower LH corner is at (1,1) inches = (72,72) points
    // upper RH corner scaled to fit on page with 1 inch margin
	
    double res = 550;
	
    double sx = res / (xmax - xmin);
    double sy = res / (ymax - ymin);
	
    double scale = fmin(sx,sy);
	
	
    double center = res/2.0;
    double width = res/2.0;
    double height = width;
	
    psout << "%!" << endl
   << "%%Title: tri2ps" << endl
   << "%%BoundingBox: " << (int)(center-width)-10 << " " << (int)(center-height)-10 << " "
   << (int)( center+width)+10 << " "
   << (int)( center+height)+10 << endl
   << "/setc { aload pop setrgbcolor } bind def" << endl
   << "/col1 { [1 0 0] setc } def" << endl
   << "/col2 { [0 0 1] setc } def" << endl
   << "/col0 { 0 setgray } def" << endl
	<< "%%Page: 1 1" << endl
   << "%%BeginPageSetup" << endl
   << "%%EndPageSetup" << endl
   << "1.0 setlinewidth" << endl;
	
    //Add the triangulation to the SVG file
    for (Triangulation::titerator it=tri->tbegin();
        it != tri->tend(); it++) {
		
        Triangle * t = *it;
		
		
        Point p0,p1,p2;
        p0 = pts->V[t->v[0]];
        p1 = pts->V[t->v[1]];
        p2 = pts->V[t->v[2]];
		
        if (p0.xy[0] >= xmin && p0.xy[0]<=xmax &&
          p0.xy[1] >= ymin && p0.xy[1]<=ymax &&
          p1.xy[0] >= xmin && p1.xy[0]<=xmax &&
          p1.xy[1] >= ymin && p1.xy[1]<=ymax &&
          p2.xy[0] >= xmin && p2.xy[0]<=xmax &&
          p2.xy[1] >= ymin && p2.xy[1]<=ymax) {
			
            psout << "newpath "
            << (int) (scale*(p0.xy[0] - xmin))+10 << " "
            << (int) (scale*(p0.xy[1] - ymin))+10 << " moveto "
            << (int) (scale*(p1.xy[0] - xmin))+10 << " "
            << (int) (scale*(p1.xy[1] - ymin))+10 << " lineto ";
            psout << "col0 stroke" << endl;
			
            psout << "newpath "
            << (int) (scale*(p0.xy[0] - xmin))+10 << " "
            << (int) (scale*(p0.xy[1] - ymin))+10 << " moveto "
            << (int) (scale*(p2.xy[0] - xmin))+10 << " "
            << (int) (scale*(p2.xy[1] - ymin))+10 << " lineto ";
            psout << "col0 stroke" << endl;
			
            psout << "newpath "
            << (int) (scale*(p2.xy[0] - xmin))+10 << " "
            << (int) (scale*(p2.xy[1] - ymin))+10 << " moveto "
            << (int) (scale*(p1.xy[0] - xmin))+10 << " "
            << (int) (scale*(p1.xy[1] - ymin))+10 << " lineto ";
            psout << "col0 stroke" << endl;
        }
    }
	
    psout << "showpage" << endl
   << "%%PageTrailer" << endl
   << "%%EndPage" << endl;
	
    psout.close();
}

/*
 Print the mesh in .node format.
 
 Arguments below are:
 string filename: the name of the file the mesh output is printed to.
 */
void PSC::PrintNodeEle(string filename) {
    string fnnode = filename + ".node";
    string fnele = filename + ".ele";
	
    ofstream fnode(fnnode.c_str());
    ofstream fele(fnele.c_str());
	
	
    fnode << pts->V.size() << " 2 0 0"  << endl;
    for (unsigned int i=0; i<pts->V.size(); i++) {
        fnode << i << " " << pts->V[i].xy[0] << " " << pts->V[i].xy[1] << endl;
    }  
	
	
    fele << tri->nt << " 3 0" << endl;
    int count = 0;
    for (Triangulation::titerator it=tri->tbegin(); it != tri->tend(); it++) {
        Triangle * t = *it;
        fele << count << " " << t->v[0] << " " << t->v[1] << " " << t->v[2] << endl;
        count++;
    }
	
}

/*
 Determines if any input curve connects the point a
 to the point b.
 
 Arguments below are:
 int a = the 'a' integer linked to a point construct
 int b = the 'b' integer linked to a point construct
 */
bool PSC::IsArc(int a,int b) {
    foreach (Curve * c, pts->l2s[a]) {
        if (c->left[a] == b || c->right[a] == b) {
            return true;
        }
    }
    return false;
}

/*
 Outputs the triangulation and curves to a 
 SVG file if the center and bounding box to 
 be printed are not given.
 
 Arguments below are:
 string filename: the name of the SVG file being written to.
 bool printTriangulation: if true, prints the triangulation created by Ruppert's
 bool noTriangulationWithArcs: if true, segments of the triangulation are removed if they overlap with a curve segment.
 bool printVertices: if true, prints small filled circles at each vertex in the mesh.
 */
void PSC::PrintSVG(string filename, bool printTriangulation, bool noTriangulationWithArcs, bool printVertices) {
    PrintSVG(filename, printTriangulation, noTriangulationWithArcs, printVertices, tri->center, tri->bbox/(config->boundingboxratio*.9));
}

/*
 Outputs the triangulation and curves to a SVG file.
 
 Arguments below are:
 string filename = the name of the SVG file being written to.
 bool printTriangulation = if true, prints the triangulation created by Ruppert's
 bool noTriangulationWithArcs = if true, segments of the triangulation are removed if they overlap with a curve segment.
 bool printVertices = if true, prints small filled circles at each vertex in the mesh.
 Point center = the center of the output printed
 double size = the distance out from the center that is printed to the SVG file
 */
void PSC::PrintSVG(string filename, bool printTriangulation, bool noTriangulationWithArcs, bool printVertices, Point center, double size){	
    ofstream svgout;
    svgout.open(filename.c_str());
	
    double xmin = center.xy[0]-size,
    ymin = center.xy[1]-size,
    xmax = center.xy[0]+size,
    ymax = center.xy[1]+size;
	
    double linewidth = (xmax-xmin)/500.0;
	
    svgout << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl 
    << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl
    << " \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl
    << " <svg width=\"16cm\" height=\"16cm\" viewBox=\"" 
    << xmin-5.0*linewidth << " " << ymin-5.0*linewidth << " " 
    << xmax-xmin+10.0*linewidth << " " << ymax-ymin+10.0*linewidth  << "\"" << endl
    << "version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">" << endl
    << "<desc>Triangular mesh produced by Reuleaux</desc>" << endl
    << "<style type=\"text/css\"><![CDATA[.Thin { fill:none; stroke:black; stroke-width:" 
    << linewidth << " }]]></style>" << endl
    << "<style type=\"text/css\"><![CDATA[.Thick { fill:none; stroke:black; stroke-width:" 
    << linewidth*3 << " }]]></style>" << endl
    << "<style type=\"text/css\"><![CDATA[.SmallAng { fill:none; stroke:red; stroke-width:" 
    << linewidth*2 << " }]]></style>" << endl;
   
    //Print the triangulation
    if (printTriangulation) {
        for (Triangulation::titerator it=tri->tbegin();
           it != tri->tend(); it++) {
			
            Triangle * t = *it;
			
            if (!IsExterior(t) && !IsInterior(t)) {
				
                Point p0,p1,p2;
                p0 = pts->V[t->v[0]];
                p1 = pts->V[t->v[1]];
                p2 = pts->V[t->v[2]];
				
                // gets the triangulation to print all the way to the bounding box
                if (p0.xy[0] >= xmin || p0.xy[0]<=xmax ||     //If triangle is inside 
                p0.xy[1] >= ymin || p0.xy[1]<=ymax || // desired bounding box.
                p1.xy[0] >= xmin || p1.xy[0]<=xmax ||
                p1.xy[1] >= ymin || p1.xy[1]<=ymax ||
                p2.xy[0] >= xmin || p2.xy[0]<=xmax ||
                p2.xy[1] >= ymin || p2.xy[1]<=ymax) {
					
					
                    if (noTriangulationWithArcs) {
                        // are v0 and v1 endpoints of an arc?
                        if (!IsArc(t->v[0],t->v[1])) {
                            svgout << "<path class=\"Thin\" d=\"M "
                            << p0.xy[0] << "," << 2.0*center.xy[1]-p0.xy[1] << " L " 
                            << p1.xy[0] << "," << 2.0*center.xy[1]-p1.xy[1] << " \"/>" << endl;			      
                        }
                        if (!IsArc(t->v[0],t->v[2])) {
                            svgout << "<path class=\"Thin\" d=\"M "
                            << p0.xy[0] << "," << 2.0*center.xy[1]-p0.xy[1] << " L " 
                            << p2.xy[0] << "," << 2.0*center.xy[1]-p2.xy[1] << " \"/>" << endl;			      
                        }
                        if (!IsArc(t->v[2],t->v[1])) {
                            svgout << "<path class=\"Thin\" d=\"M "
                            << p2.xy[0] << "," << 2.0*center.xy[1]-p2.xy[1] << " L " 
                            << p1.xy[0] << "," << 2.0*center.xy[1]-p1.xy[1] << " \"/>" << endl;			      
                        }
						
                    } else {
                        svgout << "<path class=\"Thin\" d=\"M "
                        << p0.xy[0] << "," << 2.0*center.xy[1]-p0.xy[1] << " L " 
                        << p1.xy[0] << "," << 2.0*center.xy[1]-p1.xy[1] << " L " 
                        << p2.xy[0] << "," << 2.0*center.xy[1]-p2.xy[1] << " z\"/>" << endl;
                    }
                }
            }
        }
    }
	
    //Print the curves
    for (vector<Curve * >::iterator it = curves.begin();			//Output of original input segments.
        it != curves.end(); it++) {
        //Dont need to check boundingbox area, because it must encase input shapes.
        Curve * c = *it;
        int p0, l, r;
		
        if(c->type == CURVE_SEGMENT) {
            l = c->leftEndpoint;
            r = c->rightEndpoint;
			
            // Flip y-axis to match what appears correct
            svgout << "<path class=\"Thick\" d=\"M "
            << pts->V[l].xy[0] << "," << 2.0*center.xy[1]-pts->V[l].xy[1] << " L "
            << pts->V[r].xy[0] << "," << 2.0*center.xy[1]-pts->V[r].xy[1] << "\"/>" << endl;
    
        }
		
        if(c->type == CURVE_CIRCLE) {
			
            Circle * cir = (Circle * ) c;
			
            p0 = cir->leftEndpoint;
            l = p0;
			
            do{
                r = cir->right[l];
                if (cir->name.substr(0, 8) == "SAcircle") { //Small angle protection circles
                    svgout << "<path class=\"SmallAng\" d=\"M ";
                } else { //Input circles
                    svgout << "<path class=\"Thick\" d=\"M ";
                }
				
                // Flip y-axis to match what appears correct
                svgout << pts->V[l].xy[0] << "," << 2.0*center.xy[1]-pts->V[l].xy[1] << " A "	
                << cir->radius << "," <<cir->radius << " 0 0,1 " //x-axis-rotation large-arc-flag ',' sweep-flag
                << pts->V[r].xy[0] << "," << 2.0*center.xy[1]-pts->V[r].xy[1] << "\"/>" << endl; 
                l = r;
				
            } while (l != p0);
        }
        if(c->type == CURVE_ARC){
            Arc * a = (Arc *) c;
			
            l = a->leftEndpoint;
            do{
                r = a->right[l];
                // Flip y-axis to match what appears correct
                svgout << "<path class=\"Thick\" d=\"M "
                << pts->V[l].xy[0] << "," << 2.0*center.xy[1]-pts->V[l].xy[1] << " A "	
                << a->radius << "," <<a->radius << " 0 0,1 " //x-axis-rotation large-arc-flag ',' sweep-flag
                << pts->V[r].xy[0] << "," << 2.0*center.xy[1]-pts->V[r].xy[1] << "\"/>" << endl; 
                l = r;
				
            } while (l != a->rightEndpoint);
        }
        if(c->type == CURVE_BEZIER) {
            // Flip y-axis to match what appears correct
            svgout << "<path class=\"Thick\" d=\"M "
            << pts->V[c->leftEndpoint].xy[0] << "," <<  2.0*center.xy[1]-pts->V[c->leftEndpoint].xy[1] << " Q "
            << ((Bezier2 * )c)->p1.xy[0] << "," << 2.0*center.xy[1]-((Bezier2 * )c)->p1.xy[1] << " "
            << pts->V[c->rightEndpoint].xy[0] << "," <<  2.0*center.xy[1]-pts->V[c->rightEndpoint].xy[1] << "\"/>" << endl;
        }
        if(c->type == CURVE_CUBBEZ) {
            // Flip y-axis to match what appears correct
            svgout << "<path class=\"Thick\" d=\"M "
            << pts->V[c->leftEndpoint].xy[0] << "," <<  2.0*center.xy[1]-pts->V[c->leftEndpoint].xy[1] << " C "
            << ((Bezier3 * )c)->controlPts[0].xy[0] << "," << 2.0*center.xy[1]-((Bezier3 * )c)->controlPts[0].xy[1] << " "
            << ((Bezier3 * )c)->controlPts[1].xy[0] << "," << 2.0*center.xy[1]-((Bezier3 * )c)->controlPts[1].xy[1] << " "
            << pts->V[c->rightEndpoint].xy[0] << "," <<  2.0*center.xy[1]-pts->V[c->rightEndpoint].xy[1] << "\"/>" << endl;
        }
        if (c->type == CURVE_NBEZIER) {
            if (config->printNBezierWSeg) {	//Use a segment approximation to display nth degree Bezier
                vector<Point> outpts;
                ((NBezier *) c)->getOutput(100, outpts);
				
                for(int i = 1;i < (int)outpts.size();i++) {
                    // Flip y-axis to match what appears correct
                    svgout << "<path class=\"Thick\" d=\"M "
                    << outpts[i - 1].xy[0] << "," << 2.0*center.xy[1]-outpts[i - 1].xy[1] << " L "
                    << outpts[i].xy[0] << "," << 2.0*center.xy[1]-outpts[i].xy[1] << "\"/>" << endl;
                }
            } else {	//Use quadratic Beziers to approximate the nth degree Bezier
                vector<Point> quadBezPts;
                ((NBezier *) c)->get2Beziers(quadBezPts);
				
                for (int i = 0;i < (int)quadBezPts.size() - 1;i+=2) {
                    // Flip y-axis to match what appears correct
                    svgout << "<path class=\"Thick\" d=\"M "
                    << quadBezPts[i].xy[0] << "," <<  2.0*center.xy[1]-quadBezPts[i].xy[1] << " Q "
                    << quadBezPts[i+1].xy[0] << "," << 2.0*center.xy[1]-quadBezPts[i+1].xy[1] << " "
                    << quadBezPts[i+2].xy[0] << "," <<  2.0*center.xy[1]-quadBezPts[i+2].xy[1] << "\"/>" << endl;
                }
            }
        }
    }
	
    if (printVertices) {
        //print the vertices here..
        for (int i=0; i<(int)pts->V.size(); i++) {
            if (pts->type[i] == POINT_INPUT) {
                svgout << "<circle cx=\"" << pts->V[i].xy[0] <<"\" cy=\"" 
                << 2.0*center.xy[1]-pts->V[i].xy[1] << "\" r=\"" << linewidth*2 
                << "\" stroke=\"black\" stroke-width=\"" 
                << linewidth << "\" fill=\"black\"/>" << endl;
            }
            if (pts->type[i] == POINT_PRESPLIT) {
				
                bool isSmallAng = false;
                foreach (Curve * c, pts->l2s[i]) {
                    if (c->name.substr(0,8) == "SAcircle") {
                        isSmallAng = true;
                    }
                }
				
                if (isSmallAng) {
                    svgout << "<circle cx=\"" << pts->V[i].xy[0] <<"\" cy=\"" 
                    << 2.0*center.xy[1]-pts->V[i].xy[1] << "\" r=\"" << linewidth*2 
                    << "\" stroke=\"red\" stroke-width=\"" 
                    << linewidth << "\" fill=\"red\"/>" << endl;
                } else {
                    svgout << "<circle cx=\"" << pts->V[i].xy[0] <<"\" cy=\"" 
                    << 2.0*center.xy[1]-pts->V[i].xy[1] << "\" r=\"" << linewidth*2 
                    << "\" stroke=\"blue\" stroke-width=\"" 
                    << linewidth << "\" fill=\"blue\"/>" << endl;
                }
            }
        }
    }
	
    svgout<< "</svg>" <<endl;
	
}

/*
 Prints the input given by the PSC file in SVG format.
 
 Arguments below are:
 string filename: the SVG file name
 */
void PSC::PrintInput(string filename)
{
    ofstream psout;
    psout.open(filename.c_str());
	
    double bbox = tri->bbox;
	
    double xmin = tri->center.xy[0]-bbox/config->boundingboxratio,
    ymin = tri->center.xy[1]-bbox/config->boundingboxratio,
    xmax = tri->center.xy[0]+bbox/config->boundingboxratio,
    ymax = tri->center.xy[1]+bbox/config->boundingboxratio;
	
    double res = 550;
	
    double sx = res / (xmax - xmin);
    double sy = res / (ymax - ymin);
	
    double scale = fmin(sx,sy);
	
	
    double center = res/2.0;
    double width = res/2.0;
    double height = width;
	
    psout << "%!" << endl
   << "%%Title: tri2ps" << endl
   << "%%BoundingBox: " << (int)(center-width)-10 << " " << (int)(center-height)-10 << " "
   << (int)( center+width)+10 << " "
   << (int)( center+height)+10 << endl
   << "/setc { aload pop setrgbcolor } bind def" << endl
   << "/col1 { [1 0 0] setc } def" << endl
   << "/col2 { [0 0 1] setc } def" << endl
   << "/col0 { 0 setgray } def" << endl
    << "%%Page: 1 1" << endl
   << "%%BeginPageSetup" << endl
   << "%%EndPageSetup" << endl
   << "3.0 setlinewidth" << endl;
	
	
	
    foreach(Curve * s, curves) {
        Curve::iterator it;
        for (it=s->begin(); it!=s->end(); it++) {
            int left  = *it;
			
			
            int right = s->right[left];
            if (right != -1) {
                Point p0 = pts->V[left];
                Point p1 = pts->V[right];
				
                psout << "newpath "
                << (int) (scale*(p0.xy[0] - xmin))+10 << " "
                << (int) (scale*(p0.xy[1] - ymin))+10 << " moveto "
                << (int) (scale*(p1.xy[0] - xmin))+10 << " "
                << (int) (scale*(p1.xy[1] - ymin))+10 << " lineto ";
                psout << "col0 stroke" << endl;	
            }
        }
    }
	
    psout << "showpage" << endl
   << "%%PageTrailer" << endl
   << "%%EndPage" << endl;
	
    psout.close();
}


/*
 Prints the quantitative results from the meshing to the logfile.
 */
void PSC::LogResults() {
    PrintBegin("LogResults(): Printing Mesh Summary.");
	
    int np1=0,nt1=0;
    int minhist[12];
    int maxhist[12];
	
    for (int i=0; i<12; i++) {
        minhist[i] = 0;
        maxhist[i] = 0;
    }
	
    double minangle,maxangle;
    minangle = 90.0;
    maxangle = 0.0;
	
    for (Triangulation::titerator it=tri->tbegin();
        it != tri->tend(); it++) {  
        nt1++;
		
        Triangle * t = *it;
		
        vector<double> angles;
        geom->Angles(t,angles);
		
        // make a histogram of angles
        int minind = (int)(angles[0]/5.0);
        int maxind = (int)((angles[2]-60.0)/10.0);
		
        minhist[minind]++;
        maxhist[maxind]++;
		
        if (angles[0] < minangle) {
            minangle = angles[0];
        }
        if (angles[2] > maxangle) {
            maxangle = angles[2];
        }
		
    }
	
    int pointhist[5];
    for (int i=0; i<5; i++) {
        pointhist[i] = 0;
    }
	
	
    for (Triangulation::viterator it=tri->vbegin();
        it != tri->vend(); it++) {  
        np1++;
		
        if (pts->type[*it] == POINT_BBOX)
            pointhist[0]++;
        if (pts->type[*it] == POINT_INPUT)
            pointhist[1]++;
        if (pts->type[*it] == POINT_SEGMENT)
            pointhist[2]++;
        if (pts->type[*it] == POINT_2D)
            pointhist[3]++;
        if (pts->type[*it] == POINT_PRESPLIT)
            pointhist[4]++;
		
    }
	
    if (nt1 != tri->nt) {
        *LOGFILE << "Triangle Count Mismatch: " << tri->nt << " reported,  " 
        << nt1 << " actual" << endl;
    } else {
        (*LOGFILE) <<setw(25) << left << setfill('_') << "Number of Triangles" 
        <<setw(20) << right << setfill('_') << tri->nt << endl;    
    }
    if (np1 != tri->np) {
        *LOGFILE << "Vertex Count Mismatch: " << tri->np << " reported,  " 
        << np1 << " actual" << endl;
    } else {
        (*LOGFILE) <<setw(25) << left << setfill('_') << "Number of Vertices" 
        <<setw(20) << right << setfill('_') << tri->np << endl;    
    }
	
    (*LOGFILE) << endl;
    (*LOGFILE) <<setw(25) << left << setfill('_') << "Minimum Angle" 
    <<setw(20) << right << setfill('_') << minangle << endl;
    (*LOGFILE) <<setw(25) << left << setfill('_') << "Maximum Angle" 
    <<setw(20) << right << setfill('_') << maxangle << endl;
	
    (*LOGFILE) << endl << "Number of triangles by minimum angle." << endl;
        for (int i = 0; i < 12; i++) {
        (*LOGFILE) << setfill(' ');
        (*LOGFILE) << " " << setw(3) << i * 5 << " - " << setw(3) << i * 5 +
        5;
        (*LOGFILE) << setw(11) << setfill('_') << minhist[i];
        if (i%2 == 0) {
            (*LOGFILE) << "   ";
        } else {
            (*LOGFILE) << endl;
        }
    }
	
    (*LOGFILE) << endl << "Number of triangles by maximum angle." << endl;
    for (int i = 0; i < 12; i++) {
        (*LOGFILE) << setfill(' ');
        (*LOGFILE) << " " << setw(3) << i * 10+60 << " - " << setw(3) << i * 10+60 +
        10;
        (*LOGFILE) << setw(11) << setfill('_') << maxhist[i];
        if (i%2 == 0) {
            (*LOGFILE) << "   ";
        } else {
            (*LOGFILE) << endl;
        }
    }
    (*LOGFILE) << setfill(' ');
	
    (*LOGFILE) << endl << "Number of vertices by type." << endl;
    (*LOGFILE) << setfill('_') << setw(25) << left << "POINT_BBOX" 
    << setfill('_') << setw(20) << right << pointhist[0] << endl;
    (*LOGFILE) << setfill('_') << setw(25) << left << left << "POINT_INPUT" 
    << setfill('_') << setw(20) << right << pointhist[1] << endl;
    (*LOGFILE) << setfill('_') << setw(25) << left << left << "POINT_SEGMENT" 
    << setfill('_') << setw(20) << right << pointhist[2] << endl;
    (*LOGFILE) << setfill('_') << setw(25) << left << left << "POINT_2D" 
    << setfill('_') << setw(20) << right << pointhist[3] << endl;
    (*LOGFILE) << setfill(' ');
    (*LOGFILE) << setfill('_') << setw(25) << left << left << "POINT_PRESPLIT" 
    << setfill('_') << setw(20) << right << pointhist[4] << endl;
    (*LOGFILE) << setfill(' ');
    PrintEnd("LogResults(): Printing Mesh Summary.");
}

/*
 Print the beginning information from the input
 to the logfile.
 
 Arguments below are:
 string s: the string printed at the beginning of the logfile.
 */
void PSC::PrintBegin(string s)
{
	
    cout << s << endl;
    if (logfileisopen) {
        (*LOGFILE) << endl <<
        "-------------------------------------------------" << endl;
        (*LOGFILE) << "BEGIN: " << s << endl;
        (*LOGFILE) << endl;
    }
    stats->logtime = clock();
}

/*
 Prints the final timing results and end of the logfile.
 
 Arguments below are:
 string s: the string printed at the beginning of the logfile.
 */
void PSC::PrintEnd(string s)
{
    if (logfileisopen) {
        (*LOGFILE) << endl;
        if (tri != NULL) {
            (*LOGFILE) <<setw(25) << left << setfill('_') << "Vertices in the Mesh" 
            <<setw(20) << right << setfill('_') << tri->np << endl;    
            (*LOGFILE) <<setw(25) << left << setfill('_') << "Triangles in the Mesh" 
            <<setw(20) << right << setfill('_') << tri->nt << endl;    
        }
        (*LOGFILE) <<setw(25) << left << setfill('_') << "Seconds Elapsed" 
        <<setw(20) << right << setfill('_') << ((double) (clock() - stats->logtime)) << endl;    
        (*LOGFILE) << setfill(' ');
        (*LOGFILE) << endl << "END: " << s << endl;
        (*LOGFILE) << "-------------------------------------------------"
        << endl;
    }
}

/*
 Finds angles smaller than 60 degrees and then 
 sets the radius for a circle around the point
 where a small angle is.
 */
void PSC::FixSmallAngles()
{
    PrintBegin("FixSmallAngles(): Add circles on points with angles less than 60 degrees");
	
    foreach(int vertex, inputVertices) {
		
        vector<CurveAngle> ca;
        double radius = INFINITY;
		
        if (! pts->l2s[vertex].empty()) {	//If there are curves on the input vertex
			
            foreach(Curve * s, pts->l2s[vertex]) {//Get the tangent vector for each curve off of the vertex.
                Vector tmp;
                if (s->type == CURVE_SEGMENT) {
                    int loc = s->left[vertex];
                    if (loc == -1)
                        loc = s->right[vertex];
					
                    tmp = pts->V[loc] - pts->V[vertex];
                    if (tmp.len() < radius) {
                        radius = tmp.len();
                    }
					
                    tmp.normalize();
					
                    CurveAngle tmpca;
                    tmpca.dir = tmp;
                    tmpca.c = s;
                    tmpca.curvatureInfo = 0.0;
                    if (ca.size() == 0) {
                        tmpca.angle = 0.0;
                    } else {
                        tmpca.angle = anglebtw(ca[0].dir,tmp);
                    }
                    ca.push_back(tmpca);
                }
                if (s->type == CURVE_CIRCLE) {
                    tmp = ((Circle *)s)->center - pts->V[vertex];
					
                    if (tmp.len() * 2 < radius) {
                        radius = tmp.len() * 2;
                    }
					
                    Vector v1 = Vector(tmp.xy[1], -tmp.xy[0]); // 90 degree rotation from vector to center pt of circle
                    Vector v2 = -v1;
					
                    v1.normalize();
                    v2.normalize();
					
                    CurveAngle tmpca;
                    tmpca.dir = v1;
                    tmpca.c = s;
                    tmpca.curvatureInfo = -1.0/((Circle *)s)->radius;
                    if (ca.size() == 0) {
                        tmpca.angle = 0.0;
                    } else {
                        tmpca.angle = anglebtw(ca[0].dir,v1);
                    }
                    ca.push_back(tmpca);
                    CurveAngle tmpca1;
                    tmpca1.dir = v2;
                    tmpca1.c = s;
                    tmpca1.curvatureInfo = 1.0/((Circle *)s)->radius;
                    if (ca.size() == 0) {
                        tmpca1.angle = 0.0;
                    } else {
                        tmpca1.angle = anglebtw(ca[0].dir,v2);
                    }
                    ca.push_back(tmpca1);
                }
                if (s->type == CURVE_ARC) {
                    tmp = ((Arc *)s)->center - pts->V[vertex];
					
                    if (tmp.len() * 2 < radius) {
                        radius = tmp.len() * 2;
                    }
					
                    Vector v1 = Vector(tmp.xy[1], -tmp.xy[0]); // 90 degree rotation from vector to center pt of circle
                    Vector v2 = -v1;
					
                    v1.normalize();
                    v2.normalize();
                    if (vertex != s->leftEndpoint) {
						
                        CurveAngle tmpca;
                        tmpca.dir = v1;
                        tmpca.c = s;
                        tmpca.curvatureInfo = -1.0/((Arc *)s)->radius;
                        if (ca.size() == 0) {
                            tmpca.angle = 0.0;
                        } else {
                            tmpca.angle = anglebtw(ca[0].dir,v1);
                        }
                        ca.push_back(tmpca);
                    }
                    if (vertex != s->rightEndpoint) {
                        CurveAngle tmpca1;
                        tmpca1.dir = v2;
                        tmpca1.c = s;
                        tmpca1.curvatureInfo = 1.0/((Arc *)s)->radius;
                        if (ca.size() == 0) {
                            tmpca1.angle = 0.0;
                        } else {
                            tmpca1.angle = anglebtw(ca[0].dir,v2);
                        }
                        ca.push_back(tmpca1);
                    }
                }
                //Uses a linesegment approximation of the Bezier curves to locate the nearest point
                if (s->type == CURVE_BEZIER || s->type == CURVE_CUBBEZ) {
                    vector<Vector> dirs;
                        Vector vec[4];
                    double bezrad;
                    if (s->type == CURVE_BEZIER) {
                        ((Bezier2 *) s)->findAngle(vertex, dirs);
						
                        double tvert = ((Bezier2 *) s)->findT(vertex);
                        vec[0] = pts->V[s->rightEndpoint] - pts->V[s->leftEndpoint];
                        if (tvert < 0.5) {
                            vec[1] = pts->V[s->leftEndpoint] - ((Bezier2 *) s)->p1;
                        }else {
                            vec[1] = pts->V[s->rightEndpoint] - ((Bezier2 *) s)->p1;
                        }
                        vec[2] = ((Bezier2 *) s)->p1 - pts->V[vertex];
						
                        if (radius > vec[0].len() * vec[2].len() / vec[1].len()) {
                            radius = vec[0].len() * vec[2].len() / vec[1].len();
                        }
                        Vector v1 = ((Bezier2 *) s)->findTang(tvert);
                        bezrad = v1.len2();
						
                    }else {
                        ((Bezier3 *) s)->findAngle(vertex, dirs);
                        vec[0] = pts->V[s->rightEndpoint] - pts->V[s->leftEndpoint];
                        vec[1] = ((Bezier3 *) s)->controlPts[0] - pts->V[vertex];
                        vec[2] = ((Bezier3 *) s)->controlPts[1] - pts->V[vertex];
                        if (vec[2].len() < vec[1].len()) {
                            vec[1] = vec[2];
                            vec[2] = ((Bezier3 *) s)->controlPts[1] - pts->V[s->rightEndpoint];
                        } else {
                            vec[2] = ((Bezier3 *) s)->controlPts[0] - pts->V[s->leftEndpoint];
                        }
                        if (vec[0].len() * vec[1].len() / (2 * vec[2].len()) < radius) {
                            radius = vec[0].len() * vec[1].len() / (2 * vec[2].len());
                        }
                        Vector v1 = ((Bezier3 *) s)->findTang(((Bezier3 *) s)->findT(vertex));
                        bezrad = v1.len2() / 2;
                    }
					
                    for (int i=0; i < (int)dirs.size(); i++) {
                        CurveAngle tmpca;
                        tmpca.dir = dirs[i];
                        tmpca.dir.normalize();
                        tmpca.c = s;
                        if (ca.size() == 0) {
                            tmpca.angle = 0.0;
                        } else {
                            tmpca.angle = anglebtw(ca[0].dir,dirs[i]);
                        }
                        ca.push_back(tmpca);
                    }
                }if (s->type == CURVE_NBEZIER) {
                    vector<Vector> dirs;
                    vector<Point> radPts;
                    NBezier * nb = (NBezier *) s;
                    nb->getOutput(2 * nb->degree, radPts);
                    foreach(Point pt, radPts) {
                        Vector dist = pt - pts->V[vertex];
                        if (dist.len() < radius && dist.len() > config->epsilon) {
                            radius = dist.len();
                        }
                    }
					
                    nb->findAngle(vertex, dirs);
                    for (int i=0; i < (int)dirs.size(); i++) {
                        CurveAngle tmpca;
                        tmpca.dir = dirs[i];
                        tmpca.dir.normalize();
                        tmpca.c = s;
                        if (ca.size() == 0) {
                            tmpca.angle = 0.0;
                        } else {
                            tmpca.angle = anglebtw(ca[0].dir,dirs[i]);
                        }
                        ca.push_back(tmpca);
                    }
                }
				
            }
			
            for (int ii=0; ii<(int)ca.size(); ii++) {
                // set the angles which are about 2*PI to be 0.0
                if (ca[ii].angle > 2*PI - config->epsilon) {
                    ca[ii].angle = 0.0;
                }
            }
			
            // sort ca list for insertion into small angle circle
            // bubble sort: inefficient, but negligible to overall performance.
            for (int i=0; i<(int)ca.size(); i++) {
                for (int j=0; j<(int)ca.size()-1; j++) {
					
                    if (ca[j].angle > ca[j+1].angle + config->epsilon ||
                   ( ca[j].angle > ca[j+1].angle - config->epsilon &&
                    ca[j].curvatureInfo > ca[j+1].curvatureInfo) ){
                      CurveAngle tmp = ca[j];
                      ca[j] = ca[j+1];
                      ca[j+1] = tmp;
                   }
                }
            }
			
            bool smallangle = false;
			
            for (int i=0; i<(int)ca.size()-1; i++) {
                if (ca[i+1].angle - ca[i].angle <= PI/3) {
                    smallangle = true;			    
                    break;
                }
            }
            if (!smallangle && ca.size() >= 2) {
                if (2*PI - ca[ca.size()-1].angle <= PI/3) {
                    smallangle = true;
                }
            }
			
            if (smallangle) {
                Vector v1, v2;
                double theta;
                //Check all curves not on the vertex for intersection with potential circle,
                foreach (Curve * s, curves) {				// and adjust radius accordingly
                    double rad;
                    if(pts->l2s[vertex].find(s) == pts->l2s[vertex].end()) {	
                        switch (s->type) {
                            case CURVE_SEGMENT:
                                v1 = pts->V[vertex] - pts->V[s->leftEndpoint];
                                v2 = pts->V[s->rightEndpoint] - pts->V[s->leftEndpoint];
                                theta = (PI/2) - min(anglebtw(v2, v1), anglebtw(v1, v2));
								
                                if (theta < 0) { // If the closest point on the segment is one of the endpoints
                                    v2 = pts->V[vertex] - pts->V[s->rightEndpoint];
                                    if (v1.len() < radius)	//Distance to leftEndpoint
                                        radius = v1.len();
                                    if (v2.len() < radius)	//Distance to rightEndpoint
                                        radius = v2.len();
                                } else {
                                    Vector v3 = pts->V[vertex] - pts->V[s->rightEndpoint];
                                    Vector v4 = pts->V[s->leftEndpoint] - pts->V[s->rightEndpoint];
                                    theta = (PI/2) - min(anglebtw(v3, v4), anglebtw(v4, v3));
                                    if (theta < 0) { // If the closest point on the segment is one of the endpoints
                                        v4 = pts->V[vertex] - pts->V[s->leftEndpoint];
                                        if (v3.len() < radius)	//Distance to rightEndpoint
                                            radius = v3.len();
                                        if (v4.len() < radius)	//Distance to leftEndpoint
                                            radius = v4.len();
                                    } else { 
                                        if (fabs(v3.len() * cos(theta)) < radius)
                                            radius = fabs(v3.len() * cos(theta));
                                    }
                                }
                                break;
                            case CURVE_CIRCLE:
                                v1 = ((Circle *) s)->center - pts->V[vertex];	//Vector towards center point and nearest point on circle.
                                if (fabs(v1.len() - ((Circle *) s)->radius) < radius)
                                    radius = fabs(v1.len() - ((Circle *) s)->radius);
                                break;
                            case CURVE_BEZIER:	//Tests using straight line approximations of the bezier similiar to line segments.
                                rad = ((Bezier *) s)->distToNearestPt(vertex);
                                if (radius > rad)
                                    radius = rad;
                                break;
                            case CURVE_CUBBEZ:
                                rad = ((Bezier *) s)->distToNearestPt(vertex);
                                if (radius > rad)
                                    radius = rad;
                                break;
                            case CURVE_NBEZIER:
                                rad = ((Bezier *) s)->distToNearestPt(vertex);
                                if (radius > rad)
                                    radius = rad;
                                break;
                            case CURVE_ARC:
                                v1 = ((Arc *) s)->center - pts->V[vertex];	//Nearest intersection point of circle
                                double distance = v1.len() - ((Arc *) s)->radius;
                                v1.normalize();
                                Point p = pts->V[vertex] + distance * v1;
                                if (((Arc *) s)->isOn(p)) {
                                    if (v1.len() < radius)
                                        radius = v1.len();
                                }
								
                                break;
                        }
                    }
                }
				
                //Test the distance to each input point in the mesh.
                for (int vertex1 = 0; vertex1 < (int)pts->V.size(); vertex1++) {
                    if (vertex1 != vertex) {
                        v1 = pts->V[vertex1] - pts->V[vertex];
                        if (v1.len() < radius) {
                            radius = v1.len();
                        }
                    }
                }
				
                radius/=3.0;
                vertex2Radius[vertex] = radius;
            }
        }
    }
    PrintEnd("FixSmallAngles()");
}

/*
 Adds a small angle protection circle to the mesh at
 the point given by vertex, with a radius given by
 the parameter radius.
 
 Arguments below are:
 int vertex: center point of the small angle circle.
 double radius: radius of the small angle circle.
 */
void PSC::addSmallAngleCircle(int vertex, double radius) {
    protectedVertices.insert(vertex);
    Circle * newCir = new Circle(this, pts->V[vertex], radius);
    ostringstream oss1;
    oss1 << "SAcircle-" << pts->V[vertex];
    newCir->name = oss1.str();
    vector<SmallAnglePt> insertPts;
    Point reference;
	
    foreach(Curve * c, pts->l2s[vertex]) {
        vector<Point> temp;
        c->SplitAtDist(vertex, radius, temp);
        foreach(Point p, temp) {
            int loc = pts->AddPoint(p, POINT_INPUT);
            presplitVertices.push_back(loc);
            SmallAnglePt newpt;
            newpt.pt = loc;
            if (insertPts.empty()) {
                newpt.dist = 0;
            }else {
                newpt.dist = newCir->anglebtw(pts->V[insertPts[0].pt], p);
            }
			
            if (c->left[vertex] == -1) {
                c->InsertVertex(loc, vertex);
            }
            else if(c->right[vertex] == -1) {
                c->InsertVertex(loc, c->left[vertex]);
            }
            else if(d2(p, pts->V[c->left[vertex]]) < d2(pts->V[vertex], pts->V[c->left[vertex]])) {
                c->InsertVertex(loc, c->left[vertex]);
            }
            else {
                c->InsertVertex(loc, vertex); 
            }
			
            insertPts.push_back(newpt);
        }
    }
	
    //Bubble Sort
    for (int i=0; i<(int)insertPts.size(); i++) {
        for (int j=0; j<(int)insertPts.size()-1; j++) {
            if (insertPts[j].dist > insertPts[j+1].dist) {
                SmallAnglePt tmp = insertPts[j];
                insertPts[j] = insertPts[j+1];
                insertPts[j+1] = tmp;
            }
        }
    }
	
    // new vertices have been added to existing curves...
    // then add them to the new circle
    newCir->InsertVertex(insertPts[0].pt);
    updateTriangulation(insertPts[0].pt, newCir);
	
    for(int i = 1;i < insertPts.size();i++) {
        newCir->InsertVertex(insertPts[i].pt, insertPts[i-1].pt);
        updateTriangulation(insertPts[i].pt, newCir);
    }
    curves.push_back(newCir);
    vector<int> newPts;
    newCir->PreSplit(newPts);
	
    foreach(int p, newPts) {
        newCir->InsertVertex(p);
        updateTriangulation(p, newCir);
    }
}

/*
 New vertex has been inserted into the triangulation on newCir.
 Then the queue of bad triangles and segments needs to be updated
 by removing old subsegments and triangles and adding the newly
 encroached subsegments to the mesh as well as newly created
 bad triangles.
 
 Arguments below are:
 int p: the newly inserted point.
 Circle * newCir: the new small angle circle being added to the triangulation.
 */
void PSC::updateTriangulation(int p, Circle * newCir) {
    vector<Triangle *> oldcav;
    vector<Triangle *> newcav;
    // 1. Insert vertex into triangulation
    tri->InsertPoint(p,oldcav,newcav);
	
    // 2. Delete oldcav triangles from the queue
    foreach (Triangle * t, oldcav) {
        set<Triangle *>::iterator it = triangleQueue.find(t);
        if (it != triangleQueue.end()) {
            triangleQueue.erase(it);
        }
		
        // return the memory
        triPool->DeleteTriangle(t);
    }
	
	
    // 3.  Enqueue newcav triangles as necessary
    foreach (Triangle * t, newcav) {
        if (geom->RadiusEdgeRatio(t) > config->maxradiusedgeratio && !IsProtected(t) ){
            triangleQueue.insert(t);
        }
    }
	
    // Old segments do not have to be removed, since they are stored
    // as only a leftEndpoint and a curve. However, since they may
    // no longer be encroached, this must be checked and updated as necessary.
    // Since this also requires iterating through the whole queue is easier to just
    // remove newly unencroached segments and add new encroached segments
	
    // 4. Removing no longer encroached segments from the queue
    foreach (Curve * oldSeg, pts->l2s[p]) {
        Subsegment removed;
        removed.inputCurve = oldSeg;
        removed.leftEndpoint = oldSeg->left[p];
		
        set<Subsegment>::iterator it = segmentQueue.find(removed);
        if (it != segmentQueue.end()) {
            segmentQueue.erase(it);
        }		
    }
	
    // 5.  Queue existing segments encroached by new point
    vector<Subsegment> segmentsEncroached;
    vector<int> nbs;
    tri->Neighbors(p, nbs);
    if (DoesEncroach(pts->V[p],nbs,segmentsEncroached)) {
        foreach(Subsegment ss, segmentsEncroached) {
            if (!IsProtected(ss)) {
                segmentQueue.insert(ss);
            }
        }
    }
	
    // 6. Enqueue new segments if necessary
    foreach(Curve * s, pts->l2s[p]) {
        Subsegment ss1, ss2;		
		
        ss1.inputCurve = s;
        ss1.leftEndpoint = p;
        if (!IsProtected(ss1) && IsEncroached(ss1)) {
            segmentQueue.insert(ss1);      
        }

        if (s->left[p] != -1) {
            ss2.inputCurve = s;
            ss2.leftEndpoint = s->left[p];
            if (!IsProtected(ss2) && IsEncroached(ss2)) {
                segmentQueue.insert(ss2);      
            }
        }
    }
}

/*
 Returns angle between both vectors in radians.
 
 Arguments below are:
 Vector l = the left vector
 Vector r = the right vector
 */
double PSC::anglebtw1(Vector l, Vector r) {
    l.normalize();
    r.normalize();
    double dotprod = l*r;
    if (dotprod > 1.0) {
        dotprod = 1.0;
    } else if (dotprod < -1.0) {
        dotprod = -1.0;
    }
    return acos(dotprod);
}

/*
 Returns angle between both vectors in radians.
 Vector l = the left vector
 Vector r = the right vector
 */
double PSC::anglebtw(Vector l, Vector r) {
    l.normalize();
    r.normalize();
	
    if(l == r)		//360 degree angle
        return 2*PI;
	
    double crossprod = r.xy[0] * l.xy[1] - r.xy[1] * l.xy[0];
    if (crossprod > 1.0) {
        crossprod = 1.0;
    } else if (crossprod < -1.0) {
        crossprod = -1.0;
    }
    double sinangle = asin(crossprod);  //Arcsin of magnitude of unit vector cross product
	
    double dotprod = l*r;
    if (dotprod > 1.0) {
        dotprod = 1.0;
    } else if (dotprod < -1.0) {
        dotprod = -1.0;
    }
    double cosangle = acos(dotprod);
	
	
    // 1st quaadrant
    if (sinangle > 0.0 && cosangle <PI/2.0) {
        return cosangle;
    }
    // 2nd quadrant
    if (sinangle > 0.0) {
        return cosangle;
    }
    // 3rd quadrant
    if (cosangle >= PI/2.0) {
        return PI - sinangle;
    }
    // 4th quadrant
    return 2*PI + sinangle;
}

/*
 Initialize the bounding box for the mesh.
 */
void PSC::InitBoundingBox() {
    Point center((tri->upperbound[0]+tri->lowerbound[0])/2.0,
                (tri->upperbound[1]+tri->lowerbound[1])/2.0);
	
	
    tri->center = center;
	
    double bboxsize = max(tri->upperbound[0]-tri->lowerbound[0],
                         tri->upperbound[1]-tri->lowerbound[1])/2.0;
	
    tri->bbox = bboxsize * config->boundingboxratio;
	
    tri->Initialize(); // set up the initial square
	
    tri->bbox = tri->bbox + config->epsilon;  // create a little wiggle room... 
	
    if (logfileisopen) {
        (*LOGFILE) << "Bounding box center: " << tri->center << endl;
        (*LOGFILE) << "Bounding box size:   " << tri->bbox << endl;
    }
	
    PrintEnd("PreSplit(): Reduce all curves to 90 degrees");
}

/*
 Determines if a given point is protected.
 
 Arguments below are:
 int p: the point.
 */
bool PSC::IsProtected(int p) {
    if (protectedVertices.find(p) != protectedVertices.end() ) {
        return true;
    }
    return false;
}

/*
 Determines if any point from a given subsegment is protected.
 
 Arguments below are:
 Subsegment ss: the subsegment.
 */
bool PSC::IsProtected(Subsegment ss) {
    if (IsProtected(ss.leftEndpoint) ||
       IsProtected(ss.inputCurve->right[ss.leftEndpoint])) {
        return true;
    }
    return false;
}

/*
 Determines if any vertex of a given triangle is protected.
 
 Arguments below are:
 Triangle * t: the triangle.
 */
bool PSC::IsProtected(Triangle * t) {
    if (protectedVertices.find(t->v[0])!= protectedVertices.end() || 
       protectedVertices.find(t->v[1])!= protectedVertices.end() || 
       protectedVertices.find(t->v[2])!= protectedVertices.end()) {
        return true;
    }
    return false;
}

/*
 Determines if a triangle is part of the exterior
 mesh. That is, if that triangle is not contained
 fully by any of the input features.
 
 Arguments below are:
 Triangle * t: the triangle.
 */
bool PSC::IsExterior(Triangle * t) {
    if (exteriorTriangles.find(t) == exteriorTriangles.end()){
        return false;
    }
    return true;
}

/*
 Determines if a triangle is enclosed in the interior
 of a small angle circle.
 
 Arguments below are:
 Triangle * t: the triangle.
 */
bool PSC::IsInterior(Triangle * t) {
    if (smallAngleTriangles.find(t) == smallAngleTriangles.end()){
        return false;
    }
    return true;
}


/*
 Adds triangles that are outside of all input features
 to the set exteriorTriangles so that they will not be
 output in printSVG.
 */
void PSC::SetUpExteriorTriangles() {
	
    Triangle * t_start = tri->p2t[bboxVertices[0]];
	
    exteriorTriangles.insert(t_start);
	
    set<Triangle *> TQ; // queue of triangles
	
    // breadth-first search through the exterior triangulation
    // to fill the exteriorTriangles set
	
    if (t_start->n[0] != NULL) {
        TQ.insert(t_start->n[0]);
    }
    if (t_start->n[1] != NULL) {
        TQ.insert(t_start->n[1]);
    }
    if (t_start->n[2] != NULL) {
        TQ.insert(t_start->n[2]);
    }
	
    while (TQ.size() > 0) {
		
        // pop the first element from the queue
        set<Triangle *>::iterator it = TQ.begin();
        Triangle * t_cur = *it;
        TQ.erase(it);    
		
        bool markTCur = false;
        //If a nieghbor is exterior, search through all curves that intersect with the shared vertices
        // of that triangle and check to see if a curve is on that edge, i.e. the triangle crosses a boundary
        // excempting curves that are small angle circles.
		
        if (IsExterior(t_cur->n[0])) {
            int a = t_cur->v[1];
            int b = t_cur->v[2];
            markTCur = true;
            foreach (Curve * curve, pts->l2s[a]) {
                if ((curve->left[a] == b || curve->right[a] == b) &&
                !(curve->type == CURVE_CIRCLE && curve->name.substr(0, 8) == "SAcircle")) {
                    markTCur = false;
                }
            }    
        }
        if (markTCur == false && IsExterior(t_cur->n[1])) {
            int a = t_cur->v[0];
            int b = t_cur->v[2];
            markTCur = true;
            foreach (Curve * curve, pts->l2s[a]) {
                if ((curve->left[a] == b || curve->right[a] == b) &&
                    !(curve->type == CURVE_CIRCLE && curve->name.substr(0, 8) == "SAcircle")) {
                    markTCur = false;
                }
            }
        }
        if (markTCur == false && IsExterior(t_cur->n[2])) {
            int a = t_cur->v[0];
            int b = t_cur->v[1];
            markTCur = true;
            foreach (Curve * curve, pts->l2s[a]) {
                if ((curve->left[a] == b || curve->right[a] == b) &&
                !(curve->type == CURVE_CIRCLE && curve->name.substr(0, 8) == "SAcircle")) {
                    markTCur = false;
                }
            }
        }
        //If at least one neighbor did not cross a boundary and is exterior, add to set
        if (markTCur == true) {
            exteriorTriangles.insert(t_cur);
			
            if (t_cur->n[0] != NULL && !IsExterior(t_cur->n[0]) ) {
                TQ.insert(t_cur->n[0]);
            }
            if (t_cur->n[1] != NULL && !IsExterior(t_cur->n[1]) ) {
                TQ.insert(t_cur->n[1]);
            }
            if (t_cur->n[2] != NULL && !IsExterior(t_cur->n[2]) ) {
                TQ.insert(t_cur->n[2]);
            }
			
        }
		
    }
	
}

/*
 Adds triangles that are inside of circles created
 by small angles to the set smallAngleTriangles so
 that they will not be output in printSVG.
 */
void PSC::SetUpSmallAngleTriangles() {
    Triangle * t_start;
    set<Triangle *> TQ; // queue of triangles
	
    foreach(int vert, protectedVertices) {
        t_start = tri->p2t[vert];
		
        smallAngleTriangles.insert(t_start);
        // breadth-first search  through the interior triangles
        // to fill the smallAngleTriangles set
		
        //Insert neighboring triangles into queue
        if (t_start->n[0] != NULL) {
            TQ.insert(t_start->n[0]);
        }
        if (t_start->n[1] != NULL) {
            TQ.insert(t_start->n[1]);
        }
        if (t_start->n[2] != NULL) {
            TQ.insert(t_start->n[2]);
        }
		
        while (TQ.size() > 0) {
            // pop the first element from the queue
			
            set<Triangle *>::iterator it = TQ.begin();
            Triangle * t_cur = *it;
            TQ.erase(it);    
			
            bool markTCur = false;
            //If a nieghbor is interior, search through all curves that intersect with the shared vertices
            // of that triangle and check to see if a curve is on that edge, i.e. the triangle crosses a boundary
            // and that the curve is a smallAngle Circle
            if (IsInterior(t_cur->n[0])) {
                int a = t_cur->v[1];
                int b = t_cur->v[2];
                markTCur = true;
                foreach (Curve * curve, pts->l2s[a]) {
                    if (curve->type == CURVE_CIRCLE && ((Circle *)curve)->center == pts->V[vert]) {
                        if (curve->left[a] == b || curve->right[a] == b) {
                            markTCur = false;
                        }
                        break;
                    }
                }    
            }
            if (markTCur == false && IsInterior(t_cur->n[1])) {
                int a = t_cur->v[0];
                int b = t_cur->v[2];
                markTCur = true;
                foreach (Curve * curve, pts->l2s[a]) {
                    if (curve->type == CURVE_CIRCLE && ((Circle *)curve)->center == pts->V[vert]) {
                        if (curve->left[a] == b || curve->right[a] == b) {
                            markTCur = false;
                        }
                        break;
                    }
                }
            }
            if (markTCur == false && IsInterior(t_cur->n[2])) {
                int a = t_cur->v[0];
                int b = t_cur->v[1];
                markTCur = true;
                foreach (Curve * curve, pts->l2s[a]) {
                    if (curve->type == CURVE_CIRCLE && ((Circle *)curve)->center == pts->V[vert]) {
                        if (curve->left[a] == b || curve->right[a] == b) {
                            markTCur = false;
                        }
                        break;
                    }
                }
            }
            //Insert triangle as long as at least one edge does not share a curve edge.
            if (markTCur == true) {
                smallAngleTriangles.insert(t_cur);
				
                if (t_cur->n[0] != NULL && !IsInterior(t_cur->n[0]) ) {
                    TQ.insert(t_cur->n[0]);
                }
                if (t_cur->n[1] != NULL && !IsInterior(t_cur->n[1]) ) {
                    TQ.insert(t_cur->n[1]);
                }
                if (t_cur->n[2] != NULL && !IsInterior(t_cur->n[2]) ) {
                    TQ.insert(t_cur->n[2]);
                }
            }
        }	
    }
}
