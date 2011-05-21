#include "allh.h"

#define STAR_RESERVE_SIZE 20

// faster version does not use set data structure...
void Triangulation::Neighbors(int p, vector < int >&nbs) {
    vector <Triangle * >tocheck;
    vector <Triangle *> checked2;
	
    tocheck.reserve(STAR_RESERVE_SIZE);
    checked2.reserve(STAR_RESERVE_SIZE);
	
    Triangle * startingTri = p2t[p];
	
    tocheck.push_back(startingTri);
    checked2.push_back(startingTri);
	
    startingTri->visited = 1;
	
    Points * pts = c->pts;
	
    while (tocheck.size() > 0) {
		
        Triangle *myt = tocheck[tocheck.size() - 1];
        tocheck.pop_back();
		
        for (int i=0; i<3; i++) {
            if (myt->v[i] != p && pts->found[myt->v[i]] == 0) {
                nbs.push_back(myt->v[i]);
                pts->found[myt->v[i]] = 1;
            }
        }
		
        for (int i = 0; i<3; i++) {
            if (myt->v[i] != p && myt->n[i] != NULL) {
                if (myt->n[i]->visited == 0) {
                    tocheck.push_back(myt->n[i]);
                    checked2.push_back(myt->n[i]);
                    myt->n[i]->visited = 1;
                }
            }
        }
    }
	
    for (int i=0; i<(int)checked2.size(); i++) {
    checked2[i]->visited = 0;
    }
	
    for (int i=0; i<(int)nbs.size(); i++) {
        pts->found[nbs[i]] = 0;
    }   
}

// FIXME: make efficient in 2d some day...
void Triangulation::Neighbors(int i, set < int >&nbs) {
	
    set<Triangle *> star;
	
    Triangles(i,star);
	
    foreach(Triangle * t, star) {
        if (t->v[0] != i) {
            nbs.insert(t->v[0]);
        }
        if (t->v[1] != i) {
            nbs.insert(t->v[1]);
        }
        if (t->v[2] != i) {
            nbs.insert(t->v[2]);
        }
    }
	
	
	/*
	 
	 Triangle * t = p2t[i];
	 
	 int start=-1,end=-1,index1=-1,index=-1;
	 if (i == t->v[0]) {
	 start = t->v[1];
	 end = t->v[2];
	 index = 2;
	 } else if (i == t->v[1]) {
	 start = t->v[2];
	 end = t->v[0];
	 index = 0;
	 } else if (i == t->v[2]) {
	 start = t->v[0];
	 end = t->v[1];
	 index = 1;
	 }
	 int curv = start;
	 Triangle * curt = t->n[index];
	 
	 // FIXME: issue with boundary points...
	 
	 while (curt != t) {
	 nbs.insert(curv);
	 
	 if (i == curt->v[0]) {
	 if (curv == curt->v[1]) {
	 index = 1;
	 index1 = 2;
	 } else if (curv == curt->v[2]) {
	 index = 2;
	 index1 = 1;
	 }
	 } else if (i == curt->v[1]) {
	 if (curv == curt->v[0]) {
	 index = 0;
	 index1 = 2;
	 } else if (curv == curt->v[2]) {
	 index = 2;
	 index1 = 0;
	 }
	 } else if (i == curt->v[2]) {
	 if (curv == curt->v[0]) {
	 index = 0;
	 index1 = 1;
	 } else if (curv == curt->v[1]) {
	 index = 1;
	 index1 = 0;
	 }
	 }
	 
	 curt = curt->n[index1];
	 
	 }
	 */
	
}

// faster version...
void Triangulation::Triangles(int p, vector < Triangle * >&star) {
	
    vector < Triangle * >tocheck;
	
    star.reserve(STAR_RESERVE_SIZE);
    tocheck.reserve(STAR_RESERVE_SIZE);   
	
    if (p2t.find(p) == p2t.end()) {
        return;
    }
	
    Triangle * startingTri = p2t[p];
	
    tocheck.push_back(startingTri);
    star.push_back(startingTri);
	
    startingTri->visited = 1;
	
    while (tocheck.size() > 0) {
		
        Triangle *myt = tocheck[tocheck.size() - 1];
        tocheck.pop_back();
		
        for (int i = 0; i < 3; i++) {
            if (myt->v[i] != p && myt->n[i] != NULL && myt->n[i]->visited == 0){
                tocheck.push_back(myt->n[i]);
                star.push_back(myt->n[i]);
                myt->n[i]->visited = 1;
            }
			
        }
    } 
	
    int len = star.size();
    for (int i=0; i<len; i++) {
        star[i]->visited = 0;
    }
}



// FIXME: speed this up for 2d some day...
void Triangulation::Triangles(int p, set < Triangle * >&star) {
	
    vector < Triangle * >tocheck;
	
    if (p2t.find(p) == p2t.end()) {
        return;
    }
	
    tocheck.push_back(p2t[p]);
	
    star.insert(p2t[p]);
	
    while (tocheck.size() > 0) {
		
        Triangle *myt = tocheck[tocheck.size() - 1];
        tocheck.pop_back();
		
        for (int i = 0; i < 3; i++) {
            if (myt->n[i] != NULL && myt->n[i]->contains(p)
                && star.find(myt->n[i]) == star.end()) {
                tocheck.push_back(myt->n[i]);
                star.insert(myt->n[i]);
            }
        }
    }
	
}

Triangle * Triangulation::Locate(Point p) {
    return Locate(p,p2t.begin()->first);
}

Triangle * Triangulation::Locate(int i) {
    return Locate(i,p2t.begin()->first);
}

Triangle * Triangulation::Locate(Point p, int hint) {
	
    int hintold = -1;
    double dhint, di;
    dhint = d2(p, c->pts->V[hint]);
    vector<int> nbs;
	
    Neighbors(hint, nbs);
	
    while (hintold != hint) {
        hintold = hint;
        foreach(int i, nbs) {
            di = d2(p, c->pts->V[i]);
            if (di < dhint - c->config->epsilon) {
                dhint = di;
                hint = i;
            }
        }
        nbs.clear();
        Neighbors(hint, nbs);
		
    }
	
    vector<Triangle *> star;
    Triangles(hint, star);
	
    foreach (Triangle * t, star) {
        if (c->geom->inSphereStrict(p,t)) {
            return t;
        }
    }
	
    cout <<"Triangulation::Locate() - ERROR: Locate failed." << endl;
    cout << "p=" << p << ", star.size()=" << star.size() 
    << ", hint=" << c->pts->V[hint] << endl;
    foreach(Triangle * t, star) {
        cout << t << " " << c->pts->V[t->v[0]] << " "
        << c->pts->V[t->v[1]] << " "<< c->pts->V[t->v[2]] << endl;
    }
	
	
    return NULL;
    
}

Triangle * Triangulation::Locate(int i, int hint) {
    return Locate(c->pts->V[i],hint);
}

// s is a starting simplex...
Triangle * Triangulation::Locate(Point p, Triangle * t) {
    return Locate(p,t->v[0]);
}

Triangle * Triangulation::Locate(int i, Triangle * t) {
    return Locate(i,t->v[1]);
}

void Triangulation::InsertPoint(int i) {
	
    vector<Triangle *> oldcavity;
    vector<Triangle *> newcavity;
	
    InsertPoint(i, Locate(i),oldcavity, newcavity);
}

void Triangulation::InsertPoint(int i, int hint) {
    set<Triangle *> oldcavity;
    set<Triangle *> newcavity;
	
    InsertPoint(i, Locate(i,hint),oldcavity, newcavity);
}

void Triangulation::InsertPoint(int i, Triangle * t) {
    set<Triangle *> oldcavity;
    set<Triangle *> newcavity;
	
    InsertPoint(i, t,oldcavity, newcavity);
}
void Triangulation::InsertPoint(int i, set < Triangle * >&oldcavity,
                                set < Triangle * >&newcavity) {
    InsertPoint(i, Locate(i),oldcavity, newcavity);
}

void Triangulation::InsertPoint(int i, vector < Triangle * >&oldcavity,
                                vector < Triangle * >&newcavity) {
    InsertPoint(i, Locate(i),oldcavity, newcavity);
}


void Triangulation::InsertPoint(int i, int hint, 
                                set < Triangle * >&oldcavity,
                                set < Triangle * >&newcavity) {
    InsertPoint(i, Locate(i,hint),oldcavity, newcavity);
}



void Triangulation::InsertPoint(int i, int hint, 
                                vector < Triangle * >&oldcavity,
                                vector < Triangle * >&newcavity) {
    InsertPoint(i, Locate(i,hint),oldcavity, newcavity);
}


void Triangulation::InsertPoint(int i, Triangle * t, 
                                set < Triangle * >&oldcavity,
                                set < Triangle * >&newcavity) {
	
    //cavity
    set < Face > bndry;
    Cavity(c->pts->V[i], t, oldcavity, bndry);
	
    // stitch
    Stitch(i, bndry, newcavity);
	
    // fix p2t 
	
    foreach(Triangle * t1, newcavity) {
        SetP2t(t1);
    }
	
    np++;
	
    nt += newcavity.size() - oldcavity.size();
	
}

void Triangulation::InsertPoint(int i, Triangle * t, 
                                vector < Triangle * >&oldcavity,
                                vector < Triangle * >&newcavity) {
    //get cavity
    vector < Face > bndry;
    Cavity(c->pts->V[i], t, oldcavity, bndry);
	
    // stitch new triangles
    Stitch(i, bndry, newcavity);
	
    // update info
    foreach(Triangle * t1, newcavity) {
        SetP2t(t1);
    }
    np++;
    nt += newcavity.size() - oldcavity.size();
}


void Triangulation::Cavity(Point & p, Triangle * t, vector < int >&nbs) {
	
    vector < Triangle * > cavity;
    vector < Face >  bndry;
	
    Cavity(p, t, cavity, bndry);
	
    foreach(Triangle * t1, cavity) {
        if (c->pts->found[t1->v[0]] == 0) {
            nbs.push_back(t1->v[0]);
            c->pts->found[t1->v[0]] = 1;
        }
        if (c->pts->found[t1->v[1]] == 0) {
            nbs.push_back(t1->v[1]);
            c->pts->found[t1->v[1]] = 1;
        }
        if (c->pts->found[t1->v[2]] == 0) {
            nbs.push_back(t1->v[2]);
            c->pts->found[t1->v[2]] = 1;
        }
    }
	
    foreach(int i, nbs) {
        c->pts->found[i] = 0;
    }
	
}



void Triangulation::Cavity(Point & p, Triangle * t, set < int >&nbs) {
	
    set < Triangle * > cavity;
    set < Face >  bndry;
	
    Cavity(p, t, cavity, bndry);
	
    foreach(Triangle * t1, cavity) {
        nbs.insert(t1->v[0]);
        nbs.insert(t1->v[1]);
        nbs.insert(t1->v[2]);
    }
	
}


// grab the cavity and it's boundary...
void Triangulation::Cavity(int p, Triangle * t, vector < Triangle * >&cavity,
                           vector < Face > &bndry)
{
    return Cavity(c->pts->V[p], t, cavity, bndry);
}


// Cavity
//
// This function returns the cavity and the bounary of the cavity associated with
// point p (which is not necessarily in the mesh yet...
//
//
void Triangulation::Cavity(Point & p, Triangle * t, vector < Triangle * >&cavity,
                           vector < Face > &bndry)
{
    vector < Triangle * >tocheck;
	
    tocheck.reserve(STAR_RESERVE_SIZE);
    cavity.reserve(STAR_RESERVE_SIZE);
    bndry.reserve(STAR_RESERVE_SIZE);
	
    tocheck.push_back(t);
    cavity.push_back(t);
	
    t->visited = 1;
	
    while (tocheck.size() > 0) {
		
        Triangle *myt = tocheck[tocheck.size() - 1];
        tocheck.pop_back();
		
        for (int i = 0; i < 3; i++) {
            if (myt->n[i] == NULL) {
				
                // create a dummy simplex representing the boundary...
                int a, b;
                a = myt->v[0];
                b = myt->v[1];
                if (i == 0) {
                    a = myt->v[2];
                } else if (i == 1) {
                    b = myt->v[2];
                }
				
                //start here...
                Triangle *newt = c->triPool->NewTriangle(a, b, -1);
				
                Face pr;
                pr.s = newt;
                pr.v = -1;
                // no need to link this simplex to anything...
                bndry.push_back(pr);
				
            } else {
                if (c->geom->inSphereStrict(p, myt->n[i])) {
                    if (myt->n[i]->visited == 0) {
                        // add neighbor to the cavity and queue to check...
                        cavity.push_back(myt->n[i]);
                        myt->n[i]->visited = 1;
                        tocheck.push_back(myt->n[i]);
                    }
                } else {
                    // this is on the boundary
                    Face pr;
                    pr.s = myt->n[i];
                    for (int k = 0; k < 3; k++) {
                        if (myt->n[i]->n[k] == myt) {
                            pr.v = k;
                        }
                    }
                    bndry.push_back(pr);
					
                }
            }
        }
        }
	
    for (int i=0; i<(int)cavity.size(); i++) {
        cavity[i]->visited = 0;
    }
	
    return;
}


void Triangulation::Cavity(Point & p, Triangle * t, 
                           set < Triangle * >&cavity,
                           set < Face > &bndry) {
    vector < Triangle * >tocheck;
	
    tocheck.push_back(t);
	
    cavity.insert(t);
	
    while (tocheck.size() > 0) {
		
        Triangle *myt = tocheck[tocheck.size() - 1];
        tocheck.pop_back();
		
		
        for (int i = 0; i < 3; i++) {
			
            if (myt->n[i] == NULL) {
				
                // create a dummy simplex representing the boundary...
				
                int a, b;
                a = myt->v[0];
                b = myt->v[1];
                if (i == 0) {
                    a = myt->v[2];
                } else if (i == 1) {
                    b = myt->v[2];
                }
				
                //start here...
                Triangle *newt = c->triPool->NewTriangle(a, b, -1);
				
                Face pr;
                pr.s = newt;
                pr.v = -1;
                // no need to link this simplex to anything...
                bndry.insert(pr);
				
            } else {
                if (c->geom->inSphereStrict(p, myt->n[i])) {
                    if (cavity.find(myt->n[i]) == cavity.end()) {
                        // add neighbor to the cavity and queue to check...
                        cavity.insert(myt->n[i]);
                        tocheck.push_back(myt->n[i]);
                    }
                } else {
                    // this is on the boundary
                    Face pr;
                    pr.s = myt->n[i];
                    for (int k = 0; k < 3; k++) {
                        if (myt->n[i]->n[k] == myt) {
                            pr.v = k;
                        }
                    }
                    bndry.insert(pr);
                }
            }
        }
    }

    return;
	
}



// Stitch
//
// This functions connects point number i to the boundary faces
// and returns the new simplices in newcavity.
//
// Important: does not update p2s!
//
void Triangulation::Stitch(int i, vector < Face > &bndry, vector < Triangle * >&newcavity)
{
	
    newcavity.reserve(STAR_RESERVE_SIZE);
	
    Triangle *toHook = NULL;
	
    map < int, Triangle * >toBeStitched2d;
    // make each simplex and hook it to the outside of the cavity...
    foreach(Face pr, bndry) {
        Triangle *news;
		
        if (pr.v == -1) {
            // turn the dummy simplex into the simplex to insert...
            pr.s->v[2] = i;
            pr.s->n[2] = NULL; // the opposite face is outside the mesh...
            news = pr.s;
            newcavity.push_back(pr.s);
        } else {
            int a, b;
            a = pr.s->v[0];
            b = pr.s->v[1];
            if (pr.v == 0) {
                a = pr.s->v[2];
            } else if (pr.v == 1) {
                b = pr.s->v[2];
            }
			
            news = c->triPool->NewTriangle(a, b, i);
			
            news->n[2] = pr.s;
            newcavity.push_back(news);
            pr.s->n[pr.v] = news;
        }
        SetRadiusAndCC(news);
		
        for (int i = 0; i < 2; i++) {
            map < int, Triangle * >::iterator mIt =
            toBeStitched2d.find(news->v[i]);
            if (mIt == toBeStitched2d.end()) {
                toBeStitched2d[news->v[i]] = news;
            } else {
                int k = -1, k1 = -1; // k will be set to be the missing value from i and j...
                if (i != 1) {
                    k = 1;
                } else {
                    k = 0;
                }
                toHook = mIt->second;
                if (toHook->v[0] == news->v[i]) {
                    k1 = 1;
                } else {
                    k1 = 0;
                }
                toHook->n[k1] = news;
                news->n[k] = toHook;
                toBeStitched2d.erase(mIt);
            }
        }
    }
	
    return;
	
}

void Triangulation::Stitch(int i, set < Face > &bndry, set < Triangle * >&newcavity)
{
	
    Triangle *toHook = NULL;
	
	
    map < int, Triangle * >toBeStitched2d;
    foreach(Face pr, bndry) {
        Triangle *news;
		
        if (pr.v == -1) {
            // turn the dummy simplex into the simplex to insert...
            pr.s->v[2] = i;
            pr.s->n[2] = NULL; // the opposite face is outside the mesh...
            news = pr.s;
            newcavity.insert(pr.s);
        } else {
            int a, b;
            a = pr.s->v[0];
            b = pr.s->v[1];
            if (pr.v == 0) {
                a = pr.s->v[2];
            } else if (pr.v == 1) {
                b = pr.s->v[2];
            }
			
            news = c->triPool->NewTriangle(a, b, i);
			
            news->n[2] = pr.s;
            newcavity.insert(news);
            pr.s->n[pr.v] = news;
        }
        SetRadiusAndCC(news);
		
        for (int i = 0; i < 2; i++) {
            map < int, Triangle * >::iterator mIt = toBeStitched2d.find(news->v[i]);
            if (mIt == toBeStitched2d.end()) {
                toBeStitched2d[news->v[i]] = news;
            } else {
                int k = -1, k1 = -1; // k will be set to be the missing value from i and j...
                if (i != 1) {
                    k = 1;
                } else {
                    k = 0;
                }
                toHook = mIt->second;
                if (toHook->v[0] == news->v[i]) {
                    k1 = 1;
                } else {
                    k1 = 0;
                }
                toHook->n[k1] = news;
                news->n[k] = toHook;
                toBeStitched2d.erase(mIt);
            }
        }
    }
	
}



void Triangulation::Initialize() {
	
    Point p1(center.xy[0]+bbox,center.xy[1]+bbox);
    Point p2(center.xy[0],center.xy[1]+bbox);
    Point p3(center.xy[0]-bbox,center.xy[1]+bbox);
    Point p4(center.xy[0]-bbox,center.xy[1]);
    Point p5(center.xy[0]-bbox,center.xy[1]-bbox);
    Point p6(center.xy[0],center.xy[1]-bbox);
    Point p7(center.xy[0]+bbox,center.xy[1]-bbox);
    Point p8(center.xy[0]+bbox,center.xy[1]);
	
    int l1 = c->pts->AddPoint(p1,POINT_BBOX);
    int l2 = c->pts->AddPoint(p2,POINT_BBOX);
    int l3 = c->pts->AddPoint(p3,POINT_BBOX);
    int l4 = c->pts->AddPoint(p4,POINT_BBOX);
    int l5 = c->pts->AddPoint(p5,POINT_BBOX);
    int l6 = c->pts->AddPoint(p6,POINT_BBOX);
    int l7 = c->pts->AddPoint(p7,POINT_BBOX);
    int l8 = c->pts->AddPoint(p8,POINT_BBOX);
	
    c->bboxVertices.push_back(l1);
    c->bboxVertices.push_back(l2);
    c->bboxVertices.push_back(l3);
    c->bboxVertices.push_back(l4);
    c->bboxVertices.push_back(l5);
    c->bboxVertices.push_back(l6);
    c->bboxVertices.push_back(l7);
    c->bboxVertices.push_back(l8);
	
    Triangle *t1 = c->triPool->NewTriangle(l8,l1,l2);
    Triangle *t2 = c->triPool->NewTriangle(l2,l3,l4);
    Triangle *t3 = c->triPool->NewTriangle(l4,l5,l6);
    Triangle *t4 = c->triPool->NewTriangle(l6,l7,l8);
    Triangle *t5 = c->triPool->NewTriangle(l2,l4,l6);
    Triangle *t6 = c->triPool->NewTriangle(l6,l8,l2);
	
    // link the triangles...
	
    t1->link(NULL,t6,NULL);
    t2->link(NULL,t5,NULL);
    t3->link(NULL,t5,NULL);
    t4->link(NULL,t6,NULL);
    t5->link(t3,t6,t2);
    t6->link(t1,t5,t4);
	
    SetP2t(t1);
    SetP2t(t2);
    SetP2t(t3);
    SetP2t(t4);
    SetP2t(t5);
    SetP2t(t6);
	
    SetRadiusAndCC(t1);
    SetRadiusAndCC(t2);
    SetRadiusAndCC(t3);
    SetRadiusAndCC(t4);
    SetRadiusAndCC(t5);
    SetRadiusAndCC(t6);
	
    np = 8;
    nt = 6;
	
}

void Triangulation::SetRadiusAndCC(Triangle * t) {
	
    t->cc = c->geom->circumcenter(t);
    t->rsq = d2(t->cc,c->pts->V[t->v[0]]);
}


void Triangulation::SetP2t(Triangle * t)
{
    for (int i = 0; i < 3; i++) {
        p2t[t->v[i]] = t;
    }
}

Triangulation::titerator::titerator() {
	
}

Triangulation::titerator::titerator(Triangle * t) {
    v.insert(t);
    q.push_front(t);
}

void Triangulation::titerator::operator++() {
    Triangle *s = q.front();
    q.pop_front();
    for (int i = 0; i < 3; i++) {
        if (s->n[i] != NULL) {
            if (v.find(s->n[i]) == v.end()) {
                q.push_back(s->n[i]);
                v.insert(s->n[i]);
            }
        }
    }
}

void Triangulation::titerator::operator++(int a) {
    ++(*this);
}

Triangle * Triangulation::titerator::operator*()
{
    if (q.size() > 0) {
        return q.front();
    } else {
        return NULL;
    }
}

Triangulation::titerator& Triangulation::titerator::operator=(const titerator & it) {
    q.clear();
    v.clear();
	
    foreach(Triangle * s, it.q) {
        q.push_back(s);
    }
	
    foreach (Triangle * s, it.v) {
        v.insert(s);
    }
	
    return *this;
}


bool Triangulation::titerator::operator==(const titerator & it) {
	
    if (q.size() == 0 && it.q.size() == 0) {
        return true;
    }
    if (q.size() == it.q.size() && v.size() == it.v.size() ) {
        return true;
    }
    return false;
	
}

bool Triangulation::titerator::operator!=(const titerator & it) {
    return !(*this == it);
}

Triangulation::titerator Triangulation::tbegin() {
    return titerator(p2t.begin()->second);
}

Triangulation::titerator Triangulation::tend() {
    titerator it;
    return it;
}


Triangulation::viterator::viterator() {
	
}

Triangulation::viterator::viterator(map< int, Triangle *>::iterator m) {
    mit = m;
}

void Triangulation::viterator::operator++() {
    mit++;
}

void Triangulation::viterator::operator++(int a) {
    ++(*this);
}

int Triangulation::viterator::operator*()
{
    return mit->first;
}

Triangulation::viterator& Triangulation::viterator::operator=(const viterator & it) {
    mit = it.mit;
	
    return *this;
}


bool Triangulation::viterator::operator==(const viterator & it) {
	
    if (it.mit == mit) {
        return true;
    }
    return false;
	
}

bool Triangulation::viterator::operator!=(const viterator & it) {
    return !(*this == it);
}

Triangulation::viterator Triangulation::vbegin() {
    return viterator(p2t.begin());
}

Triangulation::viterator Triangulation::vend() {
    return viterator(p2t.end());
}
