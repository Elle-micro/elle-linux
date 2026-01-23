#ifndef _movenode_h
#define _movenode_h
#include "attrib.h"
#include "check.h"
#include "convert.h"
#include "display.h"
#include "error.h"
#include "file.h"
#include "general.h"
#include "init.h"
#include "interface.h"
#include "mineraldb.h"
#include "nodes.h"
#include "polygon.h"
#include "runopts.h"
#include "stats.h"
#include "update.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <vector>
int MoveTNode(int flynn, int node1, Coords pvec, Coords *m);
int MoveDNode(int flynn, int node1, Coords pvec, Coords *m);
int GetPVector(int node, Coords *pvec, double E1, double E2, double E3,
               double E4, double distance);
double GetVectorLength(Coords vec);
double DEGCos(Coords vec1, Coords vec2);
double ChangeLength(Coords p1, Coords p2);
int GetMoveDir(int flynn, int node, double e1, double e2, double e3, double e4,
               Coords *newpos, double distance);
int GetMoveDir(int flynn, int node, Coords xyz[3], double E[3],
               Coords *movedir);
int GetNewPos(int flynn, int node, double e1, double e2, double e3, double e4,
              Coords *newpos, double distance);
int GetNewPos(int flynn, int node, Coords xyz[3], double E[3], Coords *newpos);
double GetBoundaryMobility(int flynn);

int GetPVector(int node, Coords *pvec, Coords xyz[3], double E[3]);
#define LENGTHSCALE 1e-3

/* !Only needed for debugging. Would draw the gradients and pvector on the
 * Canvas. Looks nice, but takes very long! */
extern void DrawGrad(double list[][3], Coords p1, Coords gvector);
#endif
