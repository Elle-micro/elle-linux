#ifndef _gbm_movenode_h
#define _gbm_movenode_h
#include "attrib.h"
#include "check.h"
#include "convert.h"
#include "display.h"
#include "error.h"
#include "file.h"
#include "general.h"
#include "init.h"
#include "interface.h"
#include "log.h"
#include "mineraldb.h"
#include "nodes.h"
#include "runopts.h"
#include "stats.h"
#include "update.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#define TIME_REDUCE 2

int MoveNode(int node1, Coords pvec, Coords *m);
double GetVectorLength(Coords vec);
double DEGCos(Coords vec1, Coords vec2);
int GetMoveDir(int node, double e1, double e2, double e3, double e4,
               Coords *newpos, double t);
int GetNewPos(int node, double e1, double e2, double e3, double e4,
              Coords *newpos, double t);
double GetBoundaryMobility(int node, int nb);

#endif
