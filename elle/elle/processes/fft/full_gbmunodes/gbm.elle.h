#ifndef _jkbgbm_elle_h
#define _jkbgbm_elle_h
#include "attrib.h"
#include "check.h"
#include "crossings.h"
#include "display.h"
#include "error.h"
#include "file.h"
#include "general.h"
#include "init.h"
#include "interface.h"
#include "mineraldb.h"
#include "movenode.h"
#include "nodes.h"
#include "polygon.h"
#include "runopts.h"
#include "stats.h"
#include "update.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
// #include "growthstats.h"
/*#define PI 3.1415926

#define DTOR PI/180
#define RTOD 180/PI*/

#include "triattrib.h"
#include "unodes.h"

int GBMGrowth();

int InitGrowth();

double GetNodeEnergy(int node, Coords *xy);
int GGMoveNode_all(int node, Coords *xy);

// int GGMoveNode(int node,Coords *xy);
int write_data(int stage);
int GetBodyNodeEnergy(int n, double *total_energy);

double GetNodeStoredEnergy(int node, Coords *xy);
double area_swept_gpclip(int node, Coords *xy);
double density_unodes(int node, Coords *xy);
//

#endif
