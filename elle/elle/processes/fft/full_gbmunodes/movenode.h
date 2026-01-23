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
#include "polygon.h"
#include "runopts.h"
#include "stats.h"
#include "triattrib.h"
#include "unodes.h"
#include "update.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <vector>

// Ice mineral attributes
#define GB_MOBILITY_ice 1
#define GB_MOBILITY_2_ice 2
#define DD_ENERGY_ice 3
#define SURFACE_ENERGY_ice 4
#define GB_ACTIVATION_ENERGY_ice 5

int MoveTNode(int node1, Coords pvec, Coords *m);
int MoveDNode(int node1, Coords pvec, Coords *m);
double GetVectorLength(Coords vec);
double DEGCos(Coords vec1, Coords vec2);
// int GetMoveDir( int node, double e1, double e2, double e3, double e4, Coords
// * newpos,double t );
int GetMoveDir(int node, Coords pvec, Coords *newpos, double t);
int GetNewPos(int node, double e1, double e2, double e3, double e4,
              Coords *newpos, double t);
double GetBoundaryMobility(int node, int nb);
// int GetBodyNodeEnergy2(int n, double *total_energy);
double GetMineralAttribute_ice(int attribute);
//
int SearchUnode(int rgn, Coords *xy);
int ElleGetFlynnEulerCAxis(double alpha, double beta, double gamma,
                           Coords_3D *dircos);
#endif
