
#ifndef _E_veingrowth_h
#define _E_veingrowth_h

#include "attrib.h"
#include "check.h"
#include "display.h"
#include "error.h"
#include "file.h"
#include "general.h"
#include "init.h"
#include "nodes.h"
#include "runopts.h"
#include "stats.h"
#include "update.h"
#include <cmath>
#include <cstdio>
#include <cstring>

class VeinGrowth {
public:
  double TotalTime;

  VeinGrowth();    // Constructor
  ~VeinGrowth(){}; // Destructor

  int DoVeinGrowth(int step);
  void GrowDoubleJ(int node);
  void GrowTripleJ(int node);
  int MoveDoubleJ(int node1);
  int MoveTripleJ(int node1);
  void GetRay(int node1, int node2, int node3, double *ray, Coords *movedist);
  void CheckAngles();
  int IncreaseAngle(Coords *xy, Coords *xy1, Coords *xy2, Coords *diff);
};

#endif
