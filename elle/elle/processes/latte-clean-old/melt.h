#include "check.h"
#include "convert.h"
#include "crossings.h"
#include "display.h"
#include "errnum.h"
#include "error.h"
#include "file.h"
#include "general.h"
#include "init.h"
#include "interface.h"
#include "lut.h"
#include "melt.h"
#include "mineraldb.h"
#include "nodes.h"
#include "polygon.h"
#include "runopts.h"
#include "stats.h"
#include "timefn.h"
#include "update.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

class Melt {
public:
  /********************physical constants*********************************/
  double energyofdislocations; // in Jm-1    (from mineraldb)
  double energyofsurface;      // in Jm-2    (from mineraldb)
  double mobilityofboundary;   // in m2s-1J-1    (aka fudge factor)
  double dislocationdensityscaling =
      10e13; // in m-2    (ie 1.5 in elle file = 1.5e13 in real world)
  double truetimestep =
      3.1536e10; // in s     (1000yrs here, not including leap years)
  double lengthscale = 1e-3; // in m     (sides of box)
  double R = 8.314472;       // in Jmol-1K-1    (Gas constant)
  double Qgbm =
      200; // in Jmole-1    (Activation energy for GBM  (made up number))
  /***********************************************************************/

  // using std::cout;
  // using std::endl;

  // #define DEBUG

  // UserData CurrData;   ???

  const int MeltId = F_ATTRIB_A; // melt phase

  double GetAngleEnergy(double angle);
  int GBM_MoveNode(int n, Coords *movedir, int *same, double arealiq,
                   double factor, double fudge, double areaequil,
                   double energyxlxl, double energyliqxl, double energyliqliq);
  int GetBodyNodeEnergy(int n, double *total_energy);
  int GetCSLFactor(Coords_3D xyz[3], float *factor);
  int ElleGetFlynnEulerCAxis(int flynn_no, Coords_3D *dircos);
  int GetCSLFactor(Coords_3D xyz[3], float *factor);
  int ElleGetFlynnEulerCAxis(int flynn_no, Coords_3D *dircos);
  int GetSurfaceNodeEnergy(int n, int meltmineral, Coords *xy, double *energy,
                           double arealiq, double fudge, double areaequil,
                           double energyxlxl, double energyliqxl,
                           double energyliqliq);
  double CalcAreaLiq(int meltmineral);
  void ElleNodeSameAttrib(int node, int *same, int *samemin, int attrib_id);
  extern double ElleSwitchLength();

  double incroffset = 0.002; /* this is reassigned in GBE_GrainGrowth */
}
