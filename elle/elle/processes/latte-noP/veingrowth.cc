/******************************************************
 * Spring Code Mike 2.0
 *
 * Functions for lattice class in lattice.cc
 *
 * Basic Spring Code for fracturing
 *
 *
 *
 *
 * Daniel Koehn and Jochen Arnold feb. 2002 to feb. 2003
 * Oslo, Mainz
 *
 * Daniel Koehn dec. 2003
 *
 * We thank Anders Malthe-S�renssen for his enormous help
 * and introduction to these codes
 *
 * Daniel Koehn and Till Sachau 2004/2005
 ******************************************************/

// ------------------------------------
// system headers
// ------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

// ------------------------------------
// elle headers
// ------------------------------------

#include "unodes.h" // include unode funct. plus undodesP.h
#include "veingrowth.h"
// for c++
#include "attrib.h" // enth�lt struct coord (hier: xy)
#include "attribarray.h"
#include "attribute.h"
#include "check.h"
#include "convert.h"
#include "crossings.h"
#include "display.h"
#include "error.h"
#include "file.h"
#include "general.h"
#include "interface.h"
#include "nodes.h"
#include "nodesP.h"
#include "polygon.h"
#include "runopts.h"
#include "tripoly.h"
#include "update.h"

// have to define these in the new Elle version

using std::cout;
using std::endl;
using std::vector;

// CONSTRUCTOR
/*******************************************************
 * Constructor for lattice class
 * calls the MakeLattice function and defines some variables
 * at the moment. The Constructor then builds a lattice once
 * the lattice class is called (i.e. by the user in the "elle"
 * function).
 *
 * Daniel spring 2002
 *
 * add security stop for pictdumps
 *
 * Daniel march 2003
 ********************************************************/

// ---------------------------------------------------------------
// Constructor of Lattice class
// ---------------------------------------------------------------

VeinGrowth::VeinGrowth()

// -----------------------------------------------------------------------------
// default values for some variables. Variables are declared in
// lattice header
// -----------------------------------------------------------------------------

{
  TotalTime = 0;
}

int VeinGrowth::DoVeinGrowth(int step) {
  int i, j, jj, max, hostrock, interface, count;
  int neigh[3], rgn[3]; // neighbours of the node j, 3 neighbours, NO_NB if
                        // empty
  Coords moved, new_position, old_position;
  double distance;
  // int interval=0,st_interval=0,err=0,max;
  // int *seq;
  // char fname[32];
  // FILE *fp;
  float ran_nb;

  max = ElleMaxNodes();
  // if (step%2) {

  int list[200000];
  for (j = 0; j < 200000; j++)
    list[j] = 0;

  for (jj = 0; jj < max; jj++) {
    do {
      ran_nb = rand() / (float)RAND_MAX;
      ran_nb *= max;

      j = int(ran_nb);
      // cout << j << endl;
      // cout << list[j] << endl;
    } while (list[j] == 1);
    list[j] = 1;

    if (ElleNodeIsActive(j)) {
      interface = 0;
      count = 0;
      ElleNeighbourNodes(j, neigh);
      for (i = 0; i < 3; i++) {
        if (neigh[i] != NO_NB) {
          ElleNeighbourRegion(j, neigh[i], &rgn[i]);
          ElleGetFlynnIntAttribute(rgn[i], &hostrock, EXPAND);
          interface += hostrock;
          count += 1; // counts number of neighbours,
        }
      }
      if (count == 2) // if node has two neighbours
      {
        if (interface == 1) // if neighbourRegions are different meaning fluid
                            // hostrock interface
        {
          GrowDoubleJ(j); // move double node into fluid
                          // MoveDoubleJ(j); //surface Energy graingrowth
        }
        // ElleCheckDoubleJ(j);
      } else // if count is other than two, there is a triple node
      {
        // if (interface == 1 || interface == 2) //checks wether tripleNode sits
        // at a interface
        {
          GrowTripleJ(j); // move triple node into fluid
          // MoveTripleJ(j);  //surface Energy graingrowth
        }
        // ElleCheckTripleJ(j);
      }
    }
  }
  // for (j=0;j<max;j++) {
  // if (ElleNodeIsActive(j))
  //{
  // interface = 0;
  // count = 0;
  // ElleNeighbourNodes(j,neigh);
  // for (i=0;i<3;i++)
  //{
  // if (neigh[i]!=NO_NB)
  //{
  // ElleNeighbourRegion(j,neigh[i],&rgn[i]);
  // ElleGetFlynnIntAttribute(rgn[i],&hostrock,EXPAND);
  // interface +=hostrock;
  // count +=1;  // counts number of neighbours,
  //}
  //}
  // if (count == 2)		// if node has two neighbours
  //{
  // if (interface == 1)		//if neighbourRegions are different
  // meaning fluid hostrock interface
  //{
  // ElleNodePrevPosition(j,&new_position);
  // ElleNodePosition(j,&old_position);
  // ElleSetPrevPosition(j,&old_position);
  // ElleSetPosition(j,&new_position);

  // ElleRelPosition(&old_position,j,&moved,&distance);

  // ElleCrossingsCheck(j, &moved);
  //}
  ////ElleCheckDoubleJ(j);
  //}
  // else // if count is other than two, there is a triple node
  //{
  ////if (interface == 1 || interface == 2) //checks wether tripleNode sits at a
  ///interface
  //{
  // ElleNodePrevPosition(j,&new_position);
  // ElleNodePosition(j,&old_position);
  // ElleSetPrevPosition(j,&old_position);
  // ElleSetPosition(j,&new_position);

  // ElleRelPosition(&old_position,j,&moved,&distance);

  // ElleCrossingsCheck(j, &moved);

  //}
  ////ElleCheckTripleJ(j);
  //}

  //}
  //}
  // for (j=0;j<max;j++) {
  // if (ElleNodeIsActive(j))
  //{
  // ElleNodePrevPosition(j,&old_position);
  // ElleRelPosition(&old_position,j,&moved,&distance);
  ////ElleCrossingsCheck(j, &moved);
  //}
  //}
  for (j = 0; j < max; j++) {
    if (ElleNodeIsActive(j)) {
      if (ElleNodeIsDouble(j))
        ElleCheckDoubleJ(j);
      else
        ElleCheckTripleJ(j);
    }
  }

  //// }
  ////else {
  ////for (j=max-1;j>=0;j--) {
  ////if (ElleNodeIsActive(j)) {
  ////if (ElleNodeIsDouble(j)) {
  ////MoveDoubleJ(j);
  ////ElleCheckDoubleJ(j);
  ////}
  ////else if (ElleNodeIsTriple(j)) {
  ////MoveTripleJ(j);
  ////ElleCheckTripleJ(j);
  ////}
  ////}
  ////}
  ////}

  ElleUpdate();
}

void VeinGrowth::GrowDoubleJ(int node) {
  int i, j, max, hostrock;
  int neigh[3], rgn[3];
  int count;
  double dx, dy, distance;
  Coords xyj, xyneigh, move_node, prevpos;
  Coords v1nu, v2nu; // normal unit vectors on vectors from nodes to neighbours
  Coords jnew;       // new node coordinates

  j = node;
  ElleNeighbourNodes(j, neigh);
  ElleNodePosition(j, &xyj);

  count = 0;
  for (i = 0; i < 3; i++) {
    if (neigh[i] != NO_NB) {
      ElleNeighbourRegion(j, neigh[i], &rgn[i]);
      ElleGetFlynnIntAttribute(rgn[i], &hostrock, EXPAND);
      if (hostrock == 1) // fluid is to the right of the vector
      {

        ElleNodePosition(neigh[i], &xyneigh);
        // cout << "xyj.x is: "<<xyj.x << endl;
        // cout << "xyj.y is: "<<xyj.y << endl;
        // cout << "xyneigh.x is: "<<xyneigh.x << endl;
        // cout << "xyneigh.y is: "<<xyneigh.y << endl;
        dx = xyj.x - xyneigh.x;
        dy = xyj.y - xyneigh.y;
        distance = sqrt((dx * dx) + (dy * dy));
        // cout << "dx is: " << dx << endl;
        // cout << "dy is: " << dy << endl;
        // cout << "distance is: "<< distance << endl;
        if (!count) {
          v1nu.x = dy / distance * (1);  // calculates normal unit vector
          v1nu.y = dx / distance * (-1); // calculates normal unit vector
        } else {
          v2nu.x = dy / distance * (1);  // calculates normal unit vector
          v2nu.y = dx / distance * (-1); // calculates normal unit vector
        }

      } else // fluid is to the left of the vector
      {

        ElleNodePosition(neigh[i], &xyneigh);
        dx = xyj.x - xyneigh.x;
        dy = xyj.y - xyneigh.y;
        distance = sqrt((dx * dx) + (dy * dy));
        if (!count) {
          v1nu.x = dy / distance * (-1);
          v1nu.y = dx / distance * (1);
        } else {
          v2nu.x = dy / distance * (-1); // calculates normal unit vector
          v2nu.y = dx / distance * (1);  // calculates normal unit vector
        }
      }
      count++;
    }
  }
  v1nu.x *= -0.0005;
  v1nu.y *= -0.0005;
  v2nu.x *= -0.0005;
  v2nu.y *= -0.0005;
  jnew.x = xyj.x + (v1nu.x + v2nu.x) / 2;
  jnew.y = xyj.y + (v1nu.y + v2nu.y) / 2;
  move_node.x = (v1nu.x + v2nu.x) / 2;
  move_node.y = (v1nu.y + v2nu.y) / 2;
  // cout << jnew.x << endl;
  // ElleNodePrevPosition(j,&prevpos);
  ElleSetPosition(j, &jnew);
  // ElleSetPrevPosition (j, &jnew);
  ElleCrossingsCheck(j, &move_node);
}

void VeinGrowth::GrowTripleJ(int node) {
  int i, j, max, hostrock, subinterface;
  int neigh[3], rgn[3];
  int count;
  double dx, dy, distance;
  Coords xyj, xyneigh, move_node, prevpos;
  Coords v1nu, v2nu, v1u,
      v2u;     // normal unit vectors on vectors from nodes to neighbours
  Coords jnew; // new node coordinates
  Coords tripleJhelpNeighbour1,
      tripleJhelpNeighbour2; // help neighbour points to calculate intersection
                             // of adjacent boundaries
  Coords helpNeighbourMidPoint; // Midpoint between help neighbour points
  double b;                     // distance between midpoint and help neighbours

  v1nu.x = 0.0;
  v1nu.y = 0.0;
  v2nu.x = 0.0;
  v2nu.y = 0.0;

  j = node;
  ElleNeighbourNodes(j, neigh);
  ElleNodePosition(j, &xyj);

  count = 0;

  for (i = 0; i < 3; i++) {
    if (neigh[i] != NO_NB) {
      subinterface = 0;

      ElleNeighbourRegion(neigh[i], j, &rgn[i]);
      ElleGetFlynnIntAttribute(rgn[i], &hostrock, EXPAND);
      subinterface += hostrock;

      ElleNeighbourRegion(j, neigh[i], &rgn[i]);
      ElleGetFlynnIntAttribute(rgn[i], &hostrock, EXPAND);
      subinterface += hostrock;

      if (subinterface == 1) // check if neighbour node is relevant
      {
        if (hostrock == 1) {

          ElleNodePosition(neigh[i], &xyneigh);
          // cout << "xyj.x is: "<<xyj.x << endl;
          // cout << "xyj.y is: "<<xyj.y << endl;
          // cout << "xyneigh.x is: "<<xyneigh.x << endl;
          // cout << "xyneigh.y is: "<<xyneigh.y << endl;
          dx = xyj.x - xyneigh.x;
          dy = xyj.y - xyneigh.y;
          distance = sqrt((dx * dx) + (dy * dy));
          // cout << "dx is: " << dx << endl;
          // cout << "dy is: " << dy << endl;
          // cout << "distance is: "<< distance << endl;
          if (!count) {
            v1nu.x = dy / distance * (1);  // calculates normal unit vector
            v1nu.y = dx / distance * (-1); // calculates normal unit vector
          } else {
            v2nu.x = dy / distance * (1);  // calculates normal unit vector
            v2nu.y = dx / distance * (-1); // calculates normal unit vector
          }
        }

        else // Flynn to the left is hostrock
        {

          ElleNodePosition(neigh[i], &xyneigh);
          dx = xyj.x - xyneigh.x;
          dy = xyj.y - xyneigh.y;
          distance = sqrt((dx * dx) + (dy * dy));
          if (!count) {
            v1nu.x = dy / distance * (-1);
            v1nu.y = dx / distance * (1);
            v1u.x = dx / distance;
            v1u.y = dy / distance;
          } else {
            v2nu.x = dy / distance * (-1); // calculates normal unit vector
            v2nu.y = dx / distance * (1);  // calculates normal unit vector
          }
        }
        count++;
      }
    }
  }

  tripleJhelpNeighbour1.x = xyj.x + v1u.x;
  tripleJhelpNeighbour1.y = xyj.y + v1u.y;

  tripleJhelpNeighbour2.x = xyj.x + v2u.x;
  tripleJhelpNeighbour2.y = xyj.y + v2u.y;

  helpNeighbourMidPoint.x =
      (tripleJhelpNeighbour1.x + tripleJhelpNeighbour2.x) / 2;
  helpNeighbourMidPoint.y =
      (tripleJhelpNeighbour1.y + tripleJhelpNeighbour2.y) / 2;

  dx = tripleJhelpNeighbour1.x - tripleJhelpNeighbour2.x;
  dy = tripleJhelpNeighbour1.y - tripleJhelpNeighbour2.y;
  b = sqrt((dx * dx) + (dy * dy));

  if (b != 1) {
  }

  v1nu.x *= -0.0005;
  v1nu.y *= -0.0005;
  v2nu.x *= -0.0005;
  v2nu.y *= -0.0005;
  jnew.x = xyj.x + (v1nu.x + v2nu.x) / 2;
  jnew.y = xyj.y + (v1nu.y + v2nu.y) / 2;
  // cout << "xyj.x: " << xyj.x << endl;
  // cout << "xyj.y: " << xyj.y << endl;
  // cout << "triple was moved" << endl;
  // cout << "jnew.x: " << jnew.x << endl;
  // cout << "jnew.y: " << jnew.y << endl;
  // cout << jnew.x << endl;
  // ElleNodePrevPosition(j,&prevpos);
  ElleSetPosition(j, &jnew);
  // ElleSetPrevPosition (j, &jnew);
  move_node.x = (v1nu.x + v2nu.x) / 2;
  move_node.y = (v1nu.y + v2nu.y) / 2;
  ElleCrossingsCheck(j, &move_node);
}

int VeinGrowth::MoveDoubleJ(int node1) {
  int i, nghbr[2], nbnodes[3], err;
  double maxV, gb_energy, ray, deltaT, vlen;
  double switchDist, speedUp;
  Coords xy1, movedist;

  switchDist = ElleSwitchdistance();
  speedUp = ElleSpeedup() * switchDist * switchDist * 0.02;
  maxV = ElleSwitchdistance() / 5.0;
  /*
   * allows speedUp to be 1 in input file
   */
  gb_energy = speedUp;
  deltaT = 0.0;
  /*
   * find the node numbers of the neighbours
   */
  if (err = ElleNeighbourNodes(node1, nbnodes))
    OnError("MoveDoubleJ", err);
  i = 0;
  while (i < 3 && nbnodes[i] == NO_NB)
    i++;
  nghbr[0] = nbnodes[i];
  i++;
  while (i < 3 && nbnodes[i] == NO_NB)
    i++;
  nghbr[1] = nbnodes[i];

  GetRay(node1, nghbr[0], nghbr[1], &ray, &movedist);
  if (ray > 0.0) {
    /*if (ray > ElleSwitchdistance()/100.0) {*/
    vlen = gb_energy / ray;
    if (vlen > maxV) {
      vlen = maxV;
      deltaT = 1.0;
    }
    if (vlen > 0.0) {
      movedist.x *= vlen;
      movedist.y *= vlen;
    } else {
      movedist.x = 0.0;
      movedist.y = 0.0;
    }
    TotalTime += deltaT;
    ElleUpdatePosition(node1, &movedist);
  } else {
    vlen = 0.0;
  }
}

int VeinGrowth::MoveTripleJ(int node1) {
  int i, nghbr[3], finished = 0, err = 0;
  double maxV, gb_energy[3], ray[3], deltaT, vlen[3], vlenTriple;
  double switchDist, speedUp;
  Coords xy1, movedist[3], movedistTriple;

  switchDist = ElleSwitchdistance();
  /*
   * allows speedUp to be 1 in input file
   */
  speedUp = ElleSpeedup() * switchDist * switchDist * 0.02;
  maxV = switchDist / 5.0;
  for (i = 0; i < 3; i++)
    gb_energy[i] = speedUp;
  deltaT = 0.0;
  /*
   * find the node numbers of the neighbours
   */
  if (err = ElleNeighbourNodes(node1, nghbr))
    OnError("MoveTripleJ", err);

  GetRay(node1, nghbr[0], nghbr[1], &ray[0], &movedist[0]);
  GetRay(node1, nghbr[1], nghbr[2], &ray[1], &movedist[1]);
  GetRay(node1, nghbr[2], nghbr[0], &ray[2], &movedist[2]);
  for (i = 0; i < 3; i++) {
    if (ray[i] > 0.0) {
      /*if (ray[i] > switchDist/100.0) {*/
      vlen[i] = gb_energy[i] / ray[i];
      /*if (vlen[i] > maxV) vlen[i] = maxV;*/
    } else {
      vlen[i] = 0.0;
      finished = 1;
    }
  }
  if (!finished) {
    for (i = 0; i < 3; i++) {
      if (vlen[i] < maxV) {
        movedist[i].x *= vlen[i];
        movedist[i].y *= vlen[i];
      } else {
        movedist[i].x *= maxV;
        movedist[i].y *= maxV;
      }
    }
    movedistTriple.x = movedist[0].x + movedist[1].x + movedist[2].x;
    movedistTriple.y = movedist[0].y + movedist[1].y + movedist[2].y;
    vlenTriple = sqrt(movedistTriple.x * movedistTriple.x +
                      movedistTriple.y * movedistTriple.y);
    if (vlenTriple > maxV) {
      vlenTriple = maxV / vlenTriple;
      movedistTriple.x *= vlenTriple;
      movedistTriple.y *= vlenTriple;
      deltaT = 1.0;
    }
    if (vlenTriple <= 0.0)
      movedistTriple.x = movedistTriple.y = 0.0;

    TotalTime += deltaT;
  } else {
    ElleNodePosition(node1, &xy1);
    ElleNodePlotXY(nghbr[0], &movedist[0], &xy1);
    ElleNodePlotXY(nghbr[1], &movedist[1], &xy1);
    ElleNodePlotXY(nghbr[2], &movedist[2], &xy1);
    for (i = 0; i < 3; i++) {
      movedist[i].x = movedist[i].x - xy1.x;
      movedist[i].y = movedist[i].y - xy1.y;
    }
    movedistTriple.x = (movedist[0].x + movedist[1].x + movedist[2].x) / 2.0;
    movedistTriple.y = (movedist[0].y + movedist[1].y + movedist[2].y) / 2.0;
#if XY
    vlenTriple = sqrt(movedistTriple.x * movedistTriple.x +
                      movedistTriple.y * movedistTriple.y);
    if (vlenTriple > maxV) {
      vlenTriple = maxV / vlenTriple;
      movedistTriple.x *= vlenTriple;
      movedistTriple.y *= vlenTriple;
      deltaT = 1.0;
    }
#endif
  }
  ElleUpdatePosition(node1, &movedistTriple);
}

void VeinGrowth::GetRay(int node1, int node2, int node3, double *ray,
                        Coords *movedist) {
  double dx2, dy2, dx3, dy3, tmpx, tmpy;
  double k, x0, y0;
  double switchDist;
  double eps = 1e-8;
  double r;
  Coords xy1, xy2, xy3;

  switchDist = ElleSwitchdistance();
  ElleNodePosition(node1, &xy1);
  ElleNodePlotXY(node2, &xy2, &xy1);
  ElleNodePlotXY(node3, &xy3, &xy1);
  dx2 = xy2.x - xy1.x;
  dy2 = xy2.y - xy1.y;
  dx3 = xy3.x - xy1.x;
  dy3 = xy3.y - xy1.y;
  if (dx2 == 0.0) {
    tmpx = dx2;
    tmpy = dy2;
    dx2 = dx3;
    dy2 = dy3;
    dx3 = tmpx;
    dy3 = tmpy;
  }
  *ray = 0.0;
  movedist->x = movedist->y = 0.0;

  if (dx2 > eps || dx2 < -eps) {
    k = 2.0 * dx3 * dy2 / dx2 - 2.0 * dy3;
    if (k != 0.0) {
      y0 = ((dx3 / dx2) * (dx2 * dx2 + dy2 * dy2) - dx3 * dx3 - dy3 * dy3) / k;
      x0 = (dx2 * dx2 + dy2 * dy2 - 2.0 * y0 * dy2) / (2.0 * dx2);
      r = sqrt((double)(x0 * x0 + y0 * y0));
      if (r != 0.0) {
        *ray = r;
        movedist->x = x0 / *ray;
        movedist->y = y0 / *ray;
        if (*ray < switchDist / 3.0)
          *ray = switchDist / 3.0;
      }
    }
  }
}

int VeinGrowth::IncreaseAngle(Coords *xy, Coords *xy1, Coords *xy2,
                              Coords *diff) {
  Coords xynew;
  /*
   * find the average of the 3 points
   * move 0.1 of the distance towards the average point
   */
  xynew.x = (xy->x + xy1->x + xy2->x) / 3;
  xynew.y = (xy->y + xy1->y + xy2->y) / 3;
  diff->x = xynew.x - xy->x;
  diff->y = xynew.y - xy->y;
  diff->x *= 0.1;
  diff->y *= 0.1;
}

/* #define MINANG 0.087  approx 5deg */
#define MINANG 0.14 /* approx 8deg */
void VeinGrowth::CheckAngles() {
  int moved = 1, i, j, k, max;
  int nbnodes[3];
  double currang;
  double ang;
  Coords xy[3], movedist;

  max = ElleMaxNodes();
  while (moved) {
    for (k = 0, moved = 0; k < max; k++) {
      if (ElleNodeIsActive(k)) {
        ElleNodePosition(k, &xy[0]);
        ElleNeighbourNodes(k, nbnodes);
        if (ElleNodeIsDouble(k)) {
          j = 0;
          i = 1;
          while (j < 3) {
            if (nbnodes[j] != NO_NB)
              ElleNodePlotXY(nbnodes[j], &xy[i++], &xy[0]);
            j++;
          }
          if (angle(xy[0].x, xy[0].y, xy[1].x, xy[1].y, xy[2].x, xy[2].y,
                    &currang))
            OnError("angle error", 0);
          ang = fabs(currang);
          if (ang < MINANG /*|| ang<(-M_PI+MINANG)*/) {
            do {
              IncreaseAngle(&xy[0], &xy[1], &xy[2], &movedist);
              xy[0].x += movedist.x;
              xy[0].y += movedist.y;
              if (angle(xy[0].x, xy[0].y, xy[1].x, xy[1].y, xy[2].x, xy[2].y,
                        &currang))
                OnError("angle error", 0);
              ang = fabs(currang);
            } while (ang < MINANG /*|| ang>(M_PI-MINANG)*/);
            ElleSetPosition(k, &xy[0]);
            ElleCheckDoubleJ(k);
            moved = 1;
          }
        } else if (ElleNodeIsTriple(k)) {
          for (j = 0; j < 3; j++) {
            i = (j + 1) % 3;
            ElleNodePlotXY(nbnodes[j], &xy[1], &xy[0]);
            ElleNodePlotXY(nbnodes[i], &xy[2], &xy[0]);
            if (angle(xy[0].x, xy[0].y, xy[1].x, xy[1].y, xy[2].x, xy[2].y,
                      &currang))
              OnError("angle error", 0);
            ang = fabs(currang);
            if (ang < MINANG /*|| ang>(M_PI-MINANG)*/) {
              do {
                IncreaseAngle(&xy[0], &xy[1], &xy[2], &movedist);
                xy[0].x += movedist.x;
                xy[0].y += movedist.y;
                if (angle(xy[0].x, xy[0].y, xy[1].x, xy[1].y, xy[2].x, xy[2].y,
                          &currang))
                  OnError("angle error", 0);
                ang = fabs(currang);
              } while (ang < MINANG /* || ang>(M_PI-MINANG)*/);
              ElleSetPosition(k, &xy[0]);
              ElleCheckTripleJ(k);
              moved = 1;
              j = 3;
            }
          }
        }
      }
    }
  }
}
