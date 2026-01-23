/*****************************************************
 * Copyright: (c) L. A. Evans
 * File:      $RCSfile: polyutils.h,v $
 * Revision:  $Revision: 1.2 $
 * Date:      $Date: 2006/05/18 06:30:15 $
 * Author:    $Author: levans $
 *
 ******************************************************/
#ifndef _Elle_polyutils_h
#define _Elle_polyutils_h

#include "attrib.h"
#include "triangle.h"
#include <list>
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif
#ifdef __cplusplus
}
#endif
void sortCoordsOnAngle(std::vector<Coords> &pts, Coords *ref);
void findPolygon(REAL *pointlist, REAL *normlist, int *edgelist, int nseg,
                 int *seg_count, int *segs_left, std::list<int> &processed);
#endif
