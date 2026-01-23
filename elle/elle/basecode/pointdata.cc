/*****************************************************
 * Copyright: (c) L. A. Evans
 * File:      $RCSfile: pointdata.cc,v $
 * Revision:  $Revision: 1.2 $
 * Date:      $Date: 2004/03/18 03:03:35 $
 * Author:    $Author: levans $
 *
 ******************************************************/
#include "pointdata.h"
#include <iostream>

/*****************************************************

static const char rcsid[] =
       "$Id: pointdata.cc,v 1.2 2004/03/18 03:03:35 levans Exp $";

******************************************************/
// Write object to an ostream reference
std::ostream &operator<<(std::ostream &os, const PointData &t) {
  std::cout << t.x << ", " << t.y << ", " << t.id << ", " << t.bnd << '\n';
  return os;
}

unsigned char PointData::gt(PointData &t) {
  unsigned char gt = 0;

  if (bnd > t.boundary())
    return (1);
  if (bnd == t.boundary())
    switch (bnd) {
    case 1:
      if (x > t.x)
        gt = 1;
      break;
    case 2:
      if (y > t.y)
        gt = 1;
      break;
    case 3:
      if (x < t.x)
        gt = 1;
      break;
    case 4:
      if (y < t.y)
        gt = 1;
      break;
    }
  return (gt);
}
