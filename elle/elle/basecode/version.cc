/*****************************************************
 * Copyright: (c) L. A. Evans
 * File:      $RCSfile: version.cc,v $
 * Revision:  $Revision: 1.3 $
 * Date:      $Date: 2004/03/18 03:03:35 $
 * Author:    $Author: levans $
 *
 ******************************************************/
#include "version.h"
#include "attrib.h"
#include "runopts.h"
#include "timefn.h"
#include "versiondef.h"

/*****************************************************

static const char rcsid[] =
       "$Id: version.cc,v 1.3 2004/03/18 03:03:35 levans Exp $";

******************************************************/
std::string ElleGetLibVersionString() {
  std::string vers = Version_num;
  vers += Patch_level;
  return (vers);
}

std::string ElleGetLocalTimeString() {
  std::string x;
  x += GetLocalTime();
  return (x);
}

std::string ElleGetCreationString() {
  std::string x = "# Created by ";
  x += ElleAppName();
  x += ": elle version ";
  x += ElleGetLibVersionString();
  x += "  ";
  x += GetLocalTime();
  return (x);
}

const char *ElleGetCreationCString(void) {
  static std::string cached; // static to keep storage alive after return
  cached.clear();
  cached += "# Created by ";
  cached += ElleAppName();
  cached += ": elle version ";
  cached += ElleGetLibVersionString();
  cached += "  ";
  cached += GetLocalTime();
  return cached.c_str();
}
