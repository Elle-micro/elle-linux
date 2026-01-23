/*****************************************************
 * Copyright: (c) L. A. Evans
 * File:      $RCSfile: version.h,v $
 * Revision:  $Revision: 1.2 $
 * Date:      $Date: 2003/12/19 01:15:24 $
 * Author:    $Author: levans $
 *
 ******************************************************/
#ifndef _E_version_h
#define _E_version_h
#include <string>

std::string ElleGetLibVersionString();
std::string ElleGetLocalTimeString();
std::string ElleGetCreationString();
const char *ElleGetCreationCString(void);

#endif
