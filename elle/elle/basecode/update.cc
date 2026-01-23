/*****************************************************
 * Copyright: (c) L. A. Evans
 * File:      $RCSfile: update.cc,v $
 * Revision:  $Revision: 1.3 $
 * Date:      $Date: 2007/06/28 13:15:36 $
 * Author:    $Author: levans $
 *
 ******************************************************/
#include "update.h"
#include "display.h"
#include "file.h"
#include "runopts.h"
/*****************************************************

static const char rcsid[] =
       "$Id: update.cc,v 1.3 2007/06/28 13:15:36 levans Exp $";

******************************************************/

int ElleUpdate() {
  ElleIncrementCount();
  if (ElleDisplay())
    ElleUpdateDisplay();
  /*ElleCheckFiles();*/
  return 0; // ensure defined return value
}
