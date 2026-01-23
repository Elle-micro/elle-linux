/*****************************************************
 * Copyright: (c) L. A. Evans
 * File:      $RCSfile: error.h,v $
 * Revision:  $Revision: 1.2 $
 * Date:      $Date: 2012/08/16 06:42:34 $
 * Author:    $Author: levans $
 *
 ******************************************************/
#ifndef _E_error_h
#define _E_error_h

#ifndef _E_errnum_h
#include "errnum.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif
void OnError(const char *message, int err_num);
void CleanUp(void);
#ifdef __cplusplus
}
#endif
#endif
