/*****************************************************
 * Copyright: (c) L. A. Evans
 * File:      $RCSfile: batch.h,v $
 * Revision:  $Revision: 1.1.1.1 $
 * Date:      $Date: 2002/09/18 03:17:14 $
 * Author:    $Author: levans $
 *
 ******************************************************/
#ifndef _E_display_h
#define _E_display_h
/* matches def in menus.c */
#define GRAINS 12
#define SUBGRAINS 13
#define UNITS 14
#define TRIANGLES 15

#ifdef __cplusplus
extern "C" {
#endif
void ElleUpdateDisplay();
void ElleUpdateSettings();

#ifdef __cplusplus
}
#endif
#endif
