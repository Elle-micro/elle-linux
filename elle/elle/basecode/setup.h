/*----------------------------------------------------------------
 *    Elle:   setup.h  1.0  11 September 1997
 *
 *    Copyright (c) 1997 by L.A. Evans & T.D. Barr
 *----------------------------------------------------------------*/
#ifndef _E_setup_h
#define _E_setup_h

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif

int StartApp(void); // explicit void
int Run_App(FILE *);
int SetupApp(int, char **);
void Init_Data(void); // explicit void
#ifdef __cplusplus
}
#endif
#endif
