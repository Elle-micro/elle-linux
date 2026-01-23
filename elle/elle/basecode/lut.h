/*****************************************************
 * Copyright: (c) L. A. Evans
 * File:      $RCSfile: lut.h,v $
 * Revision:  $Revision: 1.3 $
 * Date:      $Date: 2012/04/20 03:30:10 $
 * Author:    $Author: levans $
 *
 ******************************************************/
#ifndef _E_lut_h
#define _E_lut_h
#include "gz_utils.h"
#include <zlib.h>

const int LUT_MAX = 1024;

typedef struct {
  int size;
  double *data;
} EnergyLUT;

int LoadZIPGBEnergyLUT(gzFile in, char str[]);
int SaveZIPGBEnergyLUT(gzFile in);
#ifdef __cplusplus
extern "C" {
#endif
int ElleReadGBEnergyLUT(FILE *fp, char str[]);
int ElleWriteGBEnergyLUT(FILE *fp);
void ElleInitEnergyLUT(int size);
int ElleEnergyLUTSize();
void ElleSetEnergyLUT(int index, double val);
double ElleEnergyLUTValue(int index);
void ElleRemoveEnergyLUT();

#ifdef __cplusplus
}
#endif
#endif
