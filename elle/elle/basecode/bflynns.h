/*****************************************************
 * Copyright: (c) L. A. Evans
 * File:      $RCSfile: bflynns.h,v $
 * Revision:  $Revision: 1.2 $
 * Date:      $Date: 2005/07/12 07:16:25 $
 * Author:    $Author: levans $
 *
 ******************************************************/
#ifndef _E_bflynns_h
#define _E_bflynns_h

#include <zlib.h>
#ifdef __cplusplus
#include <cstdio> // C++: ensure FILE is declared
#else
#include <stdio.h> // C: ensure FILE is declared
#endif

#ifndef _E_attrib_h
#include "attrib.h"
#endif

#define GRN_FILL_COL 20

int LoadZIPParentData(gzFile in, char str[]);
int SaveZIPParentData(gzFile in);
int LoadZIPFlynnData(gzFile in, char str[]);
int SaveZIPFlynnData(gzFile in);

#ifdef __cplusplus
extern "C" {
#endif
int ElleNewGrain(int first, int last, ERegion old, ERegion *rnew);
void ElleNucleateGrainFromNode(int node, ERegion *new_rgn);
void ElleCleanArrays(void);
void ElleRemoveArrays(void);
int ElleCheckForTwoSidedGrain(int node, int *nb);
int ElleReadParentData(FILE *fp, char str[]);
int ElleWriteParentData(FILE *fp);
int ElleReadFlynnData(FILE *fp, char str[]);
int ElleWriteFlynnData(FILE *fp);
#ifdef __cplusplus
}
#endif
#endif
