/*****************************************************
 * Copyright: (c) L. A. Evans
 * File:      $RCSfile: file_utils.h,v $
 * Revision:  $Revision: 1.3 $
 * Date:      $Date: 2006/07/19 06:56:55 $
 * Author:    $Author: levans $
 *
 ******************************************************/
#ifndef _E_file_utils_h
#define _E_file_utils_h

#include "string_utils.h"
#include <stdio.h>

#ifdef __WIN32__
const char E_DIR_SEPARATOR = '\\';
#else
const char E_DIR_SEPARATOR = '/';
#endif

#ifdef __cplusplus
extern "C" {
#endif
int ElleSkipLine(FILE *fp, char *str);
int ElleSkipSection(FILE *fp, char *str, valid_terms keyterms[]);
int ElleCopyLine(FILE *fp, FILE *fpout, char *str);
int ElleCopySection(FILE *fp, FILE *fpout, char *str, valid_terms keyterms[]);
int FileExists(char *filename);
#ifdef __cplusplus
}
#endif
#endif
