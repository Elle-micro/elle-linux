#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <math.h>
#include "attrib.h"
#include "check.h"
#include "convert.h"
#include "crossings.h"
#include "display.h"
#include "errnum.h"
#include "error.h"
#include "file.h"
#include "general.h"
#include "init.h"
#include "interface.h"
#include "lut.h"
#include "mineraldb.h"
#include "nodes.h"
#include "parseopts.h"
#include "polygon.h"
#include "runopts.h"
#include "setup.h"
#include "stats.h"
#include "timefn.h"
#include "triattrib.h"
#include "unodes.h"
#include "unodesP.h"
#include "update.h"

#ifndef _E_simple_h
#define _E_simple_h

#define  SIMPLE_FLYNN  0
#define  SIMPLE_BNODE  1
#define  SIMPLE_UNODE  2

int InitElle(int argc, char **argv, char *rootname, int ProcessType);
int InitProcessFlynn();
int InitProcessBnode();
int InitProcessUnode();

int ProcessFlynns();
int ProcessBnodes();
int SetUpUnodesByFlynn();

int MyFlynnProcess(int k);
int MyBnodeProcess(int k);
int MyUnodeProcess(int k, int flynnid);
int MyFlynnProcess();
int MyBnodeProcess();
int MyUnodeProcess();

void ProcessUnodesByFlynn(int flynnid);

#endif
