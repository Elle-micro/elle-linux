
/*
 *  main.c
 */

#include "manuel.h"

#include <stdio.h>
#include <stdlib.h>

#include "error.h"
#include "init.h"
#include "parseopts.h"
#include "runopts.h"
#include "setup.h"
#include "stats.h"

#define EL_FILENAME_MAX 255
char InFile[EL_FILENAME_MAX + 1];

main(int argc, char** argv) {
  int err = 0;
  ElleRunFunc init;
  extern int InitSS(void);

  ElleInit();

  ElleSetInitFunction(InitSS);

  ElleSetOptNames("XVelocity", "unused", "unused", "unused", "unused", "unused", "unused", "unused", "unused");

  if (err = ParseOptions(argc, argv))
    OnError("", err);

  ElleSetSaveFileRoot("manuel");

  if (ElleDisplay())

    /*
     * run init and run functions
     */
    StartApp();
}
