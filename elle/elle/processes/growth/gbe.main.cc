
/*
 *  main.cc
 */

#include <stdio.h>
#include <stdlib.h>

#include "error.h"
#include "init.h"
#include "parseopts.h"
#include "runopts.h"
#include "setup.h"
#include "stats.h"

main(int argc, char** argv) {
  int err = 0;
  extern int InitGBE_Growth(void);

  ElleInit();

  ElleSetInitFunction(InitGBE_Growth);

  if (err = ParseOptions(argc, argv))
    OnError("", err);
  // ES_SetstatsInterval(100);

  ElleSetSaveFileRoot("gbegrowth");

  if (ElleDisplay())

    /*
     * run init and run functions
     */
    StartApp();

  return (0);
}
