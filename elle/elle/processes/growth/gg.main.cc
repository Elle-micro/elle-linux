
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
  extern int InitGrowth(void);

  ElleInit();

  ElleSetInitFunction(InitGrowth);

  if (err = ParseOptions(argc, argv))
    OnError("", err);
  // ES_SetstatsInterval(50);

  ElleSetSaveFileRoot("growth");

  if (ElleDisplay())

    /*
     * run init and run functions
     */
    StartApp();

  return (0);
}
