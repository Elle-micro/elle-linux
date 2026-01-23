
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
  extern int InitThisProcess(void);

  /*
   * initialise
   */
  ElleInit();

  if (err = ParseOptions(argc, argv))
    OnError("", err);

  /*
   * set the function to the one in your process file
   */
  ElleSetInitFunction(InitThisProcess);

  ElleSetDisplay(0);

  ElleSetStages(1);

  /*
   * set the interval for writing to the stats file
  ES_SetstatsInterval(100);
   */
  ElleSetSaveFileRoot("stats");

  if (ElleDisplay())

    StartApp();

  return (0);
}
