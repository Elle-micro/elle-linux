
/*
 *  main.c
 */

#include "split.h"

#include <stdio.h>
#include <stdlib.h>

#include "attrib.h"
#include "errnum.h"
#include "error.h"
#include "file.h"
#include "init.h"
#include "parseopts.h"
#include "runopts.h"
#include "setup.h"
#include "stats.h"
#include "string_utils.h"

float TotalTime;

main(int argc, char** argv) {
  int err = 0;
  ElleRunFunc init;
  extern int InitThisProcess(void);
  UserData userdata;

  ElleInit();

  ElleUserData(userdata);
  userdata[SMode] = 0; // Change default calculation mode
  ElleSetUserData(userdata);

  ElleSetOptNames("SplitMode", "unused", "unused", "unused", "unused", "unused", "unused", "unused", "unused");
  if (err = ParseOptions(argc, argv))
    OnError("", err);

  ElleSetInitFunction(InitThisProcess);

  ElleSetSaveFileRoot("split");

  if (ElleDisplay())

    StartApp();
  return (0);
}
