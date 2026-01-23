
/*
 *  main.c
 */

#include "gbdiff.h"

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

main(int argc, char** argv) {
  int err = 0;
  UserData udata;
  extern int InitVM(void);

  ElleInit();

  ElleSetInitFunction(InitVM);

  ElleUserData(udata);
  udata[Kappa] = KappaDflt;
  udata[NodeId] = NO_NB;
  udata[AttribId] = CONC_A;
  ElleSetUserData(udata);
  ElleSetOptNames("Kappa", "Bnode Id", "Attrib Value", "Reset", "Attribute", "unused", "unused", "unused", "unused");

  if (err = ParseOptions(argc, argv))
    OnError("", err);

  ElleSetSaveFileRoot("gbdiff");

  if (ElleDisplay())

    /*
     * run init and run functions
     */
    StartApp();
  CleanUp();
}
