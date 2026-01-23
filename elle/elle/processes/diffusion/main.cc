
/*
 *  main.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diffusion.h"
#include "error.h"
#include "init.h"
#include "parseopts.h"
#include "runopts.h"
#include "setup.h"

main(int argc, char** argv) {
  int err = 0;
  extern int InitDiffusion(), RunDiffusion(), CleanDiffusion();

  ElleInit();

  ElleSetInitFunction(InitDiffusion);
  ElleSetRunFunction(RunDiffusion);
  ElleSetExitFunction(CleanDiffusion);

  if (argc > 1) {
    if (err = ParseOptions(argc, argv))
      OnError("", err);
  } else
    ElleSetStages(100);

  if (strlen(ElleFile()) == 0)
    ElleSetFile("testsw1.poly");

  if (ElleDisplay())

    StartApp();

  return (0);
}
