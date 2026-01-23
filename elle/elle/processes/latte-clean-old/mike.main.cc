
/*
 *  main.cc
 */

#include <cstdio>
#include <cstdlib>

#include "error.h"
#include "file.h"
#include "init.h"
#include "lattice.h"
#include "parseopts.h"
#include "runopts.h"
#include "setup.h"

main(int argc, char** argv) {
  int err = 0;
  extern int InitSetMike(void);

  //*-----------------------------------------
  //* initialise
  //*-----------------------------------------

  ElleInit();

  //*--------------------------------------------------
  //* set the function to the one in your process file
  //*--------------------------------------------------

  ElleSetInitFunction(InitSetMike);

  if (err = ParseOptions(argc, argv))
    OnError("", err);

  //*-------------------------------------
  //* set up the X window
  //*-------------------------------------

  if (ElleDisplay())

    //*-------------------------------------------------------
    //* set the base for naming statistics and elle files
    //*------------------------------------------------------

    ElleSetSaveFileRoot("mike");

  //*--------------------------------------------------------------
  //* run your initialisation function and start the application
  //*--------------------------------------------------------------

  StartApp();

  CleanUp();

  return (0);
}
