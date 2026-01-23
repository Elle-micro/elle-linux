
/*
 *  main.cc
 */

#include <cstdio>
#include <cstdlib>

#include "error.h"
#include "experiment.h"
#include "file.h"
#include "init.h"
#include "parseopts.h"
#include "runopts.h"
#include "setup.h"

Experiment* MyExperiment;

int main(int argc, char** argv) {
  int err = 0;
  UserData userdata;

  //*-----------------------------------------
  //* initialise
  //*-----------------------------------------

  ElleInit();

  // extern int Init_Experiment(void);

  MyExperiment = new Experiment;

  //*--------------------------------------------------
  //* set the function to the one in your process file
  //*--------------------------------------------------

  // ElleSetInitFunction(InitSetMike);

  // ElleSetInitFunction(MyExperiment->Init_Experiment);

  ElleSetInitFunction(Init_Experiment);

  ElleUserData(userdata);
  userdata[0] = 0; // Change default calculation mode
                   // -u 0 read experiment type and parameters from a zip file
                   //      e.g. if input elle file is res50.elle then read res50.zip
                   // -u 1 fracturing
                   // -u 2 fracture boudinage
                   // -u 3 expanding inclusions
                   // -u 4 shrinkage cracks
                   // -u 5 viscoelastic deformation
                   // -u 6 grooves
                   // -u 7 Stylolites
                   // -u 8 combine graingrowth and fractures
                   // -u 9 solid solid phase change
                   // -u 10 heat flow
                   // -u 11 pure grain growth
                   // -u 12 Lattice gas diffusion
                   // -u 13 Lattice gas fluid flow
  ElleSetUserData(userdata);

  ElleSetOptNames("Experiment", "unused", "unused", "unused", "unused", "unused", "unused", "unused", "unused");

  // Parse CLI options; on help/syntax error, print and exit early.
  err = ParseOptions(argc, argv);
  if (err) {
    OnError("", err);
    if (err == HELP_ERR)
      return 0; // --help requested; already printed usage via OnError
    return 1;   // syntax or other option error
  }

  //*-------------------------------------------------------
  //* set the base for naming statistics and elle files
  //*------------------------------------------------------

  ElleSetSaveFileRoot("my_experiment");

  //*-------------------------------------
  //* set up the X window
  //*-------------------------------------

  if (ElleDisplay())

    //*--------------------------------------------------------------
    //* run your initialisation function and start the application
    //*--------------------------------------------------------------

    StartApp();

  CleanUp();

  return (0);
}
extern int Init_Experiment(void) {

  //-----------------------------------------------------
  // Set the run function
  //-----------------------------------------------------

  ElleSetRunFunction(Run_Experiment);
  MyExperiment->Init();
  return 0;
}
extern int Run_Experiment() {
  MyExperiment->Run();
  return 0;
}
