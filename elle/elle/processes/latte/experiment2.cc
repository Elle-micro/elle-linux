/******************************************************
 * Experiment Class for Elle/Latte
 *
 * works now mainly with the lattice spring code
 *
 * Daniel and Till 2005
 *
 * Latte Version 2.0
 * koehn_uni-mainz.de
 *
 * Latte Version 3.0
 * Erlangen 2021
 *
 * Clean Version 4.0
 * Erlangen 2025
 *
 ******************************************************/

// ------------------------------------
// system headers
// ------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <locale.h>
#include <vector>

// ------------------------------------
// elle headers
// ------------------------------------

#include "experiment.h"
#include "unodes.h" // include unode funct. plus undodesP.h

// for c++
#include "attrib.h" // enthlt struct coord (hier: xy)
#include "attribarray.h"
#include "attribute.h"
#include "check.h"
#include "convert.h"
#include "display.h"
#include "error.h"
#include "file.h"
#include "general.h"
#include "interface.h"
#include "nodes.h"
#include "nodesP.h"
#include "polygon.h"
#include "runopts.h"
#include "tripoly.h"
#include "update.h"

// have to define these in the new Elle version

using std::cout;
using std::endl;
using std::vector;

// CONSTRUCTOR
/*******************************************************
 * Dont really construct anything here at the moment
 *
 * The experiment class is directly called from the
 * Elle main function
 ********************************************************/

// ---------------------------------------------------------------
// Constructor of Experiment class
// ---------------------------------------------------------------

Experiment::Experiment()

// -----------------------------------------------------------------------------
//  first there is just some talking to the world plus
// -----------------------------------------------------------------------------

{
  cout << "Oh, my experiment starts !" << endl;
  cout << "Latte Version 4.0 2025" << endl;
  setlocale(LC_ALL, "en_US");
  experiment_time = 0;
}

double time_step;
double nonlinear_time;
double size;

/*************************************************************
 * Now we start with the initialization function
 * This function initialized some Elle basics
 * reads input file if one is there,
 * opens the interface
 * and starts the local initialization functions for the
 * desired process
 *
 *************************************************************/

void Experiment::Init() {
  //-------------------------------------------------
  // local variables
  //-------------------------------------------------

  char *infile; // input
  int err = 0;  // pass errors
  int i, j;     // counter
  float set;
  float visc;
  UserData udata; // Elle Structure for data from input
  int process;    // variable for the data (which process)
  int strength, value1, value2, value3, value4;

  //*-----------------------------------------------
  //* clear the data structures
  //*-----------------------------------------------

  ElleReinit();

  //*------------------------------------------------
  //* read the data
  //*-----------------------------------------------

  infile = ElleFile(); // input file specified by -i
  file = infile;
  ElleUserData(
      udata); // reads in data from the initial call of the program behind -u

  process = (int)udata[0]; // get this data from Elle, first number behind -u in
                           // Terminal during call

  // here we can define more udata values that are read directly into the script
  // from the interface this means for these values we do not have to recompile
  // which we normally have to do for a change in values. Here the numbers can
  // be changed during the call from the terminal directly remember that the
  // very first number after -u is the process below

  value1 = (int)udata[1]; // second number in Terminal
  value2 = (int)udata[2]; // third number in Terminal
  value3 = (int)udata[3]; // fourth number in Terminal
  value4 = (int)udata[4]; // fifth number in Terminal

  //----------------------------------------------------
  // processes are defined by simple integers
  //----------------------------------------------------
  wrapping = false; // not important right now

  if (strlen(infile) > 0) {
    //--------------------------------------------
    // only go in here if an input file is there
    //--------------------------------------------

    if (err = ElleReadData(infile))
      OnError(infile, err);

    //--------------------------------------------------
    // the switch function is used for manual input.
    // The processes can all be called from the
    // new interface.
    //
    // Without interface specify a process using the
    // -u command after calling the experiment.
    // ./my_experiment -u 1
    // -u 1 fracturing with aperture
    // -u 2 ridges with set geometry
    // -u 3 ridges with selforganized nucleation
    // -u 4 setup for grain fracturing
    //--------------------------------------------------

    switch (process) {

      /*******************************************************************************
       *
       *	PROCESSES Initialization
       *
       * so you can set up as many cases as you want for your
       * different ideas. That has the advantage that we dont have
       * to delete old ones when we want to make new stuff
       *
       * This will vary depending on the use of Latte, there should be
       * several versions for different processes, otherwise its too chaotic
       * this version has the fracturing and aperture (1) (Mathur et al., 2013)
       * the drifting for defined MORs and for the random ones (2) and (3)
       *(Hafermaas) and the new versions for fracturing grains in bands (Mathur
       *and Karimi)
       ********************************************************************************/

    case 0:
      // SetUpFromFile();
      break;

      /*******************************************************************************
       *	FRACTURING AND APERTURE (Mathur et al., 2023)
       *
       * This code simulates fracturing of a rock with or without grains and
       * measures the aperture of fractures during pure shear
       *
       * Deformation can include an initial compaction of the rock by uniaxial
       * compaction where normally nothing fractures. This simulates different
       * reservoir depth before the rock is deformed
       *
       ********************************************************************************/

    case 1: // fracturing

      cout << "Fracturing Aperture Mathur et al., 2023" << endl;
      cout << "Lattice Version 4.0, 2025" << endl;

      Activate_Lattice(); // construct the lattice

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 1);

      // function in lattice.cc that is used to set a Variation on the breaking
      // threshold of springs.
      SetVariationBreakingThreshold(0.6);

      // function to change properties of all particles (Young,
      // Viscosity,breaking_strength)
      WeakenAll(0.4, 1.0, 1.5);

      // change the relaxation threshold (<1 = more precise)
      ChangeRelaxThreshold(0.5);

      // MakeGrainBoundaries(1.0, 0.8);		// Make grain boundaries weaker

      // Allowing only tensile failure, not shear
      // typically bonds can break in shear and tension following
      // Sachau+Koehn2013

      Only_Extension(true);

      // SetFracturePlot(50,1); 	// plot fractures after 50 bonds broken

      break;

      /*******************************************************************************
       *
       *
       ********************************************************************************/

    case 2: // hydrofractures

      cout << "hydrofractures" << endl;
      cout << "Lattice Version 4.0, 2022-2025" << endl;

      Activate_Lattice(); // construct the lattice

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 1);

      UpdateElle();

      break;

      /*******************************************************************************
       *
       *
       ********************************************************************************/
    case 3: // reactions

      cout << "Reactions with PhreeqC" << endl;
      cout << "Lattice Version4.0, 2022-2025" << endl;

      Activate_Lattice(); // construct the lattice

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 1);

      Initialize_Fluid_Lattice(10, size, 0.01,
                               time_step); // setting up the fluid

      // this would set an initial fluid pressure gradient
      // Hydrostatic_Fluid(10);

      // this sets a random distribution on the porosity

      SetDistributionPorosity(0.010);

      // setting some vertical and horizontal fractures

      // variation of some elastic parameters and breaking strength

      // variation of breaking threshold
      SetVariationBreakingThreshold(0.6);

      // change overall parameters, first elastic constant, second viscosity,
      // third breaking strength
      WeakenAll(0.8, 1.0, 10.0);

      WeakenHorizontalParticleLayer(0.0, 0.4, 0.0, 0.3, 1.0, 1.0, 1.0, 1.07);

      WeakenHorizontalParticleLayer(0.6, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.07);

      // WeakenHorizontalParticleLayer(0.2, 0.95, 0.5,
      // 0.55, 1.0, 1.0, 1.0,1.07);

      // angular forces or not
      Only_Extension(true);

      // give fluid parameters, first number boundary condition, second number
      // cozeny grain size
      Fluid_Parameters(0, 1, 0.000001);

      Relaxation();
      Set_Calcite();

      cout << endl << " Phreeqc - Hydrogeochemical Modeling" << endl;
      SetResultMaps(
          "calcite,dolomite,density"); // set the desired output here. Currently
                                       // only for density and si's
      SetHydroChem("seawater.txt", "brine_2.txt",
                   "phreeqc.dat"); // first argument is pore fluid, second
                                   // arguments is the infilatrating fluid ,and

      UpdateElle();

      break;
    }
  } else
    cout << "no file open ! " << endl; // no input file in function call, can be
                                       // opened from the interface
}

/******************************************************************
 * A runfunction for the Experiment
 *
 *	In this function each process (again specified by -u in a switch
 *   command) is executed.
 *  	whereas the init function is called only ones, the run function
 * 	is called every time step. This is where the real action takes place
 *
 *  	Latte version 2.0, 2005/6
 *  	Latte version 3.0, 2021
 * 	Latte version 4.0, 2025
 ******************************************************************/

void Experiment::Run() {
  //----------------------------------------------------
  // some local variables
  //----------------------------------------------------

  int i, j, k, kk, ii, iii; // counter
  int time;                 // time
  int process;              // int for process
  UserData udata;           // elle defined user data (-u )

  ElleCheckFiles(); // Check the files

  ElleUserData(udata);     // get the usr data from elle
  process = (int)udata[0]; // process is user data 0

  ElleUpdateDisplay();

  //--------------------------------------------------
  // get the time from the interface
  //--------------------------------------------------

  time = EllemaxStages(); // number of stages

  j = 0;
  kk = 0;
  //--------------------------------------------------
  // loop through the time steps
  //--------------------------------------------------

  for (i = 0; i < time; i++) // cycle through stages
  {
    cout << "time step:" << experiment_time << endl;

    j++; // another time counter

    switch (process) {
    case 0:
      // RunFromFile(experiment_time);
      break;

    case 1: // fracturing

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 2);

      //==========================================================
      // For the first part of the simulation defined by
      // experiment time the deformation is uniaxial
      //==========================================================
      if (experiment_time < 84) {
        DeformLattice(0.0001, 0);
      }
      //==========================================================
      // For the rest of the simulation the deformation
      // is for example pure shear, so extending in x
      //==========================================================
      else {
        DeformLatticeExtend(0.0001, 0);
      }
      //==========================================================
      // function that determines the aperture of the cracks
      //==========================================================
      Get_Aperture();

      // and talk to the interface
      UpdateElle();

      break;

    case 2: // hydrofractures

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 2);

      UpdateElle();

      break;

    case 3: // reactions

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 2);
      // set fluid parameters reads back material properties
      // from the lattice (particles etc. setting porosity
      // and permeability)

      Fluid_Parameters(0, 1, 0.00001);

      if (j < 2)
        TestInitPhreeqc();

      if (j < 400) {
        nonlinear_time = Calculate_Fluid_Pressure_timeadjust(10, 1000000, 1, 1,
                                                             0.5, j * 50000);
        TransportPhreeqc(nonlinear_time, size, 1);
        GetOutputPhreeqc();
      } else {
        nonlinear_time = Calculate_Fluid_Pressure_timeadjust(10, 1000000, 1, 1,
                                                             0.5, 400 * 50000);
        TransportPhreeqc(nonlinear_time, size, 1);
        GetOutputPhreeqc();
      }

      // TESTING PHREEQC
      // TestPhreeqc(nonlinear_time*10, size);
      // Make_Dolomite(1000000, i, 0.2, 0.01);

      Relaxation();
      AccumulateSolIndex();
      // if ((j % 50) == 0)
      //{
      // DumpSolIndex (50);
      // DumpSolIndexvertical (0.33);

      // and talk to the interface.
      UpdateElle();

      break;
    }
    experiment_time++;
    Set_TimeFrac(experiment_time);
  }
}
