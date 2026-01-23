
/******************************************************
 * New Experiment Class for Elle/Latte
 *
 * works now mainly with the lattice spring code
 *
 * Daniel and Till 2005
 *
 * Latte Version 2.0
 * koehn_uni-mainz.de
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
#include "attrib.h" // enth�lt struct coord (hier: xy)
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
//
// -----------------------------------------------------------------------------

{
  cout << "Oh, my experiment starts !" << endl;
  setlocale(LC_ALL, "en_US");
  experiment_time = 0;
  // press_cal_time = cal_time;
}

double time_step;
double nonlinear_time;
double lattice_size_m;

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
  int strength;

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
  process = (int)udata[0]; // get this data from Elle
  strength = (int)udata[1];

  //----------------------------------------------------
  // processes are defined by simple integers
  //----------------------------------------------------
  wrapping = false;

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
    //--------------------------------------------------

    switch (process) {

      /*********************************************************************
       *
       *	FRACTURE PROCESSES
       *
       **********************************************************************/
    case 0:
      // SetUpFromFile();
      break;

    case 27: // Poroelastic Deformation (Hydrofracturing) new code for melt in
             // here
      cout << endl << "new fluid" << endl;
      cout << "Latte 2016" << endl;

      Activate_Lattice();

      // Loose_Springs(true);
      // Shift_Particles(0.1);
      // Shift_Particles(0.1);
      // Shift_Particles(0.1);
      // Adjust_Springs_Particles();

      Initialize_Fluid_Lattice(100, 1000, 0.1, 2000.0);

      SetDistributionPorosity(2.0);

      Fluid_Parameters(0, 1, 0.0000001);

      SetVariationBreakingThreshold(0.6);
      WeakenAll(0.8, 1.0, 1.0);

      Only_Extension(true);

      Gravity(100, 2700); // depth, density
      getTimeForCsv(experiment_time, 10);
      Relaxation();

      // WeakenHorizontalParticleLayer(0.0,1.0,0.5, 0.55, 1.0, 1.0, 1.0,1.17);

      Hydrostatic_Fluid(100);

      // Set_Calcite();

      cout << "applied gravity " << endl;
      UpdateElle(); // and update the interface of Elle

      break;

    case 32: // advection-diffusion-reaction  // large scale fluid movement with
             // phreeqc as reaction engine
      cout << endl << "new fluid react coupling with PhreeqC" << endl;
      cout << "Latte 2021" << endl;

      Activate_Lattice();

      getTimeForCsv(experiment_time, 1);

      //---------------------------------------------------------------------------------
      // this function Initialized the lattice
      //
      // first number is background fluid pressure, second number is scale (x)
      // in m third number is initial concentration of all nodes and last number
      // is time factor for the pressure diffusion.
      //---------------------------------------------------------------------------------

      time_step = 1.0;      // in seconds
      lattice_size_m = 5.0; // in meters

      Initialize_Fluid_Lattice(10, lattice_size_m, 0.01,
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

      // WeakenHorizontalParticleLayer(0.2, 0.95, 0.5,
      // 0.55, 1.0, 1.0, 1.0,1.07);

      WeakenHorizontalParticleLayer(0.0, 0.4, 0.0, 0.3, 1.0, 1.0, 1.0, 1.07);

      WeakenHorizontalParticleLayer(0.6, 1.0, 0.0, 0.3, 1.0, 1.0, 1.0, 1.07);

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
      SetHydroChem(
          "seawater.txt", "brine_2.txt",
          "phreeqc.dat"); // first argument is pore fluid, second arguments is
                          // the infilatrating fluid ,and teh last is teh
                          // databse that should be used

      UpdateElle();

      break;
    }
  } else
    cout << "no file open ! " << endl; // no input file in function call, can be
                                       // opened from the interface

  // UpdateElle();
  // ElleUpdateDisplay();
  // cout << "update display" << endl;
}

/******************************************************************
 * A runfunction for the Experiment
 *
 *	In this function each process (again specified by -u in a switch
 *   command) is executed.
 *
 *  Latte version 2.0, 2005/6
 ******************************************************************/

void Experiment::Run() {
  //----------------------------------------------------
  // some local variables
  //------------------------------

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
  //--------------------------------------------------
  // loop through the time steps
  //--------------------------------------------------

  for (i = 0; i < time; i++) // cycle through stages
  {
    cout << "time step:" << experiment_time << endl;
    j++;
    switch (process) {
    case 0:
      // RunFromFile(experiment_time);
      break;

    case 1: // fracturing
            // if (experiment_time < 5)

      getTimeForCsv(experiment_time, 2);

      if (experiment_time < 1) {
        DeformLattice(0.0001, 0);
      } else {
        DeformLatticePureShear(0.0001, 0);
      }

      Get_Aperture();

      UpdateElle();
      break;

    case 27:
      //--------------------------------------------------------
      // A is healing
      //
      // B is fracture age
      //
      //----------------------------------------------------------

      Fluid_Parameters(0, 1, 0.0000001);
      for (iii = 0; iii < 10; iii++)
        Fluid_Insert_Random_Node(0.4, 0.5, 10000000);

      Calculate_Fluid_Pressure(9800000, 0, 1, 0, 1, 0.0);

      Adjust_Gravity();

      Relaxation();

      UpdateElle();

      break;

    case 32: // Advection diffusion and reaction	large scale
      //--------------------------------------------------------
      // A is solid
      // Temperature is fluid pressure
      // B is concentration
      // C is fluid pressure gradient
      // Energy is darcy velocity in y
      // Dislocation density is porosity
      //----------------------------------------------------------

      // set fluid parameters reads back material properties
      // from the lattice (particles etc. setting porosity
      // and permeability)

      Fluid_Parameters(0, 1, 0.00001);

      if (j < 2)
        TestInitPhreeqc();

      if (j < 400) {
        nonlinear_time = Calculate_Fluid_Pressure_timeadjust(10, 1000000, 1, 1,
                                                             0, j * 50000);
        TransportPhreeqc(nonlinear_time, lattice_size_m, 1);
        GetOutputPhreeqc();
      } else {
        nonlinear_time = Calculate_Fluid_Pressure_timeadjust(10, 1000000, 1, 1,
                                                             0, 400 * 50000);
        TransportPhreeqc(nonlinear_time, lattice_size_m, 1);
        GetOutputPhreeqc();
      }

      // TESTING PHREEQC
      // TestPhreeqc(nonlinear_time*10, size);
      // Make_Dolomite(1000000, i, 0.2, 0.01);

      // Relaxation();
      AccumulateSolIndex();
      // if ((j % 50) == 0)
      //{
      // DumpSolIndex (50);
      // DumpSolIndexvertical (0.33);

      UpdateElle();

      break;
    }
    experiment_time++;
    Set_TimeFrac(experiment_time);
    // LE I moved this from lattice
    // -----------------------------------------------------------------
    // call ElleUpdate() an Elle function that updates the interface
    // then the new values will be plotted
    // if security stop is set and max number of picts is reached dont
    // plot pict anymore.
    // -----------------------------------------------------------------
    // if (!set_max_pict)
    //   {
    // if (j == 10)
    //{
    // cout << "interface" << endl;
    // ElleUpdate ();
    // j = 0;
    //}
    // }
    // else if (num_pict < max_pict)
    //  {
    //  cout << "interface" << endl;
    // ElleUpdate ();
    // max_pict = max_pict + 1;
    //}
  }
}
