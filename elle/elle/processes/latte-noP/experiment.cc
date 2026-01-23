
/******************************************************
 * New Experiment Class for Elle/Latte
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

  // here we can define more udata values that are read directly into the script
  // from the interface this means for these values we do not have to recompile
  // which we normally have to do for a change in values.

  strength = (int)udata[1];

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
    // -u 1 fracturing
    // -u 2 exchange reaction diffusion

    //--------------------------------------------------

    switch (process) {

      /*********************************************************************
       *
       *	PROCESSES Initialization
       *
       * so you can set up as many cases as you want for your
       * different ideas. That has the advantage that we dont have
       * to delete old ones when we want to make new stuff
       * At the moment only case 2 is setup for a start
       *
       **********************************************************************/
    case 0:
      // SetUpFromFile();
      break;

    case 1: // fracturing
      cout << "Fracturing" << endl;
      cout << "Lattice Version 2.0, 2004/5" << endl;

      Activate_Lattice(); // construct the lattice
      // SetPhase(0.0,0.0,1.0,1.3);	// set a distribution of breaking
      // strengths (linear)
      // SetGaussianSpringDistribution(1.0, 0.1); // set a distribution of
      // Youngs Moduli
      getTimeForCsv(experiment_time, 1);
      SetVariationBreakingThreshold(0.6);
      WeakenAll(0.4, 1.0, 1.5);
      // WeakenHorizontalBox(0.40,0.41,0.0,0.5,0.01,1.0,1.0,1.0);
      // WeakenHorizontalBox(0.59,0.60,0.5,1.0,0.01,1.0,1.0,1.0);
      ChangeRelaxThreshold(0.5);
      // MakeGrainBoundaries(1.0, 0.8);		// Make grain boundaries weaker
      Only_Extension(true);
      // Gravity(1000, 2700);
      // SetFracturePlot(50,1); 	// plot fractures after 50 bonds broken

      break;

    case 2: // exchange Daniel and Simon 2021

      cout << endl << "exchange reaction and diffusion" << endl;
      cout << "Latte 2021" << endl;

      Activate_Lattice(); // this is building up the original lattice object in
                          // C++

      //---------------------------------------------------------------------------------
      // this function Initialized the lattice
      //
      // first number is background fluid pressure, second number is scale (x)
      // in m third number is initial concentration of all nodes and last number
      // is time factor for the pressure diffusion. just to make it simpler we
      // define a background time in seconds as time_step and a background size
      // in meters potentially we should determine a time scale from the
      // explicit diffusion equation to optimize the solution and stability all
      // functions are in lattice.cc which then calls the fluid lattice
      // functions in fluid_lattice.cc
      //---------------------------------------------------------------------------------

      time_step = 0.5; // in seconds
      size = 0.01;     // in meters

      Initialize_Fluid_Lattice(10, size, 0.0,
                               time_step); // setting up the fluid

      // now we set a number of parameters, but initally a lot of these are not
      // needed

      SetDistributionPorosity(
          1.0); // this sets a distribution on the porosity, not used right now

      // variation of breaking threshold
      SetVariationBreakingThreshold(
          0.6); // for fracturing once we include that, background distribution

      // the next one is just a general adjustment of elasticity/viscosity
      // parameters, also not used right now change overall parameters, first
      // elastic constant, second viscosity, third breaking strength
      WeakenAll(0.8, 1.0, 1.0);

      // angular forces or not, this makes the solution easier for fracturing,
      // also not used right now
      Only_Extension(true);

      // This function can change the porosity in a box of xmin xmax ymin ymax,
      // change elastic/viscous parameters and the final number changes the
      // overall porosity. Also something for later
      WeakenHorizontalParticleLayer(0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.08);

      //------------------------------------------------------------
      // now comes the real stuff for the exchange and diffusion
      //------------------------------------------------------------

      // Solid state diffusion. The number is the diffusion coefficient for the
      // solid state diffusion, even though at the moment they are reset We
      // should clean this up a little bit

      Initiate_Mineral(0.0000000000001);

      //-------------------------------------------------------------------------------
      // the function set mineral defines what is garnet and what is biotite for
      // the moment. Garnet is 1 and Biotite 2 as last number. The first one is
      // the grain number from Elle, which you can see when you turn on Flynn
      // number
      //
      //--------------------------------------------------------------------------------

      Set_Mineral(107, 1);
      Set_Mineral(30, 2);

      // So this function defines the fluid diffusion coefficients, the first
      // one is the grains and the second one the grain boundary the new
      // function only uses the second one, there is no diffusion at all in the
      // grains for the grain boundary diffusion

      Set_Fluid(0.00000000000005, 0.000000005);

      // this is for relaxing the model elastically. Since nothing happens in
      // that respect this does nothing right now.
      Relaxation();

      break;

    case 3: // ridges

      cout << "ridges" << endl;
      cout << "Lattice Version 2022" << endl;

      Activate_Lattice(); // construct the lattice
      getTimeForCsv(experiment_time, 1);
      SetVariationBreakingThreshold(0.5); // Variation on breaking strings (<1),
                                          // before (0.7, 0.5), 0.3 probieren
      // SetGaussianSpringDistribution (1.0, 0.4);

      Initialize_Fluid_Lattice(1000000, 1000, 0,
                               400.0); // before (1000000, 1000, 0, 4000.0)
                                       // (pressure, scale, conc, timefactor)

      ReleaseBoundaryParticlesY(1, 0); // 1st: break? 2nd: elastic wall

      // SetWallBoundariesX(0, 1, 1, 2);

      WeakenAll(0.5, 1 * 1e20, 40); // WeakenAll(0.2,1,2); auf 0.1

      InitHeat(200, 5); // 200 (temperature, critical age)

      // InitDriftList();

      // WeakenHorizontal Slot(49, 51, 30, 70, 1.1, 1, 1) //ymin, ymax, xmin,
      // xmax, constant, viscosity, break strength

      ActivateSheet(1 * 1e20, 2 * 1e9); // das hier

      // Initialize_Heat_Box(42, 44, 0,60, 1200, 1);
      // Initialize_Heat_Box(29, 31, 0, 100, 1200, 1);

      Initialize_Heat_Box(29, 31, 75, 100, 1200, 1);
      Initialize_Heat_Box(44, 46, 50, 75, 1200, 1);
      Initialize_Heat_Box(29, 31, 25, 50, 1200, 1);
      Initialize_Heat_Box(44, 46, 0, 25, 1200, 1);

      CoolingDirectExp(0.3, 1.0, 1.0);

      WeakenHorizontalBox(0, 1, 0.0, 0.25, 1.0, 1.0, 100, 1);
      WeakenHorizontalBox(0, 1, 0.75, 1.0, 1.0, 1.0, 100, 1);

      // Initialize_Heat_Box(43, 47, 0, 100, 1200, 2);
      // Cooling(32);

      Loose_Unodes(0.2, 0.8);

      // ChangeRelaxThreshold(1.0);
      ChangeRelaxThreshold(5.0);

      // SetFracturePlot(5,1);

      break;

    case 4: // instability

      Activate_Lattice();
      getTimeForCsv(experiment_time, 1);
      SetWallBoundaries(0, 0.1);
      SetGaussianRateDistribution(2.0, 0.01);
      //----------------------------------------------------
      // change the Youngs modulus of particles to make
      // them stiffer (*4.0). Second and third number are
      // breaking strength and viscosity, 1.0 means no
      // change of these parameters
      //----------------------------------------------------
      WeakenAll(2.0, 1.0, 1000.0);
      //----------------------------------------------------
      // Set some mineral paramters, 1 means quartz as
      // mineral, sets the molecular volume and the
      // surface free energy
      //---------------------------------------------------
      Set_Mineral_Parameters(2);
      //----------------------------------------------------
      // gives the x dimension of the Elle box in meters
      //----------------------------------------------------
      Set_Absolute_Box_Size(0.001);
      //----------------------------------------------------
      // set the time for one deformation step. 6000 years
      // 4 means years
      //----------------------------------------------------
      Set_Time(0.1, 1);
      // DissolveYRow(0.0,0.2,false);
      DissolveXRow(0.90, 1.1);
      //----------------------------------------------------
      // dissolve initially one horizontal row of particles
      // in the middle of the Elle box
      // numbers are min and max y value
      //----------------------------------------------------
      Set_Fluid_Pressure(0.01);
      Set_Concentration(1, 2000);
      Make_Concentration_Box(2.0, 2000);
      Set_Dis_Time(10);
      ChangeRelaxThreshold(0.5);
      break;

    case 27: // Poroelastic Deformation (Hydrofracturing) new code for melt in
             // here
      cout << endl << "new fluid" << endl;
      cout << "Latte 2016" << endl;

      Activate_Lattice();

      Initialize_Fluid_Lattice(
          100, 1000, 0.1,
          2000.0); // pressure, scale, concentration, time factor

      SetDistributionPorosity(2.0);

      Fluid_Parameters(0, 1, 0.0000001);

      SetVariationBreakingThreshold(0.6);
      WeakenAll(1.0, 1.0, 1.0);

      Only_Extension(true);
      ChangeRelaxThreshold(0.05);
      Gravity(1000, 2700); // top depth, density
      getTimeForCsv(experiment_time, 10);
      Relaxation();

      // WeakenHorizontalParticleLayer(0.0,1.0,0.5, 0.55, 1.0, 1.0, 1.0,1.17);

      Hydrostatic_Fluid(1000);

      Set_Calcite();

      cout << "applied gravity " << endl;
      UpdateElle(); // and update the interface of Elle

      break;
    case 5: // ridges

      cout << "ridges" << endl;
      cout << "Lattice Version 2022" << endl;

      Activate_Lattice(); // construct the lattice
      getTimeForCsv(experiment_time, 1);
      SetVariationBreakingThreshold(0.5); // Variation on breaking strings (<1),
                                          // before (0.7, 0.5), 0.3 probieren
      // SetGaussianSpringDistribution (1.0, 0.4);

      Initialize_Fluid_Lattice(1000000, 1000, 0,
                               400.0); // before (1000000, 1000, 0, 4000.0)
                                       // (pressure, scale, conc, timefactor)

      ReleaseBoundaryParticlesY(1, 0); // 1st: break? 2nd: elastic wall

      // SetWallBoundariesX(0, 1, 1, 2);

      WeakenAll(0.5, 1 * 1e20, 40); // WeakenAll(0.2,1,2); auf 0.1

      InitHeat(200, 5); // 200 (temperature, critical age)

      InitDriftList(100);

      // WeakenHorizontal Slot(49, 51, 30, 70, 1.1, 1, 1) //ymin, ymax, xmin,
      // xmax, constant, viscosity, break strength

      ActivateSheet(1 * 1e20, 2 * 1e9); // das hier

      // Initialize_Heat_Box(42, 44, 0,60, 1200, 1);
      Initialize_Heat_Box(48, 52, 5, 95, 1200, 10);

      // Initialize_Heat_Box(29, 31, 75,100, 1200, 1);
      // Initialize_Heat_Box(44, 46, 50,75, 1200, 1);
      // Initialize_Heat_Box(29, 31, 25,50, 1200, 1);
      // Initialize_Heat_Box(44, 46, 0,25, 1200, 1);

      CoolingDirectExp(0.03, 0.50, 1.0);

      WeakenHorizontalBox(0, 1, 0.0, 0.25, 1.0, 1.0, 100, 1);
      WeakenHorizontalBox(0, 1, 0.75, 1.0, 1.0, 1.0, 100, 1);

      // Initialize_Heat_Box(43, 47, 0, 100, 1200, 2);
      // Cooling(32);

      Loose_Unodes(0.2, 0.8);

      // ChangeRelaxThreshold(1.0);
      ChangeRelaxThreshold(5.0);

      // SetFracturePlot(5,1);

      break;

    case 6: // fracturing
      cout << "Fracturing Deformation Bands" << endl;
      cout << "Lattice Version 2.0" << endl;

      Activate_Lattice(); // construct the lattice

      // SetGaussianSpringDistribution(1.0, 0.1); // set a distribution of
      // Youngs Moduli
      getTimeForCsv(experiment_time, 1);
      SetVariationBreakingThreshold(0.6);
      // WeakenAll(1,1.0,1.0);
      // WeakenHorizontalBox(0.40,0.41,0.0,0.5,0.01,1.0,1.0,1.0);
      // WeakenHorizontalBox(0.59,0.60,0.5,1.0,0.01,1.0,1.0,1.0);
      ChangeRelaxThreshold(0.5);
      MakeGrainBoundaries(1.0, 0.8); // Make grain boundaries weaker

      // WeakenGrain (int nb, float constant, float visc, float break_strength)
      // for (i = 1; i< 1000; i += 3)
      //{
      // WeakenGrain (i, 0.00001, 0, 1);
      //}
      WeakenGrain(14, 0.00001, 0, 1);
      WeakenGrain(0, 0.00001, 0, 1);
      WeakenGrain(10, 0.00001, 0, 1);
      WeakenGrain(12, 0.00001, 0, 1);
      WeakenGrain(9, 0.00001, 0, 1);

      // Only_Extension(true);
      //  Gravity(1000, 2700);
      SetFracturePlot(50, 1); // plot fractures after 50 bonds broken

      break;
    }
  } else
    cout << "no file open ! " << endl; // no input file in function call, can be
                                       // opened from the interface

  // UpdateElle(); // this function writes to the interface for plotting
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

      getTimeForCsv(experiment_time, 2);

      // if (experiment_time < 84)
      {
        //	DeformLattice(0.0001, 0);
      }
      // else
      {
        // DeformLatticeExtend(0.0001, 0);
      }
      DeformLatticePureShear(0.0001, 0);
      Get_Aperture();
      // else
      // DeformLatticePureShear(0.001,1);

      // DeformLatticePureShearAreaChange (0.00001, 0.00001, 1);
      // ViscousRelax(1, 1e15);
      UpdateElle();

      break;

      //---------------------------------------------------------------
      // this is the diffusion and exchange reaction function now
      // These function are all in lattice.cc
      //----------------------------------------------------------------

    case 2: // exchange reaction and diffusion

      // the first function is the grain boundary diffusion
      // the input is the time scale and length scale
      // we should let this function determine the time step
      // to optimize the process

      Calculate_Diffusion_Fluid_Explicit(time_step, size);

      // Second is the exchange, this needs input
      // it is given a reaction constant, a time setp and a size
      // I suspect the first number does not do anything
      // Can all the adjusted once we have the formulas
      // its only dumping stuff for five steps and then lets it
      // diffuse, also just a dummy for testing

      if (i < 5)
        Exchange_Solid_Fluid(1, time_step, size);

      // last comes the solid state diffusion
      // again gets the time and length scale as size

      Calculate_Diffusion_Solid(time_step, size);

      // the relaxation function is turned of, this would be for elastic
      // deformation

      // Relaxation();

      UpdateElle(); // writes to the interface

      break;

    case 3: // ridges

      getTimeForCsv(experiment_time, 2);

      // DeformLatticedriftmidnew(-0.0005, 0, 0.5, 0.4, 0.58, 0.58, 0.4);

      // DeformLatticeExtension(0.0005,0,1.0,0.001); //slower? alt:0.0005, shear
      // > 1.0? kk++; if (kk == 20) //etwas höher, 40?
      //{
      //	iii++;
      // kk = 0;
      // cout << "moved"<< endl;
      //}

      // Initialize_Heat_Box(48.5+iii, 51.5+iii, 0, 100, 1200, 20);
      // Cooling(32);
      // Initialize_Heat_Box(46+iii, 49+iii, 0, 46, 1200, 20);
      // Cooling(32);
      // Initialize_Heat_Box(42, 48, 0, 41, 800);
      // Cooling(32);
      // Initialize_Heat_Box(51+iii, 54+iii, 54, 100, 1200, 20);
      // Cooling(32);
      // Initialize_Heat_Box(58+iii, 60+iii, 54, 100, 1200);
      // Cooling(32);
      // Initialize_Heat_Box(52, 58, 59, 100, 800);
      // Cooling(32);

      // kk++;
      // if (kk == 20) //etwas höher, 40?
      //{
      //	iii++;
      //	kk = 0;
      //	cout << "moved"<< endl;
      //}

      //	if (kk == 10)
      //{
      //	Initialize_Heat_Box(40, 45, 50,100, 1200, 1);
      //	Initialize_Heat_Box(55, 60, 0,50, 1200, 1);
      //	cout << "more heat" << endl;
      //	}

      // DeformLatticedriftmidnew3(-0.001, experiment_time , 0.75, 0.4, 0.45,
      // 0.45, 0.6); //-0.0005

      // DeformLatticedriftmidnew4(-0.001, experiment_time , 0.75, 0.3,0.5,
      // 0.45,0.25, 0.3, 0.45, 1.1); DeformLatticedriftmidnew4(-0.001,
      // experiment_time , 0.75, 0.3,0.5, 0.3,0.25, 0.3, 0.3, 0.5);

      DeformLatticedriftmidnew4(-0.001, experiment_time, 0.75, 0.3, 0.5, 0.45,
                                0.25, 0.3, 0.45, 0.0);

      // DeformLatticedriftmidnew5(-0.001, experiment_time , 0.75, 0.3,0.25,
      // 0.25,0.25, 0.55, 0.45, 0.0);

      // DeformLatticedrift6(-0.001, experiment_time , 1.0, 0.3, 0.1, 0.7);

      // Relaxation();

      Find_Empty_Box(experiment_time, 1200);

      // space in meters for one box so resolution 200 = 20km for whole box
      // so with 0.001 deformation it opens 20m per step
      // for time = 10 and dif factor = 10 spreading rate is 20cm/year so fast
      // slow dif factor has to be 100

      solve_Heat(100, 10, 0.000004, 1200, 5,
                 10); // space, time (years was 20), dif const; 0.000001, temp,
                      // critAge, diffusion factor

      CoolingDirectExp(0.02, 1.0,
                       1.0); // that is the thermal weakening/strenghtening
      // ShrinkTemperature(1);
      WeakenHorizontalBox(0, 1, 0.0, 0.3, 1.0, 1.0, 100, 1);
      WeakenHorizontalBox(0, 1, 0.8, 1.0, 1.0, 1.0, 100, 1);
      SetVariationBreakingThreshold(0.5); // was 0.5

      Relaxation();
      RelaxSheet(1e2);
      cout << "viscous" << endl;
      ViscousStep(100000); // time up = less viscous?

      // RelaxSheet(1e5);
      cout << "healing" << endl;
      HealingTMol(1.8, 0.01, experiment_time, 1.0, 1.0); // vielleicht doch
                                                         // dazu?

      Dump_Young(10040, 10160);

      UpdateElle();

      // Relaxation();

      break;

    case 4:
      if (experiment_time < 1) {
        DeformLattice(0.001, 1);
        Relaxation();
        Relaxation();
      }
      GrowthDissolutionTwo(0, 1, 0, experiment_time);
      Relaxation();
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
        // Fluid_Insert_Random_Node(0.4, 0.5, 4000000);

        Calculate_Fluid_Pressure(9800000, 0, 1, 0, 1, 0.0);

      Adjust_Gravity();

      Relaxation();

      UpdateElle();

      break;
    case 5: // ridges random

      getTimeForCsv(experiment_time, 2);

      // DeformLatticedriftmidnew4(-0.001, experiment_time , 0.75, 0.3,0.5,
      // 0.45,0.25, 0.3, 0.45, 0.0);

      // DeformLatticedriftmidnew5(-0.001, experiment_time , 0.75, 0.3,0.25,
      // 0.25,0.25, 0.55, 0.45, 0.0);

      DeformLatticedrift6(-0.001, experiment_time, 1.0, 0.3, 0.1, 0.7);

      // Relaxation();

      Find_Empty_Box(experiment_time, 1200);

      // space in meters for one box so resolution 200 = 20km for whole box
      // so with 0.001 deformation it opens 20m per step
      // for time = 10 and dif factor = 10 spreading rate is 20cm/year so fast
      // slow dif factor has to be 100

      solve_Heat(100, 10, 0.000004, 1200, 5,
                 10); // space, time (years was 20), dif const; 0.000001, temp,
                      // critAge, diffusion factor

      CoolingDirectExp(0.03, 0.50,
                       1.0); // that is the thermal weakening/strenghtening
      // ShrinkTemperature(1);
      WeakenHorizontalBox(0, 1, 0.0, 0.3, 1.0, 1.0, 100, 1);
      WeakenHorizontalBox(0, 1, 0.8, 1.0, 1.0, 1.0, 100, 1);
      SetVariationBreakingThreshold(0.5); // was 0.5

      Relaxation();
      RelaxSheet(1e2);
      cout << "viscous" << endl;
      ViscousStep(1000000); // time up = less viscous?

      // RelaxSheet(1e5);
      cout << "healing" << endl;
      HealingTMol(1.8, 0.01, experiment_time, 1.0, 1.0); // vielleicht doch
                                                         // dazu?

      Dump_Young(10040, 10160);

      UpdateElle();

      // Relaxation();

      break;

    case 6: // fracturing Deformation bands

      getTimeForCsv(experiment_time, 2);

      // if (experiment_time < 84)
      { DeformLattice(0.0001, 0); }
      // else
      {
        // DeformLatticeExtend(0.0001, 0);
      }
      // DeformLatticePureShear(0.0001, 0);
      // Get_Aperture();
      // else
      // DeformLatticePureShear(0.001,1);

      // DeformLatticePureShearAreaChange (0.00001, 0.00001, 1);
      // ViscousRelax(1, 1e15);
      UpdateElle();

      break;
    }
    experiment_time++;
    Set_TimeFrac(experiment_time);
  }
}
