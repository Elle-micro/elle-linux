
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
      WeakenAll(0.8, 1.0, 1.5);
      // WeakenHorizontalBox(0.40,0.41,0.0,0.5,0.01,1.0,1.0,1.0);
      // WeakenHorizontalBox(0.59,0.60,0.5,1.0,0.01,1.0,1.0,1.0);
      ChangeRelaxThreshold(0.5);
      // MakeGrainBoundaries(1.0, 0.8);		// Make grain boundaries weaker
      Only_Extension(true);
      // Gravity(1000, 2700);
      // SetFracturePlot(50,1); 	// plot fractures after 50 bonds broken
      UpdateElle();
      break;

    case 2: // fracture boudinage
      cout << "Fracture Boudinage" << endl;
      cout << "Lattice Version 2.0, 2004/5" << endl;
      Activate_Lattice(); // contruct the lattice

      SetGaussianStrengthDistribution(
          2.0, 0.8); // set a distribution of breaking strengths (gaussian)
      WeakenAll(0.1, 1.0, 1.0); // Lower Youngs modulus of all grains

      WeakenHorizontalParticleLayer(
          0.20, 0.22, 0.0, 1.0, 10.0, 1.0, 1.0,
          1.0); // Make a layer at (ymin, ymax, Youngs Modulus, Viscosity,
                // breaking strength)
      WeakenHorizontalParticleLayer(0.9, 0.92, 0.0, 1.0, 10.0, 1.0, 1.0,
                                    1.0); // same as above

      SetFracturePlot(50, 1); // plot fractures after 50 bonds are broken

      break;

    case 3: // expanding inclusions
      cout << "Expanding Inclusions" << endl;
      cout << "Lattice Version 2.0, 2004/5" << endl;
      Activate_Lattice(); // construct the lattice
      // SetPhase(0.0,0.0,0.8,1.2); // Linear distribution of breaking strengths

      // variation of breaking threshold

      SetVariationBreakingThreshold(0.6);

      // change overall parameters, first elastic constant, second viscosity,
      // third breaking strength
      WeakenAll(0.8, 1.0, 1.0);

      // angular forces or not
      Only_Extension(false);

      ChangeRelaxThreshold(10.0);

      SetFracturePlot(20, 1); // plot fractures after 5 bonds are broken
      break;

    case 4: // shrinkage patterns
      cout << "Shrinkage Patterns" << endl;
      cout << "Lattice Version 2.0, 2004/5" << endl;
      Activate_Lattice();     // construct the lattice
      WeakenAll(0.7, 1, 3.0); // Weaken Youngs modulus of all grains

      SetVariationBreakingThreshold(0.6);

      // SetPhase(0.0,0.0,3.5,1.0); // Distribution of breaking strengths
      // (linear)
      MakeGrainBoundaries(1.0, 0.7);
      // SetWallBoundaries(1,1.0);
      // ReleaseBoundaryParticlesY(1,0);
      // ActivateSheet(1e16,1e9);
      // SetCosShrink(1.0);
      // SetCosViscSheetMod(0.8, 0.05);
      // SetSinAnisotropy(10, 1.2); // Sinusoidal anisotropy (vertical) of
      // Youngs Moduli SetGaussianStrengthDistribution(1.0,0.35);  //
      // Distribution of breaking strengths (gaussian)
      Only_Extension(true);
      SetFracturePlot(20, 1); // plot fractures after 5 bonds are broken
      ChangeRelaxThreshold(1.0);

      break;
      /**************************************************************************
       *
       *		VISCOELASTIC DEFORMATION AND FRACTURING
       *
       *************************************************************************/

    case 5: // viscoelastic
      cout << "visco-elastic" << endl;
      cout << "Lattice Version 2.0, 2004/5" << endl;
      Activate_Lattice(); // Construct the lattice
      //--------------------------------------------------
      // Set external walls around the box. No particles
      // can now leave the box.
      // External walls are only compressive
      //--------------------------------------------------
      SetWallBoundaries(1, 1.0);
      //--------------------------------------------------
      // Set a distribution on the breaking strengths of
      // all springs. Distribution can be Gaussian
      // around a mean strength or linear between
      // two endmembers. Values are distributed randomly
      //--------------------------------------------------
      // SetGaussianStrengthDistribution(2.5,0.8); // Set a Gaussian
      // distribution with a mean and a variance
      // Set a linear distribution around a mean value and with a certain width
      //--------------------------------------------------
      // Change youngs modulus, Viscosity and breaking
      // Strength of all particles
      // Youngs modulus of 1 is default, Scales as
      // 10 GPa in the model (this value is then a
      // reference).
      // Viscosity is by default e^20 if Youngs modulus 1
      // in the model is 10 GPa in reality.
      //--------------------------------------------------
      WeakenAll(0.1, 1, 1.0);
      SetPhase(0.0, 0.0, 1, 0.6);
      //----------------------------------------------------------------------------------
      // Change some initial values in order to have an anisotropy
      // a) Change values of one grain (or more)
      // b) induce a horizontal anisotropy in the form of hard flakes with
      // varying Moduli c) Insert a horizontal weaker or stronger layer
      //----------------------------------------------------------------------------------
      // WeakenGrain(7,10,10,1); // Change Youngs modulus and viscosity of grain
      // nb 7 SetAnisotropyRandom(5,5); // insert a horizontal anisotropy of
      // mica grains
      WeakenHorizontalParticleLayer(0.47, 0.52, 0.0, 1.0, 5.0, 5.0, 1.0,
                                    1.0); // Insert a horizontal layer
      //------------------------------------------------------------------------------
      // ChangeRelaxThreshold(0.1);//factor, 1=normal
      // SetFracturePlot(5,1);
      break;
      /*****************************************************************************
       *
       *		REACTIONS
       *
       ******************************************************************************/

    case 6: // grooves on free surfaces
      cout << "Dissolution Grooves" << endl;
      cout << "Phase_Lattice Version 2.0, 2004/5" << endl;
      Activate_Lattice();             // construct the lattice
      SetPhase(0.0, 0.0, 500.0, 0.8); // impossible to break (strength * 500)
      SetGaussianRateDistribution(
          2.0, 0.3); // Distribution on rate constants of reaction (gaussian)
      WeakenAll(8.0, 1.0, 1.0);  // Change Youngs modulus of all particles
      Set_Mineral_Parameters(1); // define a mineral here Quartz
      Set_Absolute_Box_Size(
          0.00008);            // Set the absolute Elle box size in meters
      Set_Time(6000.0, 4);     // set the time (here 6000 years)
      DissolveXRow(0.95, 1.1); // dissolve particles > xpos 0.95 to have a free
                               // interface on the right side
      SetWallBoundaries(0, 15.0);
      break;

    case 7: // Stylolites
      cout << "Stylolite Roughening" << endl;
      cout << "Phase_Lattice Version 2.0, 2004/5" << endl;
      Activate_Lattice();
      //----------------------------------------------------
      // give bonds very high strength (*500) to avoid
      // in this case fracturing during Stylolite growth
      //----------------------------------------------------
      // SetPhase(0.0,0.0,500.0,0.8);
      //----------------------------------------------------
      // Set a Gaussian distribution on the rate constants
      // of single particles, first mean value, second
      // deviation
      //----------------------------------------------------
      // SetGaussianRateDistribution(2.0,0.1);  // general background
      // distribution

      Set_Rate_Two_Phase(
          0.05, 0.3,
          1.0); // bimodel distribution, make some particles dissolve slower

      // Set_Rate_Three_Phase(0.05, 0.7, 1.0, 0.3, 0.9, 1.0);

      // Set_Rate_Two_Phase_Layer(1.0, 0.64, 0.9, 0.47, 0.53);		//
      // insert a layer that dissolves slower

      // Set_Rate_Two_Phase_Layer(1.0, 0.5, 0.9, 0.40, 0.43);
      // Set_Rate_Two_Phase_Layer(1.0, 0.5, 0.9, 0.3, 0.33);
      // Set_Rate_Two_Phase_Layer(1.0, 0.5, 0.9, 0.2, 0.23);
      // Set_Rate_Two_Phase_Layer(1.0, 0.5, 0.9, 0.1, 0.13);
      // Set_Rate_Two_Phase_Layer(1.0, 0.5, 0.9, 0.8, 0.83);
      // Set_Rate_Two_Phase_Layer(1.0, 0.5, 0.9, 0.6, 0.63);
      // Set_Rate_Two_Phase_Layer(1.0, 0.5, 0.9, 0.7, 0.73);
      // Set_Rate_Two_Phase_Layer(1.0, 0.5, 0.9, 0.9, 0.93);

      // Set_Rate_Two_Phase_Layer(1.0, 0.5, 0.9, 0.52, 0.54);		//
      // insert a layer that dissolves slower

      // Set_Rate_Two_Grains(0.01, 0.5);  // change larger elle grains

      //----------------------------------------------------
      // change the Youngs modulus of particles to make
      // them stiffer (*4.0). Second and third number are
      // breaking strength and viscosity, 1.0 means no
      // change of these parameters
      //----------------------------------------------------
      WeakenAll(1.0, 1.0, 500000.0);
      //----------------------------------------------------
      // Set some mineral paramters, 1 means quartz as
      // mineral; 2 is NaCl, 3 is Olivine, 4 is Calcite, sets the molecular
      // volume and the surface free energy
      //----------------------------------------------------
      Set_Mineral_Parameters(4);
      //----------------------------------------------------
      // gives the x dimension of the Elle box in meters
      //----------------------------------------------------
      Set_Absolute_Box_Size(0.1); // was 0.1
      //----------------------------------------------------
      // set the time for one deformation step. 6000 years
      // 4 means years
      //----------------------------------------------------
      Set_Time(500, 4); // was 40,4
      //----------------------------------------------------
      // dissolve initially one horizontal row of particles
      // in the middle of the Elle box
      // numbers are min and max y value
      //----------------------------------------------------
      // DissolveYRowSinus(0.49,0.5,true); //was 0.496
      DissolveYRow(0.498, 0.502, true); // was 0.496
      // ChangeRelaxThreshold(0.1);

      ChangeRelaxThreshold(10);

      // ChangeYoung(2);
      break;

    case 8:
      cout << "Combine Latte and GrainGrowth" << endl;
      cout << "Phase_Lattice Version 2.0, 2004/5" << endl;
      Activate_Lattice(); // Construct the lattice
      ElleInitNodeAttribute(ATTRIB_A);
      SetPhase(0.0, 0.0, 200.0, 1.0); // Set breaking strength distributio
      SetGaussianSpringDistribution(0.5,
                                    0.5); // set distribution on youngs moduli
      MakeGrainBoundaries(1.0, 0.5);      // define grain boundaries
      SetFracturePlot(1, 0);              // plot every fracture
      ElleAddDoubles(); // add some doubles for grain growth (Elle function)
      break;

    case 9: // phase change
      cout << "solid solid phase transformation, slow reaction, no distribution"
           << endl;
      /* Activate_MinTrans();
       SetReactions(450000.0, 12e9); // sets the activation energy (J/mol) for
       the reaction and the pressure barrier where the grain-boundary migration
       will start (Pa). Also calls setheatflowparameters() for olivine
       heat_distribution.SetHeatFlowParameters(4.2, 1005.0, 0.000001, 1000.0);
       // rho, c, diffusivity, boundary_condition (in °K) WeakenAll(20.0, 0.0,
       5000.0); AdjustConstantGrainBoundaries(); MakeGrainBoundaries(1.0,0.8);
       Set_Mineral_Parameters(3);
       Set_Absolute_Box_Size(0.005);
       Set_Time(120, 3); //(x,2):0=sek, 1=Stunden, 2=Tage, 3=Monate, 4=Jahre*/
      break;

    case 10:
      cout << "rifting" << endl;
      Activate_Lattice();

      Set_Absolute_Box(200000, 1);

      // SetPhase(1.0,0.0,1.6,1.1); //was 1.5 and 2.0

      SetVariationBreakingThreshold(0.6);

      WeakenAll(1.0, 1.0, 5.0);

      UnsetNobreakUpperRow();
      InitRealDensityYoung(3200, 1.8e11);

      SetRealDensityYoung(1.0, 0.0, 1.0, 0.9, 2750, 0.85e11);
      SetRealDensityYoung(1.0, 0.0, 0.9, 0.825, 2900, 1.08e11);

      MarkMohoParticles(0.825);

      SetRealDensityYoung(1.0, 0.0, 0.825, 0.3, 3200, 1.8e11);

      SetRealDensityYoung(0.6, 0.4, 1.0, 0.9, 2750, 2.2e11);   // Rwenzori block
      SetRealDensityYoung(0.6, 0.4, 0.9, 0.825, 2900, 2.2e11); // Rwenzori block

      // SetRealDensityYoung(0.6,0.4,0.875,0.7,2900,1.8e11);

      SetRealDensityYoung(1.0, 0.0, 0.3, 0.0, 2900, 1.08e11);

      SetRealDensityYoung(0.9, 0.1, 0.4, 0.3, 2900, 1.8e11);
      SetRealDensityYoung(0.8, 0.2, 0.5, 0.4, 2900, 1.8e11);
      SetRealDensityYoung(0.7, 0.3, 0.6, 0.5, 2900, 1.8e11);

      // SetGaussianYoungDistribution(1.0,0.1);
      InitInternalYoung(1.8e11);

      InitViscosity(1e23);
      set = 0.9;
      visc = 1e25;

      SetViscosity(1.0, 0.0, 1.0, 0.96, 1e25);   // oberste Kruste
      SetViscosity(1.0, 0.0, 0.96, 0.875, 1e23); // unterer Teil oberste Kruste
      SetViscosity(1.0, 0.0, 0.875, 0.85,
                   8e20); // oberer Teil Unterkruste was 2e21
      SetViscosity(1.0, 0.0, 0.85, 0.825,
                   8e20); // unterer Teil Unterkruste 8e21 duktil

      // SetViscosity(0.8,0.6,0.875,0.825,1e20);  // soft rift
      // SetViscosity(0.4,0.2,0.875,0.825,1e20);  // soft rift

      SetViscosity(1.0, 0.0, 0.825, 0.7, 1e21); // oberer Mantel
      // SetViscosity(0.7,0.3,0.825,0.7,5e21);
      SetViscosity(1.0, 0.0, 0.7, 0.3,
                   2e22); // unterer Mantel (Lithosphare untergrenze)
      SetViscosity(1.0, 0.0, 0.3, 0.0, 2e22); // Astenopshaere

      SetBreakingStrength(1.0, 0.0, 1.0, 0.98, 0.8);
      SetBreakingStrength(1.0, 0.0, 0.98, 0.96, 0.8);
      SetBreakingStrength(1.0, 0.0, 0.96, 0.825, 0.8); // was 0.8
      // SetBreakingStrength(1.0,0.7,1.0,0.0,1000);
      // SetBreakingStrength(0.3,0.0,1.0,0.0,1000);
      SetBreakingStrength(1.0, 0.0, 0.825, 0.5, 5);
      SetBreakingStrength(1.0, 0.0, 0.5, 0.0, 50000);

      // SetBreakingStrength(1.0,0.7,1.0,0.0,1000);

      cout << "in adjust" << endl;

      AdjustParticleConstants();

      cout << "out adjust" << endl;

      // ChangeRelaxThreshold(10);

      // SetFracturePlot(20,1);
      break;

    case 11:
      cout << "pure grain growth" << endl;
      ElleAddDoubles(); // add doubles (Elle functions)
      // ElleUpdate();      make a picture
      ElleUpdateDisplay(); // make a picture
      break;

    case 12:
      cout << " Lattice Gas" << endl;
      Activate_Lattice(); // contruct the lattice
      SetFluidLatticeGasRandom(
          0.01); // background contains some particles randomly distributed (1%)
      /*********************************************************
       * Specify some grains that have a concentration of 70 %
       *********************************************************/
      for (j = 0; j < HighestGrain(); j++) {
        SetFluidLatticeGasRandomGrain(0.4, j * 10);
      }
      break;

    case 13:
      cout << " Lattice Gas Flow" << endl;
      Activate_Lattice(); // contruct the lattice
      SetFluidLatticeGasRandom(
          0.005); // fluid density in the background (0.5 %)
      /*********************************************************
       * Specify some grains that are fracture walls
       *********************************************************/
      for (j = 0; j < HighestGrain(); j++) {
        SetWallsLatticeGas(j * 5);
      }
      break;

    case 15:
      Activate_Lattice();
      SetWallBoundaries(1, 1);
      SetGaussianRateDistribution(2.0, 0.3);
      //----------------------------------------------------
      // change the Youngs modulus of particles to make
      // them stiffer (*4.0). Second and third number are
      // breaking strength and viscosity, 1.0 means no
      // change of these parameters
      //----------------------------------------------------
      WeakenAll(8.0, 1.0, 100.0);
      //----------------------------------------------------
      // Set some mineral paramters, 1 means quartz as
      // mineral, sets the molecular volume and the
      // surface free energy
      //----------------------------------------------------
      Set_Mineral_Parameters(1);
      //----------------------------------------------------
      // gives the x dimension of the Elle box in meters
      //----------------------------------------------------
      Set_Absolute_Box_Size(0.0001);
      //----------------------------------------------------
      // set the time for one deformation step. 6000 years
      // 4 means years
      //----------------------------------------------------
      Set_Time(60.0, 4);
      DissolveYRow(0.0, 0.2, false);
      // DissolveXRow(0.8,1.1);
      //----------------------------------------------------
      //  dissolve initially one horizontal row of particles
      //  in the middle of the Elle box
      //  numbers are min and max y value
      //----------------------------------------------------
      Set_Fluid_Pressure(10);
      Set_Concentration();
      Make_Concentration_Box(1.0, 2);
      Set_Dis_Time(40);
      break;

    case 16:
      Activate_Lattice();
      SetWallBoundaries(0, 0.1);
      SetGaussianRateDistribution(2.0, 0.01);
      //----------------------------------------------------
      // change the Youngs modulus of particles to make
      // them stiffer (*4.0). Second and third number are
      // breaking strength and viscosity, 1.0 means no
      // change of these parameters
      //----------------------------------------------------
      WeakenAll(8.0, 1.0, 100.0);
      //----------------------------------------------------
      // Set some mineral paramters, 1 means quartz as
      // mineral, sets the molecular volume and the
      // surface free energy
      //---------------------------------------------------
      Set_Mineral_Parameters(2);
      //----------------------------------------------------
      // gives the x dimension of the Elle box in meters
      //----------------------------------------------------
      Set_Absolute_Box_Size(0.0001);
      //----------------------------------------------------
      // set the time for one deformation step. 6000 years
      // 4 means years
      //----------------------------------------------------
      Set_Time(6.0, 2);
      // DissolveYRow(0.0,0.2,false);
      DissolveXRow(0.8, 1.1);
      //----------------------------------------------------
      // dissolve initially one horizontal row of particles
      // in the middle of the Elle box
      // numbers are min and max y value
      //----------------------------------------------------
      Set_Fluid_Pressure(0.01);
      Set_Concentration();
      Make_Concentration_Box(1.1, 1);
      Set_Dis_Time(40);
      break;

    case 17:
      cout << "solid solid phase transformation, fast reaction, no distribution"
           << endl;
      Activate_MinTrans();
      SetReactions(
          380000.0,
          12e9); // sets the activation energy (J/mol) for the reaction and the
                 // pressure barrier where the grain-boundary migration will
                 // start (Pa). Also calls setheatflowparameters() for olivine
      heat_distribution.SetHeatFlowParameters(
          4.2, 1005.0, 0.000001,
          1000.0); // rho, c, diffusivity, boundary_condition (in °K)
      WeakenAll(20.0, 0.0, 5000.0);
      AdjustConstantGrainBoundaries();
      MakeGrainBoundaries(1.0, 0.8);
      Set_Mineral_Parameters(3);
      Set_Absolute_Box_Size(0.005);
      Set_Time(120, 3); //(x,2):0=sek, 1=Stunden, 2=Tage, 3=Monate, 4=Jahre*/
      break;

    case 18:
      cout << "solid solid phase transformation, high distribution" << endl;
      /* Activate_MinTrans();
       SetReactions(450000.0, 12e9); // sets the activation energy (J/mol) for
       the reaction and the pressure barrier where the grain-boundary migration
       will start (Pa). Also calls setheatflowparameters() for olivine
       heat_distribution.SetHeatFlowParameters(4.2, 1005.0, 0.000001, 1000.0);
       // rho, c, diffusivity, boundary_condition (in °K) SetPhase(20.0, 1.0,
       5000.0, 0.000000001); //spring constant wurde justiert (500.0 statt 0.0)
       AdjustConstantGrainBoundaries();
       MakeGrainBoundaries(1.0,0.8);
       SetGaussianYoungDistribution_2(1.0,0.5);
       Set_Mineral_Parameters(3);
       Set_Absolute_Box_Size(0.005);
       Set_Time(120, 3); //(x,2):0=sek, 1=Stunden, 2=Tage, 3=Monate, 4=Jahre*/
      break;

    case 19: // shrinkage patterns with viscous sheet
      cout << "Shrinkage Patterns, Viscoelastic sheet" << endl;
      cout << "Lattice Version 2.1, 2006" << endl;

      Activate_Lattice(); // construct the lattice

      // Loose_Springs(true);
      // Shift_Particles(0.1);
      // Shift_Particles(0.1);
      // Shift_Particles(0.1);
      // Adjust_Springs_Particles();
      // Connect_New_Lattice(true);
      // Loose_Springs(true);
      // Relaxation();
      // Connect_New_Lattice(true);

      WeakenAll(0.2, 1, 1.0); // Weaken Youngs modulus of all grains

      // SetGaussianStrengthDistribution(2.0,0.4);

      SetVariationBreakingThreshold(0.4);

      // WeakenHorizontalParticleLayer(0.47, 0.53,
      // 0.0, 1.0,  7.0, 1.0, 1.0,1.0);

      // SetPhase(0.0,0.0,1.5,1.5); // Distribution of breaking strengths
      // (linear)

      // ActivateSheet(1e19,1e10); // first visc, then young

      SetFracturePlot(20, 1); // plot fractures after 5 bonds are broken

      // SetSinAnisotropy(20, 1.2); // Sinusoidal anisotropy (vertical) of
      // Youngs Moduli SetGaussianStrengthDistribution(1.0,0.35);  //
      // Distribution of breaking strengths (gaussian)

      ActivateSheet(1e14, 1e7); // first visc, then young
      SetSinAnisotropy(10, 1.5);
      // ChangeRelaxThreshold(0.05);
      Only_Extension(true);

      break;
    case 20: // rifting with viscous sheet
      cout << "let's rift" << endl;
      cout << "Latte 2007" << endl;
      Activate_Lattice();
      // Initialize_Fluid_Lattice(1000000, 1000, 0, 4000.0);
      // SetVariationBreakingThreshold(0.7);
      // ReleaseBoundaryParticlesY(1 , 1);

      ////SetGaussianStrengthDistribution(1.0,0.2);

      // WeakenAll(0.2,1,2);

      // ActivateSheet(1*1e20,2*1e9);

      // Initialize_Heat(48, 52);

      // Adjust_break();

      // Loose_Unodes( 0.3,  0.7);

      // SetFracturePlot(20,1);

      ////Only_Extension(true);

      // ChangeRelaxThreshold(2);

      SetVariationBreakingThreshold(0.5); // Variation on breaking strings (<1),
                                          // before (0.7, 0.5), 0.3 probieren
      // SetGaussianSpringDistribution (1.0, 0.4);
      Initialize_Fluid_Lattice(1000000, 1000, 0,
                               400.0); // before (1000000, 1000, 0, 4000.0)

      ReleaseBoundaryParticlesY(1, 1); // 1st: break? 2nd: elastic wall

      // SetWallBoundaries (0, 1);

      WeakenAll(0.05, 1.5, 3.0); // WeakenAll(0.2,1,2); auf 0.1
      InitHeat(200);             // 200
      // WeakenHorizontal Slot(49, 51, 30, 70, 1.1, 1, 1) //ymin, ymax, xmin,
      // xmax, constant, viscosity, break strength
      WeakenHorizontalBox(0, 1, 0, 0.4, 1.0, 1.0, 1000, 1);
      WeakenHorizontalBox(0, 1, 0.6, 1.0, 1.0, 1.0, 1000, 1);

      // ActivateSheet(1*1e20,2*1e9); //das hier

      Initialize_Heat_Box(49, 50, 0, 52, 1200);
      Cooling(1200);
      Initialize_Heat_Box(51, 52, 50, 100, 1200);
      Cooling(1200);

      /*
                              Initialize_Heat_Box(48.5, 51.5, 0, 48, 1000);
                              Cooling(32);
                              Initialize_Heat_Box(48.5, 51.5, 52, 100, 1000);
                              Cooling(32);
                              Initialize_Heat_Box(48, 52, 48, 52, 1200);
                              Cooling(35);
                              //different temp? Adjust ref in Cooling

                              Initialize_Heat_Box(42, 45, 42, 45, 1200);
                              Cooling(32);
                              Initialize_Heat_Box(55, 58, 55, 58, 1200);
                              Cooling(32);


                              Initialize_Heat_Box(49, 51, 0, 10, 1200);
                              Cooling(32);
                              Initialize_Heat_Box(49, 51, 9, 20, 1200);
                              Cooling(32);
                              Initialize_Heat_Box(49, 51, 19, 30, 1200);
                              Cooling(32);
                              Initialize_Heat_Box(49, 51, 29, 40, 1200);
                              Cooling(32);
                              Initialize_Heat_Box(49, 51, 39, 50, 1200);
                              Cooling(32);
                              Initialize_Heat_Box(49, 51, 50, 61, 1200);
                              Cooling(32);
                              Initialize_Heat_Box(49, 51, 60, 71, 1200);
                              Cooling(32);
                              Initialize_Heat_Box(49, 51, 70, 81, 1200);
                              Cooling(32);
                              Initialize_Heat_Box(49, 51, 80, 91, 1200);
                              Cooling(32);
                              Initialize_Heat_Box(49, 51, 90, 100, 1200);
                              Cooling(32);
      */
      // Adjust_break();

      Loose_Unodes(0.2, 0.8);

      SetFracturePlot(20, 1);

      ChangeRelaxThreshold(10);

      break;
    case 21: // forest fire
      cout << "forest fire" << endl;
      /* for (i = 0; i < numParticles; i++)
      {


              runParticle->fluid_particles = 0;

        runParticle = runParticle->nextP;
      }
              Activate_Lattice();
              GrowForest(0.01);*/
      break;

    case 22: // mud crack snakes

      cout << "snakes!" << endl;

      Activate_Lattice(); // construct the lattice

      WeakenAll(0.1, 1, 1); // Weaken Youngs modulus of all grains
                            // (young,viscosity (unwichtig) ,breaking strength)

      SetGaussianStrengthDistribution(
          0.9, 0.15); //(univariante Gauss Verteilung, (mean, sigma)

      //	ActivateSheet(1e18,1e9);  // Aktiviere das Viskose Sheet
      //(viskosity, Young)

      //	SetCosViscSheetMod(1.1, 0.5);  // Cosinus Funktion auf sheet
      //(width, height) was 1.1,0.2

      // SetFracturePlot(20,1); // plot fractures after 20 bonds are broken

      break;

    case 23: // Poroelastic Deformation (Hydrofracturing)
      cout << endl << "Poroelastic Deformation (Hydrofracturing)" << endl;
      cout << "Latte 2011" << endl;

      Activate_Lattice();

      // Set_Absolute_Box(100,1);

      // SetGaussianStrengthDistribution(1.0,0.3);	// gives breaking
      // strength to springs b/w particles in full model box

      // SetPhase(1.0,0.0,1.0,1.5);

      Initialize_Pressure_Lattice(
          0,
          4000); // background pressure (0 is by default), vertical depth,
                 // particle radius(0.005 for low res. and 0.0025 for high res.)

      // Melt_Parameters(1.0e1, 30.0e-10); //viscosity (Pa s), compressibility
      // ()

      // WeakenAll(0.75,1.0,1.0);
      // // proceeds to modulus of particles in respective grain list
      //  WeakenAll(0.01,1.0,1.0);
      // Initialize_Pressure_Lattice(0, 200, 0.005);	// background pressure
      // (0 is by default), vertical depth, particle radius(0.005 for low res.
      // and 0.0025 for high res.) Pressure_init_Boundary_cond(0, 3000000, 200);
      // // 1 is flag for one boundary as continuous, ver. depth in meters

      // Make_Seal_Slot(0.15, 0.2, 0.0, 1.0, 1.2);	// reduce layer porosity
      // by increasing particle radius Make_Seal_Slot(0.25, 0.3, 0.0, 1.0, 1.2);
      // Make_Seal_Slot(0.35, 0.4, 0.0, 1.0, 1.2);
      // Make_Seal_Slot(0.45, 0.5, 0.0, 1.0, 1.2);
      // Make_Seal_Slot(0.55, 0.6, 0.0, 1.0, 1.2);
      // Make_Seal_Slot(0.65, 0.7, 0.0, 1.0, 1.2);
      // Make_Seal_Slot(0.75, 0.8, 0.0, 1.0, 1.2);

      // WeakenHorizontalParticleLayer(0.15, 0.2, 8.0, 1.0, 1.0);
      // WeakenHorizontalParticleLayer(0.25, 0.3, 8.0, 1.0, 1.0);
      // WeakenHorizontalParticleLayer(0.35, 0.4, 8.0, 1.0, 1.0);
      // WeakenHorizontalParticleLayer(0.45, 0.5, 8.0, 1.0, 1.0);
      // WeakenHorizontalParticleLayer(0.55, 0.6, 8.0, 1.0, 1.0);
      // WeakenHorizontalParticleLayer(0.65, 0.7, 8.0, 1.0, 1.0);
      // WeakenHorizontalParticleLayer(0.75, 0.8, 8.0, 1.0, 1.0);

      // WeakenHorizontalParticleLayer(0.25, 0.3, 1.2, 1.0, 3.0,0.8);
      // WeakenHorizontalParticleLayer(0.3, 0.4, 0.8, 1.0, 1.0,15);

      // WeakenHorizontalParticleLayer(0.4, 0.52, 1.8, 1.0, 1.8,1.0);
      /*
      WeakenHorizontalParticleLayer(0.52, 0.53, 0.5, 1.0, 1.0,20);
      WeakenHorizontalParticleLayer(0.55, 0.65, 1.0, 1.0, 1.0,0.9);

      WeakenHorizontalParticleLayer(0.65, 0.665, 0.8, 1.0, 1.0,5);
      WeakenHorizontalParticleLayer(0.68, 0.73, 1.0, 1.0, 1.0,0.9);
      WeakenHorizontalParticleLayer(0.73, 0.735, 0.5, 1.0, 1.0,20);
      WeakenHorizontalParticleLayer(0.735, 0.80, 1.0, 1.0, 0.8,0.85);
      WeakenHorizontalParticleLayer(0.80, 0.81, 0.5, 1.0, 1.0,20);
      WeakenHorizontalParticleLayer(0.815, 0.84, 1.0, 1.0, 0.8,0.85);
      WeakenHorizontalParticleLayer(0.84, 0.85, 0.5, 1.0, 1.0,20);
      WeakenHorizontalParticleLayer(0.855, 0.95, 1.0, 1.0, 0.8,0.85);
      WeakenHorizontalParticleLayer(0.0, 0.25, 1.0, 1.0, 3.0,0.85);

      //WeakenHorizontalParticleLayer(0.0, 0.1, 0.2, 1.0, 1.0);

      WeakenTiltedParticleLayerabsolute (-3,-0.2, 1.2,0.8, 1.0,0.006, 0.002); //
      for res 100 0.006 radius, for 200 0.003 permeable fault, weak

      WeakenTiltedParticleLayerabsolute (-3,-0.25, 1.2,0.8, 1.0,0.006, 0.0025);
      // for res 100 0.006 radius, for 200 0.003




      WeakenHorizontalParticleLayer(0.95, 1.0, 0.1, 1.0, 1.0,20);

      WeakenTiltedParticleLayerabsolute (1.05,1.55, -0.2,0.07, 1.0,0.004,
      0.00125); // this is the weak top
*/

      // WeakenHorizontalParticleLayer(0.3, 0.4, 0.5, 1.0, 1.0,5.8);
      // WeakenHorizontalParticleLayer(0.5, 0.6, 0.5, 1.0, 1.0,5.8);
      // WeakenHorizontalParticleLayer(0.7, 0.8, 0.5, 1.0, 1.0,5.8);
      // WeakenHorizontalParticleLayer(0.1, 0.2, 0.5, 1.0, 1.0,5.8);
      // WeakenHorizontalParticleLayer(0.3, 0.6, 0.5, 1.0, 1.0,0.8);

      SetVariationBreakingThreshold(0.6);
      WeakenAll(0.6, 1.0, 2.0);
      // SetDistributionPorosity(5);

      /**************************************************************
       * Particles are scaled to real radius for Gravity
       * last number is width of whole simulation box in meters
       * ***********************************************************/

      Real_Density_Young(0.45, 0.55, 2300, 2500, 1.0e+10, 5.0e+10, 600);

      // ChangeRelaxThreshold(0.05);
      /*******************************
       * Fluid Pressure Part
       * ******************************/

      Particle_Fluid_parameters();

      Pressure_initialize_Hydrostatic_gradient(
          600, 1000); // box size, depth in meters

      Calculate_Pressure_1st_time(1, 1e-09);

      Relaxation(); // and do a full relaxation

      cout << "applied hydrostatic gradient " << endl;
      /***********************************
       * Apply Gravity
       ***********************************/

      Gravity(1000, 2700); // depth, density
      Relaxation();
      cout << "applied gravity " << endl;

      // ChangeRelaxThreshold(0.5);
      // Relaxation();
      cout << "applied gravity 0.5" << endl;
      // ChangeRelaxThreshold(0.05);
      // cout << "applied gravity 0.05" << endl;

      // ChangeRelaxThreshold(0.05);
      // Relaxation();

      // UpdateElle ();

      UpdateElle(); // and update the interface of Elle

      SetFracturePlot(10, 1);

      break;
    case 24: // Poroelastic Deformation (Hydrofracturing)
      cout << endl << "Mobile Hydros" << endl;
      cout << "Latte 2013" << endl;

      Activate_Lattice();

      Initialize_Pressure_Lattice(
          10,
          1000); // background pressure (0 is by default), vertical depth,
                 // particle radius(0.005 for low res. and 0.0025 for high res.)
      SetDistributionPorosity(4.0);

      //	Initialize_Concentration_Lattice(5.0);

      // WeakenHorizontalParticleLayer(0.3, 0.4, 0.5, 1.0, 1.0,5.8);
      // WeakenHorizontalParticleLayer(0.5, 0.6, 0.5, 1.0, 1.0,5.8);
      // WeakenHorizontalParticleLayer(0.7, 0.8, 0.5, 1.0, 1.0,5.8);
      // WeakenHorizontalParticleLayer(0.1, 0.2, 0.5, 1.0, 1.0,5.8);
      // WeakenHorizontalParticleLayer(0.0, 1.0, 1.0, 1.0, 1.0,1.17);
      // WeakenHorizontalParticleLayer(0.0, 0.5, 0.0, 1.0,  1.0, 1.0, 1.0,1.13);
      // WeakenHorizontalParticleLayer(0.5, 0.7, 0.0, 1.0, 1.0, 1.0, 1.0,1.17);

      // WeakenHorizontalParticleLayer(0.3, 0.32, 0.9, 1.0, 1.0,1.1);

      SetVariationBreakingThreshold(0.6);
      // WeakenAll(0.5,1.0,1.0);

      /**************************************************************
       * Particles are scaled to real radius for Gravity
       * last number is width of whole simulation box in meters
       * ***********************************************************/

      Real_Density_Young(0.45, 0.55, 2300, 2500, 1.0e+10, 5.0e+10, 1000);

      // ChangeRelaxThreshold(0.05);
      /*******************************
       * Fluid Pressure Part
       * ******************************/

      Particle_Fluid_parameters();

      Pressure_initialize_Hydrostatic_gradient(1000,
                                               10); // box size, depth in meters

      // Calculate_Pressure_1st_time(1, 1e-09);

      // Calculate_Concentration();
      Only_Extension(true);

      Relaxation(); // and do a full relaxation

      cout << "applied hydrostatic gradient " << endl;
      /***********************************
       * Apply Gravity
       ***********************************/

      Gravity(1000, 2700); // depth, density
      Relaxation();
      cout << "applied gravity " << endl;

      // ChangeRelaxThreshold(0.1);
      // Relaxation();
      // cout << "applied gravity 0.5" << endl;
      // ChangeRelaxThreshold(0.05);
      // cout << "applied gravity 0.05" << endl;

      // ChangeRelaxThreshold(0.05);
      // Relaxation();

      // UpdateElle ();

      UpdateElle(); // and update the interface of Elle

      // SetFracturePlot(10,1);

      break;

    case 25: // Poroelastic Deformation (Hydrofracturing)
      cout << endl << "diffusion-advection" << endl;
      cout << "Latte 2014" << endl;

      Activate_Lattice();

      Initialize_Pressure_Lattice(
          0,
          1000); // background pressure (0 is by default), vertical depth,
                 // particle radius(0.005 for low res. and 0.0025 for high res.)
      SetDistributionPorosity(5.0);

      Initialize_Concentration_Lattice(0.0);

      // WeakenHorizontalParticleLayer(0.3, 0.4, 0.5, 1.0, 1.0,5.8);
      // WeakenHorizontalParticleLayer(0.5, 0.6, 0.5, 1.0, 1.0,5.8);
      // WeakenHorizontalParticleLayer(0.7, 0.8, 0.5, 1.0, 1.0,5.8);
      // WeakenHorizontalParticleLayer(0.1, 0.2, 0.5, 1.0, 1.0,5.8);
      // WeakenHorizontalParticleLayer(0.0, 1.0, 1.0, 1.0, 1.0,1.17);
      // WeakenHorizontalParticleLayer(0.0, 1.0, 1.0, 1.0, 1.0,1.08);

      // WeakenHorizontalParticleLayer(0.3, 0.32, 0.5, 1.0, 1.0,0.5);

      // WeakenHorizontalParticleLayer(0.0, 1.0, 1.0, 1.0, 1.0,1.17);

      SetVariationBreakingThreshold(0.6);
      WeakenAll(0.6, 1.0, 80.0);

      /**************************************************************
       * Particles are scaled to real radius for Gravity
       * last number is width of whole simulation box in meters
       * ***********************************************************/

      Real_Density_Young(0.45, 0.55, 2300, 2500, 1.0e+10, 5.0e+10, 200);

      // ChangeRelaxThreshold(0.05);
      /*******************************
       * Fluid Pressure Part
       * ******************************/

      Particle_Fluid_parameters();
      // Pressure_initialize_Hydrostatic(200, 200, 0); // boxsize, depth, add
      // pressure
      Pressure_initialize_Hydrostatic_gradient(
          200, 200); // box size, depth in meters

      Calculate_Pressure_1st_time(1, 1e-09);

      Calculate_Concentration(100000, 0.00000025);
      Interpolation_Concentration();

      // Relaxation();		// and do a full relaxation

      cout << "applied hydrostatic gradient " << endl;
      /***********************************
       * Apply Gravity
       ***********************************/

      // Gravity(1000,2700); //depth, density
      Relaxation();
      cout << "applied gravity " << endl;

      // ChangeRelaxThreshold(0.1);
      // Relaxation();
      // cout << "applied gravity 0.5" << endl;
      // ChangeRelaxThreshold(0.05);
      // cout << "applied gravity 0.05" << endl;

      Set_Calcite();

      // ChangeRelaxThreshold(0.05);
      // Relaxation();

      // UpdateElle ();

      UpdateElle(); // and update the interface of Elle

      // SetFracturePlot(10,1);

      break;

    case 26: // Poroelastic Deformation (Hydrofracturing)
      cout << endl << "new fluid" << endl;
      cout << "Latte 2016" << endl;

      Activate_Lattice();

      // Loose_Springs(true);
      // Shift_Particles(0.1);
      // Shift_Particles(0.1);
      // Shift_Particles(0.1);
      // Adjust_Springs_Particles();

      Initialize_Fluid_Lattice(
          10, 1, 1, 1.0); // first background pressure, scale in m,
                          // concentration, time factor (background sec)

      SetDistributionPorosity(2.0);

      Fluid_Parameters(0, 1, 0.000001);

      SetVariationBreakingThreshold(0.4); // was 0.2
      WeakenAll(1.0, 1.0, 0.8);

      Only_Extension(true);

      Gravity(100, 2700); // depth, density
      ChangeRelaxThreshold(0.1);

      Relaxation();

      // WeakenHorizontalParticleLayer(0.6, 0.7, 0.0, 0.5, 1.0, 1.0, 1.0,1.17);
      // WeakenHorizontalParticleLayer(0.65, 0.75,
      // 0.5, 1.0, 1.0, 1.0, 1.0,1.17); WeakenHorizontalParticleLayer(0.3, 0.4,
      // 0.0, 0.5, 1.0, 1.0, 1.0,1.17); WeakenHorizontalParticleLayer(0.35,
      // 0.45, 0.5, 1.0, 1.0, 1.0, 1.0,1.17);

      // Hydrostatic_Fluid(1000);

      // Set_Calcite();

      cout << "applied gravity " << endl;
      // UpdateElle ();	// and update the interface of Elle

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

      Hydrostatic_Fluid(1000);

      Set_Calcite();

      cout << "applied gravity " << endl;
      UpdateElle(); // and update the interface of Elle

      break;
    case 28: // advection-diffusion-reaction
      cout << endl << "new fluid react" << endl;
      cout << "Latte 2018" << endl;

      Activate_Lattice();

      // this makes plus minus random lattices, but then the elasticity is
      // strange
      // Loose_Springs(true);
      // Shift_Particles(0.1);
      // Shift_Particles(0.1);
      // Shift_Particles(0.1);
      // Adjust_Springs_Particles();

      //---------------------------------------------------------------------------------
      // this function Initialized the lattice
      //
      // first number is background fluid pressure, second number is scale (x)
      // in m third number is initial concentration of all nodes and last number
      // is time in seconds for the pressure diffusion.
      //---------------------------------------------------------------------------------

      time_step = 1.0; // in seconds
      size = 1;        // in meters
      int iii;

      // Initialize_Pressure_Lattice(0,1000);

      Initialize_Fluid_Lattice(10, size, 0.0,
                               time_step); // setting up the fluid

      // this would set an initial fluid pressure gradient
      // Hydrostatic_Fluid(10);

      //---------------------------------------------------------------------------------
      // change porosity of the matrix. first value ymin, then yax, then xmin
      // and x max then Youngs modulus, viscosity, breaking strength and fluid
      // radius
      //
      // fluid radius larger than 1.0 makes it less permeable
      //---------------------------------------------------------------------------------

      // WeakenHorizontalParticleLayer(0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0,1.06);
      // // was 1.08

      // WeakenHorizontalParticleLayer(0.0, 0.1, 0.48, 0.52, 1.0, 1.0, 1.0,0.7);

      // WeakenTiltedParticleLayerabsolute(0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.04);

      // this sets a random distribution on the porosity

      SetDistributionPorosity(0.3);

      // setting some vertical and horizontal fractures

      // WeakenTiltedParticleLayerabsolute(0.18, 0.22, 0.0, 1.0, 1.0, 1.0,
      // 0.0025); WeakenTiltedParticleLayerabsolute(0.38, 0.42,
      // 0.0, 1.0, 1.0, 1.0, 0.0025); WeakenTiltedParticleLayerabsolute(0.58,
      // 0.62, 0.0, 1.0, 1.0, 1.0, 0.0025);
      // WeakenTiltedParticleLayerabsolute(0.78, 0.82, 0.0, 1.0, 1.0, 1.0,
      // 0.0025);

      // WeakenTiltedParticleLayerabsolute(0.3, 0.35, 0.5, 1.0, 1.0, 1.0,
      // 0.0025);

      // WeakenHorizontalParticleLayer(0.0, 1.0, 0.1, 0.11, 1.0, 1.0, 1.0,0.9);
      // //first y then x WeakenHorizontalParticleLayer(0.0, 1.0, 0.2,
      // 0.21, 1.0, 1.0, 1.0,0.9); //first y then x
      // WeakenHorizontalParticleLayer(0.0, 1.0, 0.3, 0.31, 1.0, 1.0, 1.0,0.9);
      // //first y then x WeakenHorizontalParticleLayer(0.0, 1.0, 0.4,
      // 0.41, 1.0, 1.0, 1.0,0.9); //first y then x
      // WeakenHorizontalParticleLayer(0.0, 1.0, 0.6, 0.61, 1.0, 1.0, 1.0,0.9);
      // //first y then x WeakenHorizontalParticleLayer(0.0, 1.0, 0.8,
      // 0.81, 1.0, 1.0, 1.0,0.9); //first y then x

      // WeakenHorizontalParticleLayer(0.3, 0.31, 0.0, 1.0, 1.0, 1.0, 1.0,0.9);
      // //first y then x WeakenHorizontalParticleLayer(0.6, 0.61,
      // 0.0, 1.0, 1.0, 1.0, 1.0,0.9); //first y then x

      // variation of some elastic parameters and breaking strength

      // variation of breaking threshold

      SetVariationBreakingThreshold(0.6);

      // change overall parameters, first elastic constant, second viscosity,
      // third breaking strength
      WeakenAll(0.8, 1.0, 1.0);

      // angular forces or not
      Only_Extension(true);

      // make the grain boundary permeable
      // Permeable_Grain_Boundary(0.06); // was 0.01

      /*Permeable_Grain(0.05,65);
      Permeable_Grain(0.05,61);
      Permeable_Grain(0.05,43);
      Permeable_Grain(0.05,68);
      Permeable_Grain(0.05,59);
      Permeable_Grain(0.05,47);
      Permeable_Grain(0.05,17);
      Permeable_Grain(0.05,3);
      Permeable_Grain(0.05,5);
      Permeable_Grain(0.05,14);
      Permeable_Grain(0.05,18);
      Permeable_Grain(0.05,39);
      Permeable_Grain(0.05,42);
      Permeable_Grain(-0.01,63);
      Permeable_Grain(-0.01,67);
      Permeable_Grain(0.05,66);
      Permeable_Grain(0.05,62);
      Permeable_Grain(0.05,48);


      //Permebale Layer
      Permeable_Grain(0.05,31);
      Permeable_Grain(0.05,38);
      Permeable_Grain(0.05,25);
      Permeable_Grain(0.05,19);
      Permeable_Grain(0.05,28);
      Permeable_Grain(0.05,40);
      Permeable_Grain(0.05,41);
      Permeable_Grain(0.05,29);
      Permeable_Grain(0.05,30);
      Permeable_Grain(0.05,37);
      Permeable_Grain(0.05,34);
      Permeable_Grain(0.05,27);
      Permeable_Grain(0.05,20);
      Permeable_Grain(0.05,21);
      Permeable_Grain(0.05,32);

      Permeable_Grain(0.05,26);
      Permeable_Grain(0.05,36);
      Permeable_Grain(0.05,22);
      Permeable_Grain(0.05,23);
      Permeable_Grain(0.05,24);
      Permeable_Grain(0.05,35);
      Permeable_Grain(0.05,46);

      //Impermeable layer
      Permeable_Grain(-0.01,44);
      Permeable_Grain(-0.01,45);
      Permeable_Grain(-0.01,7);
      Permeable_Grain(-0.01,9);
      Permeable_Grain(-0.01,11);
      Permeable_Grain(-0.01,0);
      Permeable_Grain(-0.01,15);
      Permeable_Grain(-0.01,12);
      Permeable_Grain(-0.01,1);
      Permeable_Grain(-0.01,8);
      Permeable_Grain(-0.01,4);
      Permeable_Grain(-0.01,13);
      Permeable_Grain(-0.01,6);
      Permeable_Grain(-0.01,16);
      Permeable_Grain(-0.01,2);

      Permeable_Grain(-0.01,50);
      Permeable_Grain(-0.01,51);
      Permeable_Grain(-0.01,56);
      Permeable_Grain(-0.01,57);
      Permeable_Grain(-0.01,73);
      Permeable_Grain(-0.01,49);
      Permeable_Grain(-0.01,58);
      Permeable_Grain(-0.01,74);
      Permeable_Grain(-0.01,55);
      Permeable_Grain(-0.01,69);
      Permeable_Grain(-0.01,72);
      Permeable_Grain(-0.01,60);
      Permeable_Grain(-0.01,53);
      Permeable_Grain(-0.01,71);
      Permeable_Grain(-0.01,70);
      Permeable_Grain(-0.01,64);
      Permeable_Grain(-0.01,54);*/

      /*Permeable_Grain(0.06,5);
      Permeable_Grain(0.06,6);
      Permeable_Grain(0.06,4);
      Permeable_Grain(0.06,3);
      Permeable_Grain(0.06,2);
      Permeable_Grain(0.06,21);
      Permeable_Grain(0.06,23);
      Permeable_Grain(0.06,22);
      Permeable_Grain(0.06,24);*/

      // give fluid parameters, first number boundary condition, second number
      // cozeny grain size
      Fluid_Parameters(0.0000000001, 1, 0.00001);

      // cout << "applied gravity " << endl;
      // UpdateElle ();	// and update the interface of Elle

      /*
      break;
      // give fluid parameters, first number boundary condition, second number
      cozeny grain size Fluid_Parameters(1, 0.00001);



      for (iii = 0; iii<100; iii++) // for max and min x value

      {
              for (jjj = 28; jjj < 68; jjj++) // for y coordinates, how thick is
      permable layer?

              {
                      Add_Concentration(iii, jjj, 0.5);

              }

      }

*/

      Relaxation();

      Set_Calcite();

      // cout << "applied gravity " << endl;
      // UpdateElle ();	// and update the interface of Elle

      break;

    case 29:
      cout << "Andrew rifting and dykes" << endl;

      Activate_Lattice(); // constructor

      Set_Absolute_Box(
          1000, 0); // model size in meters, sets real radius of particles

      SetVariationBreakingThreshold(0.3); // variation on breaking threshold

      WeakenAll(3.5, 1.0, 0.15); // make everything different (young, viscosity,
                                 // breaking threshold)

      UnsetNobreakUpperRow(); // release upper row

      InitRealDensityYoung(
          2300, 0.4e11); // density of all plus Young, set for single particles

      MarkMohoParticles(
          0.825); // mark a row of particles for output, good for plots

      WeakenHorizontalParticleLayer(0.0, 0.3, 0.0, 1.0, 1.0, 1.0, 80.0, 1.0);
      // WeakenHorizontalParticleLayer(0.2, 0.3, 0.6, 1.0, 4.0 ,1.0);
      // WeakenHorizontalParticleLayer(0.3, 0.4, 0.8, 1.0, 2.0 ,1.0);

      // WeakenHorizontalParticleLayer(0.9, 1.0, 1.0, 1.0, 0.5 ,1.0);

      // InitInternalYoung(1.0e11);  // converts real Young back to
      // non-dimensional one

      InitViscosity(1e23); // sets background viscosity

      // visosity setting for layers or blocks

      // SetViscosity(1.0,0.0,1.0,0.80,1e25);  //oberste Kruste
      // SetViscosity(0.6,0.4,0.4, 0.0,1e20);  //unterer Teil oberste Kruste

      // and set breaking strength for layers or blocks

      // SetBreakingStrength(0.49,0.51,0.5,0.45,0.4);
      // SetBreakingStrength(1.0,0.0,0.98,0.96,0.8);

      cout << "in adjust" << endl;

      AdjustParticleConstants();

      cout << "out adjust" << endl;

      Gravity_noload_nofluid();

      // ChangeRelaxThreshold(10);

      SetFracturePlot(10, 1);

      // Only_Extension(true);

      break;

    case 30:
      cout << "Andrew rifting and dykes with dykes!" << endl;

      Activate_Lattice(); // constructor

      Set_Absolute_Box(
          1, 0); // model size in meters, sets real radius of particles

      Initialize_Fluid_Lattice(1, 1000, 0.0, 2.0);

      SetDistributionPorosity(4.0);

      Fluid_Parameters(0, 1, 0.0000001);

      SetVariationBreakingThreshold(0.3); // variation on breaking threshold

      WeakenAll(1.0, 1.0, 0.15); // make everything different (young, viscosity,
                                 // breaking threshold)

      UnsetNobreakUpperRow(); // release upper row

      InitRealDensityYoung(
          2300, 0.4e11); // density of all plus Young, set for single particles

      MarkMohoParticles(
          0.825); // mark a row of particles for output, good for plots

      // WeakenHorizontalParticleLayer(0.0, 0.05, 0.3, 1.0, 50.0 ,1.0);
      // WeakenHorizontalParticleLayer(0.05, 0.2, 0.8, 1.0, 5.0 ,1.0);
      // WeakenHorizontalParticleLayer(0.3, 0.4, 0.8, 1.0, 2.0 ,1.0);

      // WeakenHorizontalParticleLayer(0.9, 1.0, 1.0, 1.0, 0.5 ,1.0);

      // InitInternalYoung(1.0e11);  // converts real Young back to
      // non-dimensional one

      InitViscosity(1e23); // sets background viscosity

      // visosity setting for layers or blocks

      // SetViscosity(1.0,0.0,1.0,0.80,1e25);  //oberste Kruste
      // SetViscosity(0.2,0.0,1.0, 0.0,1e20);  //unterer Teil oberste Kruste

      // and set breaking strength for layers or blocks

      // SetBreakingStrength(0.49,0.51,0.5,0.45,0.4);
      // SetBreakingStrength(1.0,0.0,0.98,0.96,0.8);

      cout << "in adjust" << endl;

      AdjustParticleConstants();

      cout << "out adjust" << endl;

      Gravity_noload_nofluid();

      // ChangeRelaxThreshold(10);

      SetFracturePlot(10, 1);

      // Only_Extension(true);

      break;

    case 31: // advection-diffusion-reaction  // large scale fluid movement
      cout << endl << "new fluid react large scale" << endl;
      cout << "Latte 2020" << endl;

      Activate_Lattice();

      //---------------------------------------------------------------------------------
      // this function Initialized the lattice
      //
      // first number is background fluid pressure, second number is scale (x)
      // in m third number is initial concentration of all nodes and last number
      // is time factor for the pressure diffusion.
      //---------------------------------------------------------------------------------

      time_step = 5.0; // in seconds
      size = 100.0;    // in meters

      Initialize_Fluid_Lattice(10, size, 0.0,
                               time_step); // setting up the fluid

      // this would set an initial fluid pressure gradient
      // Hydrostatic_Fluid(10);

      //---------------------------------------------------------------------------------
      // change porosity of the matrix. first value ymin, then yax, then xmin
      // and x max then Youngs modulus, viscosity, breaking strength and fluid
      // radius
      //
      // fluid radius larger than 1.0 makes it less permeable
      //---------------------------------------------------------------------------------

      SetDistributionPorosity(1.0);

      // variation of some elastic parameters and breaking strength

      // variation of breaking threshold
      SetVariationBreakingThreshold(0.6);

      // change overall parameters, first elastic constant, second viscosity,
      // third breaking strength
      WeakenAll(0.8, 1.0, 1.0);

      // angular forces or not
      Only_Extension(true);

      // make the grain boundary permeable
      // Permeable_Grain_Boundary(0.3);

      // give fluid parameters, first number boundary condition, second number
      // cozeny grain size
      Fluid_Parameters(0, 1, 0.0001);

      Relaxation();

      Set_Calcite();

      cout << "applied gravity " << endl;
      // UpdateElle ();	// and update the interface of Elle

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

      time_step = 1.0; // in seconds
      size = 5.0;      // in meters

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
        // DeformLatticeExtend(0.0001, 0);
        DeformLatticePureShear(0.0001, 0);
      }
      // DeformLatticePureShear(0.0001, 0);
      Get_Aperture();
      // else
      // DeformLatticePureShear(0.001,1);

      // DeformLatticePureShearAreaChange (0.00001, 0.00001, 1);
      // ViscousRelax(1, 1e15);
      UpdateElle();
      break;

    case 2: // fracture boudinage
      DeformLatticePureShear(0.001, 1);
      break;

    case 3: // expanding inclusions

      // ShrinkGrain(50, 0.002, 1);
      ShrinkGrain(28, 0.00008, 0);
      ShrinkGrain(15, 0.00008, 0);
      // ShrinkGrain(69, 0.001, 1);
      ShrinkGrain(47, 0.00008, 0);

      Relaxation();

      break;

    case 4: // shrinkage patterns
      ShrinkBox(0.001, 1, 0);
      // RelaxSheet(1.5e8);
      break;

    case 5: // viscos relax
      DeformLatticePureShear(0.001, 1);
      ViscousRelax(1, 1e11); // strain 0.001 and 1e11 are strain rate 10^-12
      break;

    case 6: // grooves
      DeformLattice(0.002, 1);
      Dissolution_Strain(20);
      break;

    case 7: // Stylolites
      //-------------------------------------------------
      // Deform the lattice from upper and lower
      // boundaries. Upper and lower part of the lattice
      // are now just pressed together assuming
      // there is no resistance. Stresses build up once
      // the two sides meet. Deformation steps have to
      // be very small
      // The side walls are fixed. Movement is vertical
      // in steps of 0.00005 * y size (1.0)
      // second number means make a picture after the
      // movement. 0 means take no picture.
      //-------------------------------------------------
      DeformLatticeNewAverage2side(0.0001, 1);
      //---------------------------------------------------
      // Dissolution routine. Particles are just dissolved
      // depending on stress, elastic and surface
      // energies. One particle is dissolved in a time
      // step. 100 means take a picture after 100 particles
      // have dissolved.
      //---------------------------------------------------
      Dissolution_StylosII(1, 4, 0, 0, 1);

      // Dumplayer (i);

      DumpStatisticSurface(i);

      break;

    case 8: // combined grain growth and fracturing
      // DeformLatticePureShear(0.001,1);  // pure shear deformation strain in y
      // direction is 0.1 %
      DoGrowth(i); // grain growth step, calls growth function in graingrowth.cc
      GetNewElleStructureTwo(); // reread the Elle structure for Latte after
                                // grain growth
      break;

    case 9:                    // phase change
      DeformLattice(0.002, 1); // deform and relax //(0.0002,x):
                               // deformationsstep, ursprünglich 0.002
      heat_distribution.SetHeatFlowEnabled(1);
      //***float Lattice::DeformLattice(float move, int plot)***
      DumpStatisticStressBox(0.2, 0.8, 0.2, 0.8, 0.002);
      // Start_Reactions();
      break;

    case 10: // Rift

      if (i == 0) {
        // ChangeRelaxThreshold(0.1);
        ApplyGravity(200000);
        Relaxation();
        UpdateElle();
        // for (j=0;j<10;j++)
        // ViscousRelax(1,300 * 365*24*60*60);
        // ChangeRelaxThreshold(100);
      }
      // ChangeRelaxThreshold(0.5);
      time_step = 1000; // years , have to use 2000 times steps to get 2 mio
                        // years if time is 1000
      // if (i == 100)
      // SetRealDensityYoung(0.8,0.3,0.72,0.25,3100,1.8e11);
      // if (i == 150)
      // SetRealDensityYoung(0.8,0.3,0.72,0.25,3000,1.8e11);
      // if (i == 200)
      // SetRealDensityYoung(0.85,0.25,0.72,0.25,2900,1.8e11);
      // if (i == 250)
      // SetRealDensityYoung(0.85,0.25,0.72,0.25,2800,1.8e11);

      // if (i < 2000)
      DeformLatticeXMidTwo_b(-SpreadingStrain(0.003, 200000, time_step), 0.8,
                             1);

      // else if (i < 4000)
      // DeformLatticeXMidTwo_b(-SpreadingStrain(0.004,200000,time_step),0.8,1);

      // else
      // DeformLatticeXMidTwo_b(-SpreadingStrain(0.006,200000,time_step),0.8,1);

      // second term, 0.85 is the height of the localized versus distributed
      // stretching in the model below is distributed and above localized.
      // spreading strain calculates move of the right wall boundary
      // move is first number in spreadstrain times years divided by 200 000.
      // that multiplied by 200 000 again gives movement in mm per year. so
      // 0.004 in 1000 years is 4 mm per year

      // if (i < 600)
      {
        // DeformLatticeXMid(-SpreadingStrain(0.002,200000,time_step), 1, 0.8,
        // 0.4, 0.6);
      }
      // else
      {
        // DeformLatticeXMid(-SpreadingStrain(0.003,200000,time_step), 1, 0.8,
        // 0.4, 0.6); //0.004 is 4mm/year Last model dec 2012
      }

      ViscousRelax(1, time_step * 365 * 24 * 60 * 60);
      // DamageModel(0.25,0.025,2);

      OutputTopParticles((i + 1 * time_step));
      OutputMohoParticles((i + 1 * time_step));

      break;

    case 11:       // pure grain growth
      DoGrowth(i); // do only grain growth (in graingrowth.cc)
      break;

    case 12:                   // diffusion
      UpdateFluidLatticeGas(); // run one lattice gas step (transport +
                               // collisions)
      break;

    case 13:                            // fluid flow
      InsertFluidLatticeGas(0.03, 0.5); // pump in fluid at left boundary
      UpdateFluidLatticeGas(); // one lattice gas step (transport + collisions)
      RemoveFluidLatticeGas(0.97); // suck out fluid at right boundary
      break;

    case 15:
      // DeformLattice(0.001,1);
      DeformLatticeNoAverage(0.001, 1);
      GrowthDissolution(20, 1, 1, experiment_time);
      break;

    case 16:
      DeformLattice(0.001, 1);
      GrowthDissolution(1, 1, 0, experiment_time);
      break;

    case 17:
      DeformLattice(0.002, 1); // deform and relax //(0.0002,x):
                               // deformationsstep, ursprünglich 0.002
      heat_distribution.SetHeatFlowEnabled(1);
      //***float Lattice::DeformLattice(float move, int plot)***
      DumpStatisticStressBox(0.2, 0.8, 0.2, 0.8, 0.002);
      // Start_Reactions();
      break;

    case 18:
      DeformLattice(0.002, 1); // deform and relax //(0.0002,x):
                               // deformationsstep, ursprünglich 0.002
      heat_distribution.SetHeatFlowEnabled(1);
      //***float Lattice::DeformLattice(float move, int plot)***
      DumpStatisticStressBox(0.2, 0.8, 0.2, 0.8, 0.002);
      // Start_Reactions();
      break;

    case 19:

      // DeformLatticeExtension(0.0001,1,1.0,-0.0);

      ShrinkBox(0.001, 1, 0);

      RelaxSheet(1e18);

      Relaxation();

      break;

    case 20:

      // DeformLatticeExtension(0.0005,1,1.0,-0.0);

      // solve_Heat(1000, 10000, 0.000001);

      // RelaxSheet(1e20);

      // Relaxation();

      ////Healing(1.1,0.03,experiment_time,1.0,0.3);

      // Find_Empty_Box(experiment_time);

      // Relaxation();

      DeformLatticeExtension(0.001, 1, 1.0,
                             0.0); // slower? alt:0.0005, shear > 1.0?
      // DeformLatticeXMid(0.005, 1, 0.0, 0.4, 0.6);

      // RelaxSheet(1e20);

      // Relaxation();

      Initialize_Heat_Box(49 + int(i / 20), 51 + int(i / 20), 0, 52, 1200);
      Cooling(1200);
      Initialize_Heat_Box(51 + int(i / 20), 53 + int(i / 20), 50, 100, 1200);
      Cooling(1200);

      solve_Heat(1000, 1000, 0.000004); // space, time, dif const; 0.000001

      Find_Empty_Box(experiment_time, 1200);

      Cooling(1200);

      Relaxation();

      Healing(1.1, 0.03, experiment_time, 1.0, 0.3); // vielleicht doch dazu?

      // Relaxation();

      break;

    case 21: // forest fire

      // GrowForest(0.09);
      // case 21:�
      //  ForestFire(0.04);
      // GrowForest(0.01);
      /* for (i = 0; i < numParticles; i++)
            {

                    if (runParticle->fluid_particles == -1)
                            runParticle->fluid_particles = 0;

              runParticle = runParticle->nextP;
            }*/
      // ForestFireStat();
      break;

    case 22: // mud crack snakes

      ShrinkBox(0.001, 1, 0); // shrink the box from all sides (shrink by amount
                              // radius * x, plot=1 noplot=0, 4=channel)

      // RelaxSheet(1.5e8);       // viscous Relaxation step (time step in
      // seconds)

      break;

    case 23:

      /*************************************************************
       * Model to fracture and heal with fluid pressure
       * For plotting:
       * U_ATTRIB_A = healing age
       * U_ATTRIB_B = fracture age
       * U_ATTRIB_C = fluid pressure gradient
       * U_FRACTURES = active fractures
       * U_Temperature = fluid pressure
       * DISLOCATION_DENSITY = average porosity
       * VISCOSITY = crack_seal_counter cycle
       *
       *************************************************************/

      /***********************************************************
       * Calculate Fluid Pressure Gradient in the Run Function
       * **********************************************************/

      Particle_Fluid_parameters(); // valid for all functions

      Pressure_initialize_Hydrostatic(600, 200,
                                      i * 500); // boxsize, depth, add pressure

      // Pressure_xGauss_RandSource_area(0, 0.4, 0, 1.0, 20000000, 12.0);

      for (iii = 0; iii < 100; iii++)
        Calculate_Pressure_1st_time(
            1, 1e-09); // note second term is not used - takes the "time_a" set
                       // in pressure_lattice

      cout << "out" << endl;
      Relaxation();

      /**************************************************************
       * Deform the Lattice
       * ***********************************************************/

      DeformLatticeExtend(0.0000001,
                          0); // strain per step, plot yes (1), no (0)

      UpdateElle();

      // Fracture_Reaction(8,0.2,7.0);

      /****************************************************
       * Apply the Healing Function
       * ************************************************/

      Healing(1.1, 0.03, experiment_time, 1.0, 0.3);
      Relaxation();

      /****************************************************
       * Write out some Statistics
       * ************************************************/

      DumpStatisticGravityStress(50);
      DumpStatisticStressBox(0.1, 0.8, 0.2, 0.8, i);
      DumpStat_Boxes_Stress(10, 25);
      DumpStat_Boxes_Pore_Perm(30, 30);
      DumpStatisticPorosity(0.3, 0.7, 0.3, 0.7, i);

      /**************************
       * Update Interface
       * ************************/

      UpdateElle();

      break;
    case 24:

      /*************************************************************
       * Model to fracture and heal with fluid pressure
       * For plotting:
       * U_ATTRIB_A = healing age
       * U_ATTRIB_B = fracture age
       * U_ATTRIB_C = fluid pressure gradient
       * U_FRACTURES = active fractures
       * U_Temperature = fluid pressure
       * DISLOCATION_DENSITY = average porosity
       * VISCOSITY = crack_seal_counter cycle
       *
       *************************************************************/

      /***********************************************************
       * Calculate Fluid Pressure Gradient in the Run Function
       * **********************************************************/
      //	if (i==2000)
      //	Gravity_no_Move(2700,2700);
      //	if (i==4000)
      //	Gravity_no_Move(2400,2700);
      // if (i==6000)
      //	Gravity_no_Move(2100,2700);

      // if (i==0)
      //	{
      //	Particle_Fluid_parameters();
      //	Pressure_initialize_Hydrostatic(1000, 3000, 2000); // boxsize,
      //depth, add pressure 	for (iii=0;iii<30;iii++)
      //	Calculate_Pressure_1st_time(1, 1e-09);
      //	}
      // else
      //{
      Particle_Fluid_parameters(); // valid for all functions

      //	Pressure_initialize_Hydrostatic(1000, 3000, 2000 + (i*1000)); //
      //boxsize, depth, add pressure

      Pressure_initialize_Hydrostatic(1000, 1000, 20);

      // Pressure_xGauss_RandSource_area(0, 0.4, 0, 1.0, 20000000, 12.0);
      // Pressure_Insert_Random_Node(100, 5000000);
      Pressure_Insert_Random_Nodey(0.01, 0.5, 100, 40000000);

      for (iii = 0; iii < 10; iii++)
        Calculate_Pressure_1st_time(
            1, 1e-09); // note second term is not used - takes the "time_a" set
                       // in pressure_lattice
      //}
      Adjust_Gravity();
      cout << "out" << endl;
      Relaxation();

      /**************************************************************
       * Deform the Lattice
       * ***********************************************************/

      // DeformLatticeExtend(0.000001,0); // strain per step, plot yes (1), no
      // (0)

      // UpdateElle ();

      // Fracture_Reaction(8,0.2,7.0);

      /****************************************************
       * Apply the Healing Function
       * ************************************************/

      // Healing(1.1,0.03,experiment_time,1.0,1.0);
      // Relaxation();

      /****************************************************
       * Write out some Statistics
       * ************************************************/

      // DumpStatisticGravityStress(50);
      DumpStatisticStressBox(0.3, 0.7, 0.3, 0.7, i);
      // DumpStat_Boxes_Stress(10, 25);
      // DumpStat_Boxes_Pore_Perm(30, 30);
      DumpStatisticPorosity(0.3, 0.7, 0.3, 0.7, i);

      /**************************
       * Update Interface
       * ************************/

      UpdateElle();

      break;
    case 25:

      /*************************************************************
       * Model to fracture and heal with fluid pressure
       * For plotting:
       * U_ATTRIB_A = healing age
       * U_ATTRIB_B = fracture age
       * U_ATTRIB_C = fluid pressure gradient
       * U_FRACTURES = active fractures
       * U_Temperature = fluid pressure
       * DISLOCATION_DENSITY = average porosity
       * VISCOSITY = crack_seal_counter cycle
       *
       *************************************************************/

      /***********************************************************
       * Calculate Fluid Pressure Gradient in the Run Function
       * **********************************************************/

      Particle_Fluid_parameters(); // valid for all functions

      Pressure_initialize_Hydrostatic(200, 200,
                                      200); // boxsize, depth, add pressure

      // Pressure_xGauss_RandSource_area(0, 0.4, 0, 1.0, 20000000, 12.0);
      // Pressure_Insert_Random_Node(100, 5000000);
      // Pressure_Insert_Random_Nodey(0.0,0.3,100, 5000000);

      for (iii = 0; iii < 2; iii++)
        Calculate_Pressure_1st_time(
            1, 1e-09); // note second term is not used - takes the "time_a" set
                       // in pressure_lattice

      SetDarcyforAdv();

      for (iii = 0; iii < 5; iii++)
        Calculate_Concentration(100000,
                                0.00000025); // first is velocity of fluid,
                                             // second is diffusion coefficient

      // Get_Concentration_Directly();

      Interpolation_Concentration();

      Make_Dolomite(1, i, 0, 1, 1); // number is reaction rate

      // Adjust_Gravity();
      cout << "out" << endl;
      // Relaxation();

      /**************************************************************
       * Deform the Lattice
       * ***********************************************************/

      // DeformLatticeExtend(0.000001,0); // strain per step, plot yes (1), no
      // (0)

      // Fracture_Reaction(8,0.2,7.0);

      /****************************************************
       * Apply the Healing Function
       * ************************************************/

      // Healing(1.1,0.03,experiment_time,1.0,3.0);
      // Relaxation();

      /****************************************************
       * Write out some Statistics
       * ************************************************/

      // DumpStatisticGravityStress(50);
      DumpStatisticStressBox(0.3, 0.7, 0.3, 0.7, i);
      // DumpStat_Boxes_Stress(10, 25);
      // DumpStat_Boxes_Pore_Perm(30, 30);
      DumpStatisticPorosity(0.3, 0.7, 0.3, 0.7, i);

      /**************************
       * Update Interface
       * ************************/

      UpdateElle();

      break;

    case 26:
      //--------------------------------------------------------
      // A is healing
      //
      // B is fracture age
      //
      //----------------------------------------------------------

      Fluid_Parameters(0, 1, 0.000001);

      // for (iii=0;iii<4;iii++)
      {
        // Fluid_Insert_Random_Nodexy(0.45, 0.55, 0.45, 0.55,  15000000);
        // Fluid_Insert_Random_Nodexy(0.20, 0.40, 0.5, 1.0,  1500000);
      }
      // Fluid_Insert_Random_Nodexy(0.51, 0.53, 0.49, 0.51,  1500000);
      // Fluid_Insert_Random_Nodexy(0.49, 0.51, 0.51, 0.53,  1500000);
      // Fluid_Insert_Random_Nodexy(0.51, 0.53, 0.51, 0.53,  1500000);

      // for (iii=0;iii<15;iii++)

      // Fluid_Insert_Random_Nodexy(0.45, 0.65, 0.5, 1.0,  1500000);

      // for (iii=0;iii<10;iii++)

      //  Fluid_Insert_Random_Nodexy(0.3, 0.5, 0.1, 0.4,  1000000);

      // for (iii=0;iii<10;iii++)

      // Fluid_Insert_Random_Nodexy(0.2, 0.4, 0.4, 0.6,  100000000);

      Input_Fluid_Lattice(200000, 10, 10);

      // for (iii=0;iii<10;iii++)

      // Fluid_Insert_Random_Nodexy(0.1, 0.3, 0.6, 0.9,  1000000);

      // Fluid_Insert_Random_Node(0.0, 0.5, 2000000);

      // Input_Fluid_Lattice(70000, iii, 0);
      // Input_Fluid_Lattice(9000000, 101, 100);
      // Input_Fluid_Lattice(10000000, 100, 101);
      // Input_Fluid_Lattice(14000000, 101, 101);

      // UpdateElle ();

      // Calculate_Fluid_Pressure(9800000, (9800000*2), 1, 1, 1, 0.0); //9800000

      Calculate_Fluid_Pressure(10, 10, 1, 1, 1, 0.0); // 9800000

      Adjust_Gravity();

      Relaxation();

      // Healing(1.1,0.5,experiment_time,1.0,1.0);

      // Calculate_Con(100000, 0.00000025);

      // Make_Dolomite(1, i);

      DumpStatisticStressBox(0.1, 0.6, 0.3, 0.7, i);
      // DumpStat_Boxes_Stress(10, 25);
      // DumpStat_Boxes_Pore_Perm(30, 30);
      DumpStatisticPorosity(0.1, 0.6, 0.3, 0.7, i);

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

    case 28: // Advection diffusion and reaction
      //--------------------------------------------------------
      // A is solid
      // Temperature is fluid pressure
      // B is concentration
      // C is fluid pressure gradient
      // Energy is darcy velocity in y
      // Dislocation density is porosity
      //----------------------------------------------------------

      // set fluid parameters
      Fluid_Parameters(0.00000000001, 1, 0.00001);

      // calculate the fluid pressure
      // first upper boundary, second lower boundary, the itterations,
      // then sides wrapping then use average concentration
      // for (iii = 0; iii < 100; iii++)
      // Input_Fluid_Lattice(300*j, iii, 50);
      // for (iii = 0; iii < 100; iii++)
      //	Input_Fluid_Lattice(300*j, 50, iii);
      // for (iii = 0; iii < 50; iii++)
      //	Input_Fluid_Lattice(300*j, 25, iii);
      // for (iii = 0; iii < 50; iii++)
      // Input_Fluid_Lattice(300*j, iii, 25);

      // Calculate the reaction

      cout << "calculate fluid pressure" << endl;

      nonlinear_time =
          Calculate_Fluid_Pressure_timeadjust(10, 10, 1, 1, 0.5, 1000 * j);

      // for (iii=0; iii<5; iii++)
      {
        // Calculate advection, then diffusion
        // if (j<400000)
        cout << "calculate advection" << endl;

        Calculate_Advection(nonlinear_time * 100, size, 2, 1);
      }

      cout << "calculate diffusion" << endl;

      Calculate_Diffusion(0.00001, nonlinear_time * 800, size); // was 5000000

      // Calculate_Diffusion_Metal(0.00001, nonlinear_time*800, size);

      // Calculate_Diffusion_Temperature(0.00001, nonlinear_time, size);

      // Calculate the reaction

      // Make_Metal(nonlinear_time*100, i, 0.0000, -0.0001, size);
      // Make_Dolomite(time_step*50, i, 0.0000, -0.0001, size);

      cout << "calculate reaction" << endl;

      Reaction_Arrhenius(nonlinear_time, size);

      Relaxation();

      // Dump_concentration_growth (1);

      UpdateElle();

      break;

    case 29: // Rift and dykes

      if (i == 0) {
        ChangeRelaxThreshold(0.01);
        ApplyGravity(0);
        Relaxation();
        //	UpdateElle();
        // for (j=0;j<10;j++)
        // ViscousRelax(1,300 * 365*24*60*60);
        // ChangeRelaxThreshold(5);
      }

      // time_step = 100; //years , have to use 2000 times steps to get 2 mio
      // years if time is 1000

      DeformLatticeXMid(-0.0000003, 1, 0.2, 0.3, 0.7);
      DeformLatticeXMid(-0.0000003, 1, 0.4, 0.4, 0.6);
      DeformLatticeXMid(-0.0000004, 1, 0.6, 0.45, 0.55);

      // DeformLatticeExtend (0.0000001, 1);

      Relaxation();

      Save_Particle_Velocity();

      // ViscousRelax(1,time_step*365*24*60*60);

      UpdateElle();

      OutputTopParticles((i + 1 * time_step));
      OutputMohoParticles((i + 1 * time_step));
      // DumpStatisticStressBox (0.95, 1.0, 0.45, 0.55, i);
      // DumpStatisticStressBox (0.90, 0.95, 0.45, 0.55, i);
      // DumpStatisticStressBox (0.85, 0.90, 0.45, 0.55, i);
      // DumpStatisticStressBox (0.8, 0.85, 0.45, 0.55, i);
      // DumpStatisticStressBox (0.75, 0.8, 0.45, 0.55, i);
      // DumpStatisticStressBox (0.7, 0.75, 0.45, 0.55, i);
      // DumpStatisticStressBox (0.65, 0.7, 0.45, 0.55, i);
      // DumpStatisticStressBox (0.6, 0.65, 0.45, 0.55, i);
      // DumpStatisticStressBox (0.55, 0.60, 0.45, 0.55, i);
      // DumpStatisticStressBox (0.5, 0.55, 0.45, 0.55, i);

      break;

    case 30: // Rift and dykes

      if (i == 0) {
        ChangeRelaxThreshold(0.05);
        ApplyGravity(0);
        Relaxation();
        //	UpdateElle();
        // for (j=0;j<10;j++)
        // ViscousRelax(1,300 * 365*24*60*60);
        // ChangeRelaxThreshold(5);
      }

      Fluid_Parameters(0, 1, 0.000001);

      // time_step = 100; //years , have to use 2000 times steps to get 2 mio
      // years if time is 1000

      // DeformLatticeXMid(-0.0000002, 1, 0.2, 0.4, 0.6);
      // DeformLatticeXMid(-0.0000002, 1, 0.4, 0.43, 0.57);
      // DeformLatticeXMid(-0.0000002, 1, 0.6, 0.48, 0.52);

      // DeformLatticeExtend (0.0000001, 1);

      // ViscousRelax(1,time_step*365*24*60*60);

      for (iii = 0; iii < 80; iii++)
        Input_Fluid_Lattice(10000, 50, iii);

      Calculate_Fluid_Pressure(1, 0, 1, 1, 1, 0.0);

      Relaxation();

      Save_Particle_Velocity();

      UpdateElle();

      OutputTopParticles((i + 1 * time_step));
      OutputMohoParticles((i + 1 * time_step));
      DumpStatisticStressBox(0.95, 1.0, 0.45, 0.55, i);
      DumpStatisticStressBox(0.90, 0.95, 0.45, 0.55, i);
      DumpStatisticStressBox(0.85, 0.90, 0.45, 0.55, i);
      DumpStatisticStressBox(0.8, 0.85, 0.45, 0.55, i);
      DumpStatisticStressBox(0.75, 0.8, 0.45, 0.55, i);
      DumpStatisticStressBox(0.7, 0.75, 0.45, 0.55, i);
      DumpStatisticStressBox(0.65, 0.7, 0.45, 0.55, i);
      DumpStatisticStressBox(0.6, 0.65, 0.45, 0.55, i);
      DumpStatisticStressBox(0.55, 0.60, 0.45, 0.55, i);
      DumpStatisticStressBox(0.5, 0.55, 0.45, 0.55, i);

      break;

    case 31: // Advection diffusion and reaction	large scale
      //--------------------------------------------------------
      // A is solid
      // Temperature is fluid pressure
      // B is concentration
      // C is fluid pressure gradient
      // Energy is darcy velocity in y
      // Dislocation density is porosity
      //----------------------------------------------------------

      // for (iii=0;iii<30;iii++)
      // Fluid_Insert_Random_Nodexy(0.0, 0.1, 0.45, 0.55, 40000);

      // Fluid_Insert_Random_Nodexy(0.8, 1.0, 0.4, 0.6, 20000);

      // set fluid parameters
      Fluid_Parameters(0, 1, 0.00005);

      // calculate the fluid pressure
      // first upper boundary, second lower boundary, the itterations,
      // then sides wrapping then use average concentration

      // Calculate_Fluid_Pressure(1, 10100, 1, 0, 0.5, 0.0);
      Calculate_Fluid_Pressure(10, 10000, 1, 0, 0, j * 10000);
      // Calculate advection, then diffusion

      Calculate_Advection(time_step * 200, size, 2, 1);

      Calculate_Diffusion(0.000000001, time_step * 100, size); // was 5000000

      // Calculate_Advection(time_step, size, 2);
      // Calculate_Diffusion(0.0001, time_step, size);

      // Calculate the reaction

      // Make_Dolomite(1000000, i, 0.2, 0.01);

      // Relaxation();

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
