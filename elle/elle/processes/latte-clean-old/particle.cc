/******************************************************
 * Spring Code Latte 2.0
 *
 * Functions for particle class in particle.cc
 *
 * Class lists defined in Lattice header
 *
 * Daniel Koehn and Jochen Arnold
 * Feb. 2002 to Feb. 2003
 *
 * Daniel Dec. 2003
 *
 * We thank Anders Malthe-Srenssen for his enormous help
 * and introduction to these codes
 *
 * new version daniel/Till 2004/5
 * koehn_uni-mainz.de
 *
 * Clean Version 4.0
 * Erlangen 2025
 *
 ******************************************************/
#include <iostream>

#include "particle.h"

using std::cout;
using std::endl;

// CONSTRUCTOR
/****************************************************
 * Does not do much , defines radius
 * sets flags and some variables to 0
 *
 * Daniel spring 2002
 **************************************************/

Particle::Particle()
    :

      //-------------------------------------------
      // variable definitions
      //-------------------------------------------

      radius(0.005) // set radius default only for x = 100 ! now scaled
                    // automatically

{
  angle_of_internal_friction = 3.141526 / 6.0; // is 30 degrees
  delta_surfE = 0.0;
  damage = 0;
  phase = 0;
  healing = 0;
  surfE = 0;
  disconst = 0.0;
  fracture_time = 0;
  grain = 0;
  fluid_pressure_gradient = 0.0;
  average_porosity = 0.0;
  average_permeability = 0.0;
  no_reaction = false;
  crack_seal_counter = 0.0;
  neig_of_two_spinel_grains = false;
  react_time = 0;

  particle_age = 0;

  newly_overgrown = 0;
  aperture = 0;
  overgrown_area = 0;

  age = 2000;

  neig3 = NULL;
  is_moho_particle = false;
  prob = 0.0;

  fix_y = false; // particle can move in y direction
  fix_x = false; // particle can move in x direction

  done = false; // flag for done Relaxation of particle

  mineral = 1;    // default number for mineral
  conc = 44000.0; // default concentration
  fluid_P = 0;    // default fluid pressure

  merk = 0;
  fracture_reaction = 0.0;
  fluid_particles = 0; // for lattice gas

  isUpperBox = true; // for stylolites

  rate_factor = 1.0; // for distribution in dissolution

  isHole = false;                 // flag for dissolution etc.
  isHoleBoundary = false;         // flag for hole boundaries
  isHoleInternalBoundary = false; // particle of fluid along interface

  box_pos = -1; // no boxposition is -1

  draw_break = false; // particle has no broken bond
  neigh = 0;          // no neighbour breaks

  next_inBox = NULL; // pointer  to  next  in  box

  young = 1.0;                   // default youngs modulus
  springconst = sqrt(3.0) / 2.0; // formula for spring constant

  grain = -1; // no grain is -1
  grain_attribute = -1;

  mV = 1; // just in case for Energy

  is_boundary = false;         // default not a grainboundary particle
  is_lattice_boundary = false; // boundary of box for boundary conditions
  is_left_lattice_boundary = false;
  is_right_lattice_boundary = false;
  is_top_lattice_boundary = false;

  movement = 0;
  temperature = 0;

  nofluid = false; // for lattice gas fracture walls

  //---------------------------------------------
  // if 1.0 is 10 GPa then viscosity 1e11 means
  // viscosity of 1e21 Pas
  //---------------------------------------------

  viscosity = 1e11; // set viscosity default
  previous_viscosity = viscosity;

  //---------------------------------------------
  // set all neighbours to zero
  //---------------------------------------------

  for (i = 0; i < 8; i++) {
    neigP[i] = 0;        // pointer  to  neighbour
    break_Str[i] = 0.0;  // tensile normal stress for spring
    break_Sp[i] = false; // does spring break ?
    no_break[i] = false; // springs can break
    novisc = false;
    fluid_particle[i] = 0;     // lattice gas fluid particles
    fluid_particle_new[i] = 0; // lattice gas fluid particles
  }

  // set list for elle nodes to -1 = empty slot

  for (i = 0; i < 32; i++) // maximum of 32 Elle_nodes per particle at moment
  {
    elle_Node[i] = -1;
  }

  // the value of the endpoint of every spring, used for the viscous deformation
  for (i = 0; i < 6; i++) {
    xy[i][0] = cos(double(i) / 3.0 * 3.1415927);
    xy[i][1] = sin(double(i) / 3.0 * 3.1415927);
  }
  // initial x/y-position coordinates in conjugate radii, uses the unit-circle
  // (for angular springs)
  ax = 0;
  ay = 1.0;
  bx = 1.0;
  by = 0;

  right_lattice_wall_pos = 1.0;
  left_lattice_wall_pos = 0.0;

  transform = false;

  visc_flag = false; // for debugging

  wrapping = false;

  cluster_nb = 0;

  use_gravitation = false;

  fg = 0;
}

// RELAX (double,Particle**,int,bool)
/*****************************************************
 * relaxation routine
 * returns true if the particle is relaxed, i.e. not
 * moving anymore.
 * gets the relaxationthreshold from the lattice
 * and the boxlist and x dimension
 * the routine calculates forces from neighbours
 * and repulsive forces from nonconnected neighbours
 * and calculates how much particle will move to find
 * equilibrium. Then the particle is overrelaxed.
 *
 * Daniel spring, summer 2002
 *
 * added interactive repulsion constant, caluclated
 * from the youngs moduli of particles
 *
 * daniel March, 2003
 ****************************************************/

//----------------------------------------------------------------------------
// function Relax in particle class
//
// receives relaxationthreshold, box list from lattice, size of x direction
//  positions of the walls and wall constant, whether or not code is viscous
//
// called form Lattice::FullRelax()
//-----------------------------------------------------------------------------

bool Particle::Relax(double relaxthresh, Particle **list, int size, bool wall,
                     double rightwall, double leftwall, double lowerwall,
                     double upperwall, double wallconstant, bool length,
                     int debug_nb, int visc_rel, bool sheet, bool grav,
                     bool grav_Press, int part_nr, int numb) {
  //------------------------------------------------
  // local variables
  //------------------------------------------------

  int i, j;            // counter
  int x, y, z;         // for list
  double dx, dy;       // x and y differences
  double dd, alen;     // distances (displacement) between particles, where
                       // alen=radius+radius
  double uxi, uyi;     // unitdistance (x and y)
  double fn, fx, fy;   // forces(strain) (normal, x and y)
  double anc;          // move xx and yy
  double fover = 1.1;  // overrelaxation factor
  double rep_constant; // repulsion constant
  int box, pos;        // box  position helper
  int zaehl;
  double fac;

  double nx, ny, unx, uny, f, usx, usy, ux, uy;

  //----------------------------------------------------
  // set some to Zero
  //----------------------------------------------------

  xx = 0;
  yy = 0;
  anc = 0;

  //------------------------------------------------------
  // dont relax particle itself in Rep list
  //------------------------------------------------------

  done = true; // that is myself (everything that is done is not used in repbox

  //-----------------------------------------------------------
  // loop through all connected neighbours of particle
  //-----------------------------------------------------------

  z = 0;

  for (i = 0; i < 8; i++) // try all springs
  {
    if (neigP[i]) // if neighbour is connected
    {

      zaehl++; // count

      //--------------------------------------------------------------
      // get distance to neighbour
      //--------------------------------------------------------------

      dx = neigP[i]->xpos - xpos; // distance in x
      dy = neigP[i]->ypos - ypos; // distance in y

      dd = sqrt((dx * dx) + (dy * dy)); // diagonal without sign

      //--------------------------------------------------------------
      // get equilibrium length
      // neighpos memorizes the distance to neighbour i in x and y
      // this is used as equilibrium length
      // good to calculate also angular forces
      //--------------------------------------------------------------

      ux = dx - neigpos[i][0]; // difference beteen distance now and equilibrium
      uy = dy - neigpos[i][1]; // difference beteen distance now and equilibrium

      if (dd != 0.0) // safety check
      {
        nx = dx / dd; // unit vector x
        ny = dy / dd; // unit vector y
      } else {
        nx = 0;
        ny = 0;
      }

      f = ux * nx + uy * ny; // overall strain difference times direction

      unx = nx * f; // normal strain in x
      uny = ny * f; // normal strain in y

      usx = ux - unx; // shear strain in x
      usy = uy - uny; // shear strain in y

      // *** normal force - strain times spring constant
      fx = unx * springf[i];
      fy = uny * springf[i];

      // add up forces
      xx = xx + fx;
      yy = yy + fy;
      anc = anc + springf[i];

      // *** shear force - strain times shear constant
      fx = usx * springs[i];
      fy = usy * springs[i];

      // add up forces
      xx = xx + fx;
      yy = yy + fy;
      anc = anc + springs[i];

      neigP[i]->done = true; // this one was done
    }
  }

  //-----------------------------------------------------------------------
  // now do the Repulsion box, try all next neighbours,
  // if they are already done (i.e. had connection) and if they apply
  // a repulsion on particle
  //-----------------------------------------------------------------------

  if (!list[box_pos]) // means I am not in the list, which should not happen
  {
    cout << " problem replist " << endl;
    ;
  } else // all is fine, go in the list
  {

    //---------------------------------------------------
    // nine positions in 2d box 3x3
    //---------------------------------------------------

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {

        //------------------------------------------
        //  case  upper  or  lower  row (+-size)
        //  we now define the 9 box postions.
        //  we are looking in our own and the
        //  surrounding postions for particles that
        //  are not connected with springs but may
        //  still push into us
        //------------------------------------------

        if (j == 0)
          box = 2 * size * (-1);
        if (j == 1)
          box = 0;
        if (j == 2)
          box = 2 * size;
        //-----------------------------------------
        //  now  check  this  position
        //  all  particles  in  this  position !!
        //-----------------------------------------

        pos = box_pos + (i - 1) + box;

        //-----------------------------------------
        // repbox should always be positive
        // the box has direct pointers to particles
        //-----------------------------------------

        if (list[pos]) // if there is somebody
        {
          neig = list[pos]; // helppointer  for  while  loop

          while (neig) // as long as there is somebody
          {

            if (!neig->done) // if not done already
            {

              if (neig->xpos) // just in case
              {
                dx = neig->xpos - xpos;

                dy = neig->ypos - ypos;

                dd = sqrt((dx * dx) + (dy * dy));

                //------------------------------------------------
                // get equilibrium length
                //------------------------------------------------

                alen = (radius + neig->radius);

                //----------------------------------------
                // determine the unitlengths in x and y
                //----------------------------------------

                if (dd != 0.0) {
                  uxi = dx / dd;
                  uyi = dy / dd;
                } else {
                  cout << " zero divide in Relax " << endl;
                  uxi = 0.0;
                  uyi = 0.0;
                }

                //-----------------------------------------------------------
                // determine the force (strain) on particle in RepBox = Young's
                // modulus (grain) * strain
                //-----------------------------------------------------------

                rep_constant = (springconst + neig->springconst) / 2.0;

                if (springconst == 0.0)
                  rep_constant = 0.0;
                if (neig->springconst == 0.0)
                  rep_constant = 0.0;

                fn = rep_constant * (dd - alen); // normal force

                if (fn <
                    0.0) // only compressive forces, only now action is taken
                {
                  fx = fn * uxi; // force (strain) x
                  fy = fn * uyi; // force (strain) y

                  //------------------------------------------------------------
                  // add forces (strain) in xx and yy and springconstants
                  // repulsion constant 1.0 at moment
                  //------------------------------------------------------------

                  xx = xx + fx;
                  yy = yy + fy;
                  anc = anc + rep_constant;
                }
              }
            }
            if (neig->next_inBox) // if there is anotherone in box
            {
              neig = neig->next_inBox; // shift the help pointer
            } else {
              break; // if not go out to next box
            }
          }
        }
      }
    }
  }

  //--------------------------------------------------
  // Get forces (strain) from the viscous sheet
  //--------------------------------------------------

  if (sheet) {
    dx = sheet_x - xpos;
    dy = sheet_y - ypos;
    dd = sqrt((dx * dx) + (dy * dy));

    if (dd != 0.0) {
      uxi = dx / dd;
      uyi = dy / dd;
    } else {
      uxi = 0.0;
      uyi = 0.0;
    }

    fn = rep_constant * 0.00001 * dd;

    fx = fn * uxi; // force (strain) x
    fy = fn * uyi; // force (strain) y

    //------------------------------------------------------------
    // add forces (strain) in xx and yy and springconstants
    // repulsion constant 1.0 at moment
    //------------------------------------------------------------

    xx = xx + fx;
    yy = yy + fy;

    anc = anc + 1;
  }

  //---------------------------------------------------------
  // set all particle flags back to not done
  //---------------------------------------------------------

  done = false; // myself

  for (i = 0; i < 8; i++) // and all my neighbours
  {
    if (neigP[i]) {
      neigP[i]->done = false;
    }
  }

  //========================================================
  // Gravity
  //
  // if gravity is defined the particles can have an
  // extra force due to gravity. there is two functions
  // for this. Either fg as a gravity force that is applied
  // for each particle or using the Add Gravitation function
  //
  //========================================================

  if (fg != 0.0) {
    yy = yy - fg;
    anc += springconst;
  }

  // second gravity function

  if (grav)
    yy = yy - (Add_Gravitation() / size);

  //===========================================================
  // Here the fluid pressure gradients in x and y can be added
  // these have to be rescaled from Pascal to the
  //===========================================================

  xx += (F_P_x) / (10000000000 * size);
  yy += (F_P_y) / (10000000000 * size);

  //-----------------------------------------------------------------
  // now devide forces (strain) again by springconstants (anc)
  // to get movement of particle
  //-----------------------------------------------------------------

  if (anc != 0.0) {
    anc = 1.0 / anc;
    xx = xx * anc; // here strain is calculated from accumulated force strain =
                   // stress / Young
    yy = yy * anc;
  } else {
    xx = 0.0;
    yy = 0.0;
  }

  //=====================================================================
  // This applies to boundaries where x or y are fixed, so they are zero
  //=====================================================================

  if (fix_x)
    xx = 0;
  if (fix_y)
    yy = 0;

  //=====================================================================
  // Now determine the movement that is used for the relaxation
  //=====================================================================

  movement = sqrt(xx * xx + yy * yy);

  //============================================================
  // if movement is too large slow it down
  // a WORKING slow-down mechanism)
  //============================================================

  if (movement > 1.0 * radius) {
    fac = sqrt(0.1 * radius * radius) / movement;
    xx *= fac;
    yy *= fac;
    movement = sqrt(xx * xx + yy * yy);
  }

  //---------------------------------------------------------
  // now return flag if movement larger than relaxthreshold
  //---------------------------------------------------------

  if (relaxthresh < sqrt((xx * xx) + (yy * yy))) // we are still relaxing
  {
    xpos = xpos + xx * fover; // overrelax
    ypos = ypos + yy * fover; // overrelax

    return false;
  } else // particle is relaxed
  {
    return true;
  }
}

double Particle::Add_Gravitation() {

  double g = 9.81, yy = 0.0, norm_stress, strain;

  norm_stress = real_density * g * 2.0 * real_radius * (1 - average_porosity);
  strain = norm_stress / real_young; // E = stress/strain

  yy = strain;
  yy = strain * 1.5;

  return (yy);
}

double Particle::Press_gravitation(int curr_nr, int par_num) {
  double stress, strain;

  // stress = real_density * 9.81 * real_radius;		// stress =
  // force/area = mg/A = rho*g*(r) * 2, bz here consider weight of unit particle
  // with diameter=2*r

  if (curr_nr > par_num)
    stress = real_density * 9.8 * (real_radius + 1000.0);
  else
    stress = real_density * 9.8 * real_radius;
  // cout << curr_nr << " "<< par_num << " ";

  strain = stress / real_young; // E = stress/strain
  // strain /= 1000.0;
  // strain /= par_num;

  // cout << nb << " " << stress << endl;

  // strain = strain * 1.5;
  return (strain);
}

// RELAX (particle **, int)
/*****************************************************
 * Routine that calculates stresses and which
 * bonds will break
 * gets Repulsion box from lattice plus x size
 *
 *
 * Daniel spring and summer 2002
 *
 * Daniel and Jochen scale bug fixed
 * Feb. 2003
 *
 * added active repulsion constant calculated from
 * the youngs modulus of particles
 *
 * daniel March 2003
 ****************************************************/

//-----------------------------------------------------------
// function Relax (bool) in particle class
//
//
//  gets RepulsionBox and x size
//
//
// called form Lattice::FullRelax()
//------------------------------------------------------------

bool Particle::Relax(Particle **list, int size, int visc_rel, bool wall,
                     double rightwall, double leftwall, double lowerwall,
                     double upperwall, double wallconstant) {
  //------------------------------------------------
  // local variables
  //------------------------------------------------

  int i, j, k;          // counter
  double dx, dy;        // x and y differences
  double dd, alen;      // distances between particles
  double rep_constant;  // repulsion constant
  double fn, fx, fy;    // forces
  double uxi, uyi;      // unitdistance
  double area;          // area of particle
  double ten_break;     // help for breaking strenght
  int box, pos, spring; // help for Repulsion

  double lambda, my, poisson;

  bool bbreak = false;

  double nx, ny, unx, uny, f, usx, usy, ux, uy, ratio;

  //----------------------------------------------------
  // set some to Zero
  //----------------------------------------------------

  xx = 0.0;
  yy = 0.0;
  sxx = 0.0;
  syy = 0.0;
  sxy = 0.0;
  exx = 0.0;
  eyy = 0.0;
  exy = 0.0;

  mbreak = 0.0;

  done = true;

  spring = 0;

  // TODO: add gravitation

  //-----------------------------------------------------------
  // loop through all connected neighbours of particle
  //-----------------------------------------------------------

  for (i = 0; i < 8; i++) // try all
  {
    if (neigP[i]) // if neighbour is connected
    {
      spring++;
      // get distance to neighbour
      //--------------------------------------------------------------

      dx = neigP[i]->xpos - xpos;
      dy = neigP[i]->ypos - ypos;

      dd = sqrt((dx * dx) + (dy * dy));

      //--------------------------------------------------------------
      // get equilibrium length
      //--------------------------------------------------------------

      neig = neigP[i];

      ux = dx - neigpos[i][0];
      uy = dy - neigpos[i][1];

      // dd = sqrt((dx*dx) + (dy*dy));

      alen =
          sqrt(neigpos[i][0] * neigpos[i][0] + neigpos[i][1] * neigpos[i][1]);

      ratio = alen / (1.0 / size);

      // cout << ratio << endl;

      // ux = ux * alen/0.01;
      // uy = uy * alen/0.01;

      double nx, ny;

      if (dd != 0.0) {
        nx = dx / dd;
        ny = dy / dd;
      } else {
        nx = 0;
        ny = 0;
      }

      f = ux * nx + uy * ny;
      unx = nx * f * ratio;
      uny = ny * f * ratio;

      usx = ux - unx;
      usy = uy - uny;

      // *** normal force
      fx = unx * springf[i];
      fy = uny * springf[i];

      xx = xx + fx;
      yy = yy + fy;

      sxx += nx * fx;
      syy += ny * fy;
      sxy += nx * fy;

      // *** shear force
      fx = usx * springs[i];
      fy = usy * springs[i];

      xx = xx + fx;
      yy = yy + fy;

      sxx += nx * fx;
      syy += ny * fy;
      sxy += nx * fy;

      neigP[i]->done = true;
    }
  }

  if (!list[box_pos]) {
    // cout << " problem replist " << endl;
    ;
  } else {
    //------------------------------------------------------
    // local box again 3 time 3
    //------------------------------------------------------

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        //----------------------------------------------
        // define box positions for different cases
        //----------------------------------------------

        if (j == 0)
          box = 2 * size * (-1);
        if (j == 1)
          box = 0;
        if (j == 2)
          box = 2 * size;

        pos = box_pos + (i - 1) + box; // define box position

        if (pos > 0) {
          neig = 0;
          neig = list[pos]; // get particle at that position

          while (neig) {

            if (!neig->done) // if not already done
            {
              if (neig->xpos) // just in case
              {
                //--------------------------------------------------------------
                // get distance to neighbour
                //--------------------------------------------------------------

                dx = neig->xpos - xpos;
                dy = neig->ypos - ypos;

                dd = sqrt((dx * dx) + (dy * dy));

                //--------------------------------------------------------------
                // get equilibrium length
                //--------------------------------------------------------------
                if (visc_rel == 1) {
                  alen = Alen(1, dx, dy) + Alen(0, dx, dy);
                  // alen = radius + neig->radius;
                } else
                  alen = radius + neig->radius;

                //--------------------------------------------------------------
                // determine the unitlengths in x and y
                //--------------------------------------------------------------

                if (dd != 0.0) {
                  uxi = dx / dd;
                  uyi = dy / dd;
                } else {
                  cout << " zero divide in Relax " << endl;
                  uxi = 0.0;
                  uyi = 0.0;
                }
                //-----------------------------------------------------------
                // determine the force on particle = constant times strain
                //-----------------------------------------------------------

                rep_constant = (young + neig->young) / 2.0;

                if (young == 0.0)
                  rep_constant = 0.0;
                if (neig->young == 0.0)
                  rep_constant = 0.0;

                // and rescale to spring constant

                rep_constant = rep_constant * sqrt(3.0) / 2.0;

                fn = rep_constant * (dd - alen); // normal force
                if (fn < 0.0)                    // only compressive forces
                {
                  fx = fn * uxi; // force x
                  fy = fn * uyi; // force y

                  //------------------------------------------------------------
                  // add forces in xx and yy and springconstants
                  // repulsion constant 1.0 at moment
                  //------------------------------------------------------------

                  xx = xx + fx;
                  yy = yy + fy;
                  sxx = sxx + fx * radius * uxi;
                  syy = syy + fy * radius * uyi;
                  sxy = sxy + fx * radius * uyi;
                  // anc = anc + rep_constant;
                }
              }
            }
            if (neig->next_inBox) // if there is another one in box
            {
              neig = neig->next_inBox; // move pointer
            } else {
              break;
            }
          }
        }
      }
    }
  }

  // Gravity
  if (fg != 0.0) {

    // syy = syy - (fg*radius);
    // yy = yy - gravforce;
    // anc += springconst;
  }

  // if (grav)
  // yy = yy - (Add_Gravitation( )/size);

  // xx += (F_P_x)/10000000000;
  // yy += (F_P_y)/10000000000;

  // anc += springconst;

  // done = false;

  // for (i=0;i<8;i++) {
  // if (neigP[i]) {
  // neigP[i]->done = false;
  //}
  //}

  //----------------------------------------------------------------
  // rescale stress by area of particle
  //----------------------------------------------------------------

  // area wird in viscous routine constant gehalten, daher kann sie weiterhin
  // so benutzt werden.

  // area = radius*radius*3.1415926;
  // area = average_radius*average_radius*3.1415926;

  // area = (1/(size*2.0))*(1/(size*2.0))*3.14;

  // area = radius*radius*3.1415926;
  area = 2.0 * radius;

  // area = 0.005 * 0.005 * 3.14;

  // area = area * 6.0/spring;
  // area = 2.0*radius *1.0;

  sxx = sxx / area;
  syy = syy / area;
  sxy = sxy / area;

  //~ cout << sxx << " " << syy << " " << sxy << endl;
  // now calculate strain tensor from stress tensor
  poisson = 1.0 / 3.0;
  lambda = -(poisson * young) / (2.0 * poisson * poisson + poisson - 1.0);
  my = young / (2.0 * poisson + 2.0);

  exx = ((sxx - syy) * lambda + 2.0 * my * sxx) /
        (4.0 * my * lambda + 4.0 * my * my);
  eyy = -((sxx - syy) * lambda - 2.0 * my * syy) /
        (4.0 * my * lambda + 4.0 * my * my);
  exy = sxy / (2.0 * my);

  // and the principal strains
  //  e1 = 0.5 * (exx + eyy) - ( exy*exy + 0.25 * (exx - eyy) * (exx - eyy) );
  //  e2 = 0.5 * (exx + eyy) + ( exy*exy + 0.25 * (exx - eyy) * (exx - eyy) );

  // and the volumetric strain
  ev = exx + eyy;

  return bbreak;
}

// CHANGEBOX (particle **, int)
/**********************************************************************
 * routine that changes the box for a particle if it moved out of
 * its box during the relaxation routine or a deformation
 * First have to remove particle from box depending on where it is
 * in the box list and reconnect the box list of remaining particles
 * and then we have to add the particle to the end of the list
 * in the new box.
 *
 * Daniel spring 2002
 **********************************************************************/

//--------------------------------------------------------------------------
// function ChangeBox in particle class
//
// receives box list from Lattice and the x size.
//
// called from relaxation routines in Lattice and from deformation routines
//--------------------------------------------------------------------------

void Particle::ChangeBox(Particle **list, int size) {
  int x, y; // some variables

  //------------------------------------------------
  // again find current position of particle in box
  // is a double check
  //------------------------------------------------

  x = int(xpos * size);
  y = int(ypos * size);
  x = (y * size * 2) + x;

  // if (x<0) x = 0;

  if (x != box_pos) // if current position not the same as saved position
  {
    // cout << "not boxpos"<< endl;
    //--------------------------------------------------------------
    //   first remove particle from old position
    //--------------------------------------------------------------

    if (list[box_pos]) {
      neig = neig2 = list[box_pos];

      if (neig->nb == nb) // if we are the first
      {
        list[box_pos] = neig->next_inBox; // the next is the first
      } else                              // if not loop
      {
        // while(neig!=this)  // loop until particle is found
        while (neig->nb != nb) {
          if (neig->next_inBox) // if not the particle is not there in box
          {
            neig2 = neig;            // have a pointer to the last particle
            neig = neig->next_inBox; // set pointer to the next
          } else {
            cout << "Problem in RepBox, changeBox  " << nb << endl;
            break;
          }
        }
        //-------------------------------------------------------
        // in here means we found the particle. Now the previous
        // particle is made to point to the next one that is
        // we just take our particle out of the list
        //-------------------------------------------------------

        neig2->next_inBox = neig->next_inBox;
      }
      //--------------------------------------------------------------------
      // next inBox is empty because we put particle to the end of the list
      //--------------------------------------------------------------------

      next_inBox = NULL;
    } else {
      cout << "problem in repBox, ChangeBox  " << nb << endl;
    }

    //----------------------------------------------------------------
    //  put  in  new  position  at  end  of  list
    //  neig  now  points  at  me =  at this particle
    //----------------------------------------------------------------

    // *** 8.1.04: Line 843: Verbesserung eines bugs aus alter mike-version: ***

    box_pos = x;

    if (!list[box_pos]) {
      list[box_pos] = neig; // in case the box is empty
    } else {
      neig2 = list[box_pos]; // for the last particle

      while (neig2->next_inBox) // go to end of list
      {
        neig2 = neig2->next_inBox;
      }

      neig2->next_inBox = neig; // and add particle
    }
  }
}

// the distance in X between 2 Particles
// takes wrapping into account
float Particle::GetXDist(Particle *center, Particle *peripher) {
  if (!((center->is_right_lattice_boundary &&
         peripher->is_left_lattice_boundary) ||
        (center->is_left_lattice_boundary &&
         peripher->is_right_lattice_boundary))) {
    return (peripher->xpos - center->xpos);
  } else if (center->is_right_lattice_boundary &&
             peripher->is_left_lattice_boundary) {
    return ((peripher->xpos + right_lattice_wall_pos) - (center->xpos));
  } else if (peripher->is_right_lattice_boundary &&
             center->is_left_lattice_boundary) {
    return (peripher->xpos - (right_lattice_wall_pos + center->xpos));
  } else
    cout << "f..."; // meraculous work ;)
}

// SET_SPRINGS
/***********************************************************
 * function that gives the springs a constant
 * and a breaking strength
 * checks which neighbours are there.
 * gives just one constant at the moment,
 * should be overwritten later to make grain boundaries
 *
 * Daniel spring 2002
 * fixed breaking strength bug
 * Daniel and Jochen Feb. 2003
 *********************************************************/

void Particle::SetSprings(double E, double poisson) {
  int i;

  Particle::poisson_ratio = poisson;
  Particle::E = E;

  average_radius = radius;

  for (i = 0; i < 8; i++) {

    if (neigP[i]) // if the neighbour is there
    {
      springf[i] =
          sqrt(3.0) * E / (3.0 * (1.0 - poisson)); // set spring constant
      springs[i] = (1.0 - 3.0 * poisson) / (1.0 + poisson) *
                   springf[i]; // set spring to 1.0

      springv[i] = 1.0; // viscosity

      break_Str[i] = 0.0017;

      rad[i] = radius; // equilibrium length

      spring_boundary[i] = false;
    }
  }

  for (i = 0; i < 9; i++) {
    rep_rad[i] = radius;
  }
}

// FIND_NEIGHBOUR_F
/*****************************************************
 * Simple Routine to catch a pointer to another
 * particle object in the list of particles
 * that is neigDist from the current particle
 * in the list, returns the pointer
 *
 * Daniel spring 2002
 *****************************************************/

Particle *Particle::FindNeighbourF(int neigDist) {
  neig = nextP;
  for (i = 0; i < (neigDist - 1); i++) {
    neig = neig->nextP;
  }
  return neig;
}

// FIND_NEIGHBOUR_B
/******************************************************
 * same as above but runs backwards in list
 * returns pointer
 *
 * Daniel spring 2002
 *****************************************************/

Particle *Particle::FindNeighbourB(int neigDist) {
  neig = prevP;
  for (i = 0; i < (neigDist - 1); i++) {
    neig = neig->prevP;
  }
  return neig;
}

// CONNECT
/******************************************************
 * This Routine connects all the particles with their
 * Neighbours. Each particle receives pointers to the
 * neighbours.
 * function receives type of lattice and number of
 * particles in x and y direction from lattice
 * function is called from lattice class
 * from the MakeLattice functiion of the lattice object
 * which builds the complete lattice
 *
 * Daniel spring 2002
 *****************************************************/

void Particle::Connect(int lType, int lParticlex, int lParticley,
                       bool wrapping) {

  Particle::wrapping = wrapping;

  if (!wrapping) {

    if (lType == 1) // if triagonal

    /***************************************************
     * now we connect all particles to each other
     * even and uneven rows have different neighbours
     * and all the boundary particles also !!
     ************************************************/

    {
      // left lower corner
      if (nb == 0) {
        neigP[0] = nextP;
        neigP[1] = FindNeighbourF(lParticlex);
      }
      // right lower corner
      else if ((nb + 1) == lParticlex) {
        neigP[1] = FindNeighbourF(lParticlex);
        neigP[2] = FindNeighbourF(lParticlex - 1);
        neigP[3] = prevP;
      }
      // left upper corner
      else if (nb == ((lParticlex * lParticley) - lParticlex)) {
        // check if last row is even or not
        if (lParticley == (2 * (lParticley / 2))) {
          neigP[4] = FindNeighbourB(lParticlex);
          neigP[5] = FindNeighbourB(lParticlex - 1);
          neigP[0] = nextP;
          cout << "even" << endl;
        } else // uneven
        {
          neigP[5] = FindNeighbourB(lParticlex);
          neigP[0] = nextP;
          cout << "uneven" << endl;
        }
      }
      // right upper corner
      else if ((nb + 1) == (lParticlex * lParticley)) {
        // check if last row is even or not
        if (lParticley == (2 * (lParticley / 2))) {
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex);
        } else // uneven
        {
          neigP[3] = prevP;
          neigP[5] = FindNeighbourB(lParticlex);
          neigP[4] = FindNeighbourB(lParticlex + 1);
        }
      }
      // rest of lower row
      else if ((nb + 1) < lParticlex) {
        neigP[0] = nextP;
        neigP[1] = FindNeighbourF(lParticlex);
        neigP[2] = FindNeighbourF(lParticlex - 1);
        neigP[3] = prevP;
      }
      // rest of upper row
      else if ((nb + 1) > ((lParticlex * lParticley) - lParticlex)) {
        if (lParticley == (2 * (lParticley / 2))) {
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex);
          neigP[5] = FindNeighbourB(lParticlex - 1);
          neigP[0] = nextP;
        } else {
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex + 1);
          neigP[5] = FindNeighbourB(lParticlex);
          neigP[0] = nextP;
        }
      }
      // left boundary row
      else if ((nb + 1) == (((nb / lParticlex) * lParticlex) + 1)) {
        if (xpos == 0.0) // rows 0,2,etc.
        {
          neigP[5] = FindNeighbourB(lParticlex);
          neigP[0] = nextP;
          neigP[1] = FindNeighbourF(lParticlex);
        } else // rows 3,5,etc.
        {
          neigP[4] = FindNeighbourB(lParticlex);
          neigP[5] = FindNeighbourB(lParticlex - 1);
          neigP[0] = nextP;
          neigP[1] = FindNeighbourF(lParticlex + 1);
          neigP[2] = FindNeighbourF(lParticlex);
        }
      }
      // right boundary row
      else if ((nb + 1) == (((nb + 1) / lParticlex) * lParticlex)) {
        // rows 0,2 etc.
        if (((nb + 1) / lParticlex) == ((((nb + 1) / lParticlex) / 2) * 2)) {
          neigP[2] = FindNeighbourF(lParticlex);
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex);
        } else // rows 1,3 etc.
        {
          neigP[1] = FindNeighbourF(lParticlex);
          neigP[2] = FindNeighbourF(lParticlex - 1);
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex + 1);
          neigP[5] = FindNeighbourB(lParticlex);
        }
      }
      // the rest
      else {
        // first even numbers (row 1,3,5 etc)
        // since now even numbers are uneven rows !
        if (((nb + 1) / lParticlex) == (2 * (((nb + 1) / lParticlex) / 2))) {
          neigP[0] = nextP;
          neigP[1] = FindNeighbourF(lParticlex);
          neigP[2] = FindNeighbourF(lParticlex - 1);
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex + 1);
          neigP[5] = FindNeighbourB(lParticlex);
        } else // rows 2,4,etc.
        {
          neigP[0] = nextP;
          neigP[1] = FindNeighbourF(lParticlex + 1);
          neigP[2] = FindNeighbourF(lParticlex);
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex);
          neigP[5] = FindNeighbourB(lParticlex - 1);
        }
      }
    }
  } else {
    if (lType == 1) // if triagonal

    /***************************************************
     * now we connect all particles to each other
     * even and uneven rows have different neighbours
     * and all the boundary particles also !!
     ************************************************/

    {
      // left lower corner
      if (nb == 0) {
        neigP[0] = nextP;
        neigP[1] = FindNeighbourF(lParticlex);
        neigP[2] = FindNeighbourF(2 * lParticlex - 1);
        neigP[3] = FindNeighbourF(lParticlex - 1);
      }
      // right lower corner
      else if ((nb + 1) == lParticlex) {
        neigP[1] = FindNeighbourF(lParticlex);
        neigP[2] = FindNeighbourF(lParticlex - 1);
        neigP[3] = prevP;
        neigP[0] = FindNeighbourB(lParticlex - 1);
      }
      // left upper corner
      else if (nb == ((lParticlex * lParticley) - lParticlex)) {
        // check if last row is even or not
        if (lParticley == (2 * (lParticley / 2))) {
          neigP[4] = FindNeighbourB(lParticlex);
          neigP[5] = FindNeighbourB(lParticlex - 1);
          neigP[0] = nextP;
          neigP[3] = FindNeighbourF(lParticlex - 1);
          cout << "even" << endl;
        } else // uneven
        {
          neigP[5] = FindNeighbourB(lParticlex);
          neigP[0] = nextP;
          neigP[3] = FindNeighbourF(lParticlex - 1);
          neigP[4] = prevP;
          cout << "uneven" << endl;
        }
      }
      // right upper corner
      else if ((nb + 1) == (lParticlex * lParticley)) {
        // check if last row is even or not
        if (lParticley == (2 * (lParticley / 2))) {
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex);
          neigP[5] = FindNeighbourB(2 * lParticlex - 1);
          neigP[0] = FindNeighbourB(lParticlex - 1);
        } else // uneven
        {
          neigP[3] = prevP;
          neigP[5] = FindNeighbourB(lParticlex);
          neigP[4] = FindNeighbourB(lParticlex + 1);
          neigP[0] = FindNeighbourB(lParticlex - 1);
        }
      }
      // rest of lower row
      else if ((nb + 1) < lParticlex) {
        neigP[0] = nextP;
        neigP[1] = FindNeighbourF(lParticlex);
        neigP[2] = FindNeighbourF(lParticlex - 1);
        neigP[3] = prevP;
      }
      // rest of upper row
      else if ((nb + 1) > ((lParticlex * lParticley) - lParticlex)) {
        if (lParticley == (2 * (lParticley / 2))) {
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex);
          neigP[5] = FindNeighbourB(lParticlex - 1);
          neigP[0] = nextP;
        } else {
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex + 1);
          neigP[5] = FindNeighbourB(lParticlex);
          neigP[0] = nextP;
        }
      }
      // left boundary row
      else if ((nb + 1) == (((nb / lParticlex) * lParticlex) + 1)) {
        if (xpos == 0.0) // rows 0,2,etc.
        {
          neigP[5] = FindNeighbourB(lParticlex);
          neigP[0] = nextP;
          neigP[1] = FindNeighbourF(lParticlex);
          neigP[3] = FindNeighbourF(lParticlex - 1);
          neigP[2] = FindNeighbourF(2 * lParticlex - 1);
          neigP[4] = prevP;
        } else // rows 3,5,etc.
        {
          neigP[4] = FindNeighbourB(lParticlex);
          neigP[5] = FindNeighbourB(lParticlex - 1);
          neigP[0] = nextP;
          neigP[1] = FindNeighbourF(lParticlex + 1);
          neigP[2] = FindNeighbourF(lParticlex);
          neigP[3] = FindNeighbourF(lParticlex - 1);
        }
      }
      // right boundary row
      else if ((nb + 1) == (((nb + 1) / lParticlex) * lParticlex)) {
        // rows 1,3 etc.
        if (((nb + 1) / lParticlex) == ((((nb + 1) / lParticlex) / 2) * 2)) {
          neigP[2] = FindNeighbourF(lParticlex);
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex);
          neigP[1] = nextP;
          neigP[5] = FindNeighbourB(2 * lParticlex - 1);
          neigP[0] = FindNeighbourB(lParticlex - 1);
        } else // rows 0,2 etc.
        {
          neigP[1] = FindNeighbourF(lParticlex);
          neigP[2] = FindNeighbourF(lParticlex - 1);
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex + 1);
          neigP[5] = FindNeighbourB(lParticlex);
          neigP[0] = FindNeighbourB(lParticlex - 1);
        }
      }
      // the rest
      else {
        // first even numbers (row 1,3,5 etc)
        // since now even numbers are uneven rows !
        if (((nb + 1) / lParticlex) == (2 * (((nb + 1) / lParticlex) / 2))) {
          neigP[0] = nextP;
          neigP[1] = FindNeighbourF(lParticlex);
          neigP[2] = FindNeighbourF(lParticlex - 1);
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex + 1);
          neigP[5] = FindNeighbourB(lParticlex);
        } else // rows 2,4,etc.
        {
          neigP[0] = nextP;
          neigP[1] = FindNeighbourF(lParticlex + 1);
          neigP[2] = FindNeighbourF(lParticlex);
          neigP[3] = prevP;
          neigP[4] = FindNeighbourB(lParticlex);
          neigP[5] = FindNeighbourB(lParticlex - 1);
        }
      }
    }
  }
}

// BREAKBOND()
/****************************************************************
 * routine that breaks bonds between particle
 * the routine breaks the bond of particle i and also the bond
 * of its connected neighbour (same bond just specified in both
 * particles)
 *
 *
 * Daniel spring 2002
 *
 * Particle with broken bonds is now part of grain boundary,
 * is_boundary flag is set true for both particles
 *
 * Daniel March 2003
 ***************************************************************/

//-----------------------------------------------------------------
// function BreakBond in particle class
//
//-----------------------------------------------------------------

double Particle::BreakBond() {
  int i; // counter
  int j; // Reference for neighbour

  //-----------------------------------------------------------------
  // first break the neighbour
  // find neighbour and break
  //-----------------------------------------------------------------

  // cout << endl;
  // cout << "neigh : " << neigh << endl;
  // cout << "neigP[neigh]->nb : " << neigP[neigh]->nb << endl;

  for (i = 0; i < 8; i++) {
    if (neigP[neigh]->neigP[i]) {
      if (neigP[neigh]->neigP[i]->nb == nb) {
        // cout << "neigP[neigh]->neigP[i]->nb : " << neigP[neigh]->neigP[i]->nb
        // << endl; cout << "nb : " << nb << endl; cout << "i : " << i << endl;
        j = i; // thats the neighbour

        cout << " with neighbour : " << i << endl;
      }
    }
  }

  neigP[neigh]->neigP[j] = 0;      // neighbour gone
  neigP[neigh]->draw_break = true; // particle has a broken bond
  neigP[neigh]->springf[j] = 0.0;  // springconstant is zero
  neigP[neigh]->springs[j] = 0.0;
  neigP[neigh]->break_Str[j] = 0.0;  // breakingstrength is zero
  neigP[neigh]->break_Sp[j] = false; // spring is broken
  neigP[neigh]->is_boundary = true;  // particle is now part of grain boundary
  neigP[neigh]->crack_seal_counter++;
  //-----------------------------------------------------
  // now the particle itself
  //-----------------------------------------------------

  neigP[neigh] = 0;
  draw_break = true;
  springf[neigh] = 0.0;
  springs[neigh] = 0.0;
  break_Str[neigh] = 0.0;
  break_Sp[neigh] = false;
  is_boundary = true;
  mbreak = 0.0;
  crack_seal_counter++;

  //------------------------------------------------------
  // get Moment tensor for fracture
  //------------------------------------------------------

  //-------------------------------------------------------
  // for the particle movement you have to take
  // an initial position that takes into account any
  // deformation related homogeneous movements of particles
  //-------------------------------------------------------

  return 0.0;
}

/***************************************************
 *
 * 	Visco elastic routine changes lengths of springs
 *
 ****************************************************/

double Particle::Alen(int par, double dx, double dy) {

  double l, m, ag, x, y;

  // 	if ((par == 0 && is_lattice_boundary) || (par == 1 &&
  // neig->is_lattice_boundary)) 		return ( radius );

  switch (par) {

  case 0:

    dx *= -1.0;
    dy *= -1.0;

    ag = atan(dy / dx);

    if (dx > 0) {
      ag = ag + 3.14159265;
    } else if (dx < 0 && dy > 0) {
      ag = ag + 2 * 3.14159265;
    }

    if (dx == 0.0) {
      if (dy < 0.0) {
        ag = 3.14159265 / 2.0;
      } else {
        ag = 3.14159265 * 1.5;
      }
    }

    x = bx * cos(ag) - ax * sin(ag);
    y = by * cos(ag) + ay * sin(ag);

    break;

  case 1:

    ag = atan(dy / dx);

    if (dx > 0) {
      ag = ag + 3.14159265;
    } else if (dx < 0 && dy > 0) {
      ag = ag + 2 * 3.14159265;
    }

    if (dx == 0.0) {
      if (dy < 0.0) {
        ag = 3.14159265 / 2.0;
      } else {
        ag = 3.14159265 * 1.5;
      }
    }

    x = neig->bx * cos(ag) - neig->ax * sin(ag);
    y = neig->by * cos(ag) + neig->ay * sin(ag);

    break;
  }

  l = sqrt(x * x + y * y) * radius;
  // 	if (l<(radius/2.0))
  // 		cout << "heil eris" << endl;

  return (l);
}
