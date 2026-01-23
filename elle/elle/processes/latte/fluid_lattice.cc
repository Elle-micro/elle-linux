#include "fluid_lattice.h"
#include <iostream>
#include <stdio.h>

using namespace std;

//-------------------------------------------------------------------------------
// fluid lattice February 2016, Daniel, added Advection/Diffusion 2018 Daniel
//
// based on pressure lattice by Irfan Ghani
// ADI from Till Sachau's temperature diffusion code
//
// Derivation of Pressure diffusion equation based on Renauld Toussaint's work
//
// This code deals with a compressible fluid, has a square lattice that is
// linked to the particle lattice of latte. Fluid lattice cell width is twice
// that of the repulsion box of the particle lattice. Repulsion box is used to
// link both
//
// Fluid lattice deals with variations of fluid pressure per node assuming Darcy
// flow of fluids in the background. Reads in porosity from the solid code,
// converts that into a permeability. Gives back fluid pressure gradients as
// Fluid forces in X and Y that act on the solid.
//
//-------------------------------------------------------------------------------

// ---------------------------------------------------------------
// Constructor of Fluid_Lattice class,
// is called by initialize_fluid_lattice function in lattice.cc
//
// at the moment the constructor does not do much, just defines some
// basic parameters
//
// reads in the lattice size that is passed on by the solid lattice
// = amount of particles along X, this is also the width of the
// repulsion box within the initiial configuration.
// Careful that repulsion box width in the solid is twice X in order
// to leave space on the right hand side for extension.
// ---------------------------------------------------------------

Fluid_Lattice::Fluid_Lattice(int lattice_size)
    :

      // -----------------------------------------------------------------------------
      // default variable values. Variables are declared in Pressure_lattice
      // header
      // -----------------------------------------------------------------------------

      dim(100), // the default grid width, but this is read in from the lattice
                // directly

      box(1) // ?
{
  pl_size = lattice_size; // read in the actual lattice size from the solid
                          // (amount of particles in x)

  dim = lattice_size / 2; // fluid lattice size (dim) is half of lattice size

  width = 1.0 / dim; // width of fluid lattic box, non-dimensional

  time_a = 1; // in sec, one day is 86400, default time 1 second

  mu = 1.0e-03; // viscosity default

  compressibility =
      4.5e-10; // water compressibility = 4.5e10 m.m/N at 25 C°, default

  compressibility_c = -1;
  compressibility_T = 210e-6;

  rho_fluid = 1.0e+03; // water mass density (1000 kg/m³), default

  scale_constt = 100; // box size in meters, default, is changed later on
}

//--------------------------------------------------------------------------------------------
// Boundary and initialization conditions, functions that build up the initial
// config.
//--------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------
// this function just sets the background pressure for the fluid lattice in
// cases where there is no initial gradient.
//
// two fluid pressure matrices to calculate delta P in times
//
// called from lattice (void Lattice::Initialize_Fluid_Lattice(double pressure,
// double scale))
//--------------------------------------------------------------------------------------------

void Fluid_Lattice::Background_Pressure(double pressure, double concentration) {
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      oldPf[i][j] = pressure;
      Pf[i][j] = pressure;
      Con[i][j] = concentration;
      Metal_A[i][j] = concentration;
      Temperature[i][j] = concentration;
      Si[dimX][dimY] = 0;
    }
  }
}

// this was just for debugging, maybe remove

void Fluid_Lattice::Background_Concentration(double concentration) {
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      Con[i][j] = concentration;
    }
  }
}

void Fluid_Lattice::Set_Concentration(int x, int y, double concentration) {
  Con[x][y] = concentration;
}

//---------------------------------------------------------------------------------------------
// Setting the scale (and non-dimensional area of node)
//
// now also includes a change in time so that this has not to be changed in
// fluid lattice
//
// called from lattice (void Lattice::Initialize_Fluid_Lattice(double pressure,
// double scale))
//---------------------------------------------------------------------------------------------

void Fluid_Lattice::Set_Scale(double scale, double time_factor) {
  scale_constt = scale;
  area_node = width * width;
  time_a = time_a * time_factor;
}

//--------------------------------------------------------------------------------------------------
// setting a boundary condition at the start with hydrostatic conditions within
// the box
//--------------------------------------------------------------------------------------------------

void Fluid_Lattice::Set_hydrostatic_gradient(int depth) {
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {

      Pf[i][j] = (depth + (scale_constt * (float((dim - 1) - j) / (dim - 1)))) *
                 1000 * 9.81;
    }
  }
}

//------------------------------------------------------------------------------------------------
// Sending information to the fluid lattice, need porosity or solid fraction per
// node and the solid velocity.
//
// function receives the particle list, which is a pointer to the repulsion box,
// so that we can find the particles easily. Each of the boxes contains pointers
// to the actual particle/particles in the box, so that parameters in the
// particle can be easily accessed.
//
// fracture effect not yet used
//
// The solid fraction is determined by summing up the radii of particles, but
// the fluid radii. Particles have a radius for the elastic part and a DIFFERENT
// radius for the fluid representing the "area" of the particle
//
// called from lattice: void Lattice::Fluid_Parameters()
//-------------------------------------------------------------------------------------------------

void Fluid_Lattice::Read_Porosity_Movement_Matrix(double fracture_effect,
                                                  Particle **list,
                                                  int boundary) {
  int i, j, k, particle_in_box, count;
  double solid_average, vel_x_av, vel_y_av;

  // loop through the matrix, i is x, and j is y (note that Irfan did the
  // opposite in pressure lattice) Each fluid box contains 4 repulsion boxes.

  // calculate an average for the boundary conditions. These boundary conditions
  // turn out to be a major pain there are three different ones at the moment,
  // the average is number 3

  solid_average = 0.0;
  vel_x_av = 0.0;
  vel_y_av = 0.0;
  count = 0;

  for (i = 1; i < dim - 1; i++) // loop in x = i and y = j except for boundaries
  {
    for (j = 1; j < dim - 1; j++) {

      particle_in_box = 0; // count sum of particles in the box (FLUID box)
      rho = 0;             // solid fraction sum up

      vel_x[i][j] = 0; // velocity x matrix
      vel_y[i][j] = 0; // velocity y matrix

      solid_fraction[i][j] = 0.9; // solid fraction matrix

      //------------------------------------------------------------------------------------------------
      // In order to get the right solid fraction and to avoid grid effects
      // converting a triangular solid lattice into the square fluid lattice we
      // use a smoothing function and consider solid in the four repulsion boxes
      // of the fluid node plus the surrounding 12 repulsion boxes. The solid is
      // then weight using the smoothing function depending on how far the
      // particle is from the centre of the pressure node.
      //
      // this is also applied to the velocities.
      //-------------------------------------------------------------------------------------------------

      //-------------------------------------------------------------------------------------------------
      // Because of the above the outer rows of the matrix have to be treated
      // separately, otherwise we look into boxes that dont exist.
      //-------------------------------------------------------------------------------------------------

      // if((i > 0 && i < dim-1) && (j > 0 && j < dim-1))
      {
        //--------------------------------------------------------------------------------------------
        // Now define counters for all 16 repulsion boxes that we have to look
        // into to define the solid fraction for the pressure node as a function
        // of i and j. start in the lower left hand corner with k = 0, then go
        // to the right. So 0,1,2,3 are below the actual pressure node, 5,6,9,10
        // are in the actual pressure node.
        //
        // the repulsion box is a one-dimensional list starting in the lower
        // left hand corner of the lattice with 0 and running towards the right
        // with TWICE the lattice width (so for a 100 lattice up to 199). The
        // box in the second horizontal row above 0 is then 200 and so on.
        //
        // at the moment this routine is not using the boundaries because the
        // boxes there are either not completely full or are empty. There is
        // already the try to double up full boxes for the boundaries, but that
        // does not work yet, needs further debugging
        //--------------------------------------------------------------------------------------------

        for (k = 0; k < 16; k++) {
          if (k == 0) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 1 +
                  ((i - 1) * 2); // pl_size = resolution of particle lattice
            if (i == 0)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 3 + ((i - 1) * 2); // take 2
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 1 +
                    ((i - 1) * 2); // take 8
            if (j == 0 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10

          } else if (k == 1) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 + ((i - 1) * 2);
            if (i == 0)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 4 + ((i - 1) * 2); // take 3
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9
            if (j == 0 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10

          } else if (k == 2) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 3 + ((i - 1) * 2);
            if (i == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 1 + ((i - 1) * 2); // take 0
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10
            if (j == 0 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9

          } else if (k == 3) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 + ((i - 1) * 2);
            if (i == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 2 + ((i - 1) * 2); // take 1
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 4 +
                    ((i - 1) * 2); // take 11
            if (j == 0 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9

          }

          else if (k == 4) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 1 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 3 +
                    ((i - 1) * 2); // take 14
            if (j == 0 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10

          } else if (k == 5) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 4 +
                    ((i - 1) * 2); // take 7
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 2 +
                    ((i - 1) * 2); // take 13
            if (j == 0 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10
          } else if (k == 6) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 1 +
                    ((i - 1) * 2); // take 4
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 3 +
                    ((i - 1) * 2); // take 14
            if (j == 0 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9

          } else if (k == 7) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 4 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 4 +
                    ((i - 1) * 2); // take 15
            if (j == 0 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9

          } else if (k == 8) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 1 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10
            if (j == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 1 + ((i - 1) * 2); // take 0
            if (j == dim - 1 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6

          } else if (k == 9) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 4 +
                    ((i - 1) * 2); // take 11
            if (j == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 2 + ((i - 1) * 2); // take 1
            if (j == dim - 1 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6

          } else if (k == 10) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 1 +
                    ((i - 1) * 2); // take 8
            if (j == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 3 + ((i - 1) * 2); // take 2
            if (j == dim - 1 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5

          } else if (k == 11) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 4 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9
            if (j == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 4 + ((i - 1) * 2); // take 3
            if (j == dim - 1 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5

          } else if (k == 12) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 1 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 3 +
                    ((i - 1) * 2); // take 14
            if (j == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 1 +
                    ((i - 1) * 2); // take 4
            if (j == dim - 1 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6

          } else if (k == 13) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 2 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 4 +
                    ((i - 1) * 2); // take 15
            if (j == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5
            if (j == dim - 1 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6

          } else if (k == 14) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 3 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 1 +
                    ((i - 1) * 2); // take 12
            if (j == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6
            if (j == dim - 1 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5

          } else if (k == 15) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 4 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 2 +
                    ((i - 1) * 2); // take 13
            if (j == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 4 +
                    ((i - 1) * 2); // take 7
            if (j == dim - 1 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5
          }

          //-----------------------------------------------------------
          // If the box is filled deal with the particle
          //-----------------------------------------------------------

          if (list[box]) {
            help = list[box]; // help now points at particle

            //------------------------------------------------------------------------------
            // first the smoothening function, get distance from centre node to
            // particle in x and y
            //------------------------------------------------------------------------------

            x_diff = abs(help->xpos - ((double(i) * width) + (width / 2.0)));
            y_diff = abs(help->ypos - ((double(j) * width) + (width / 2.0)));

            // adjust if too small

            if (y_diff < (help->radius / 2.0))
              y_diff = help->radius / 2.0;

            if (x_diff < (help->radius / 2.0))
              x_diff = help->radius / 2.0;

            // should not be larger than the width

            if (x_diff > width)
              x_diff = width;
            if (y_diff > width)
              y_diff = width;

            // final smooth function

            smooth_func = (1.0 - (x_diff / width)) * (1.0 - (y_diff / width));

            // add up solid fraction

            rho = rho + smooth_func * help->area_par_fluid;

            // add up velocities

            vel_x[i][j] =
                vel_x[i][j] +
                (help->velx * smooth_func); // velocity (flux) per node
            vel_y[i][j] = vel_y[i][j] + (help->vely * smooth_func);

            if (fracture_effect > 0.0) {
              if (help->is_boundary) {
                fracture[i][j] = fracture_effect;
              }
            }

            // count particles

            particle_in_box++;

            //-----------------------------------------------------------------------
            // check if there is an additional particle in the repulsion box
            // the additional particle would be attached to the help particle
            // with the next_inBox pointer.
            //-----------------------------------------------------------------------

            while (help->next_inBox) {
              // and do the same stuff as above.

              x_diff = abs(help->next_inBox->xpos -
                           ((double(i) * width) + (width / 2.0)));
              y_diff = abs(help->next_inBox->ypos -
                           ((double(j) * width) + (width / 2.0)));

              if (y_diff < (help->radius / 2.0))
                y_diff = help->radius / 2.0;
              if (x_diff < (help->radius / 2.0))
                x_diff = help->radius / 2.0;

              if (x_diff > width)
                x_diff = width;
              if (y_diff > width)
                y_diff = width;

              smooth_func = (1.0 - (x_diff / width)) * (1.0 - (y_diff / width));

              rho = rho + smooth_func * help->next_inBox->area_par_fluid;

              vel_x[i][j] =
                  vel_x[i][j] + (help->next_inBox->velx * smooth_func);
              vel_y[i][j] =
                  vel_y[i][j] + (help->next_inBox->vely * smooth_func);

              if (fracture_effect > 0.0) {
                if (help->is_boundary) {
                  fracture[i][j] = fracture_effect;
                }
              }

              particle_in_box++;

              help = help->next_inBox;
            }
          }
        }

        // calculate solid fraction

        solid_fraction[i][j] = rho / area_node;

        solid_average =
            solid_average + solid_fraction[i][j]; // for boundary condition
        count++;

        if (solid_fraction[i][j] >= 1.0) // cannot be more than 1
          solid_fraction[i][j] = 0.999;

        if (solid_fraction[i][j] <= 0.01) // not too small, might check
          solid_fraction[i][j] = 0.01;

        vel_x[i][j] = vel_x[i][j] / rho; // average velocity(flux) per node
        vel_y[i][j] = vel_y[i][j] / rho;

        vel_x_av += vel_x[i][j]; // for boundary condition
        vel_y_av += vel_y[i][j];
      }
    }
  }

  if (boundary == 1) {
    // boundary conditions, apply after matrix filled,
    // condition 1: simply copy the next row from the matrix
    // to the boundary row on sides as well as at bottom and
    // top. Works with random background noise. Becomes
    // problematic when Elle grains are used (boundaries then
    // seem to slow things down because they are mainly oriented 90 degrees
    // to the boundary.
    // condition 2: copy neighbour row to boundary row for upper and lower
    // boundary and wrap rows for left and right boundary. Again can be
    // problematic when grain boundaries are used, but does help avoid
    // "doubling" effects during dynamic permeability development condition 3:
    // take an average of the whole lattice and apply that to boundaries. Needs
    // testing.

    for (j = 1; j < dim - 1; j++) // copy neighbour row
    {
      solid_fraction[0][j] = solid_fraction[1][j];
      vel_x[0][j] = vel_x[1][j];
      vel_y[0][j] = vel_y[1][j];

      solid_fraction[dim - 1][j] = solid_fraction[dim - 2][j];
      vel_x[dim - 1][j] = vel_x[dim - 2][j];
      vel_y[dim - 1][j] = vel_y[dim - 2][j];
    }

    for (i = 0; i < dim; i++) // copy neighbour row
    {
      solid_fraction[i][0] = solid_fraction[i][1];
      vel_x[i][0] = vel_x[i][1];
      vel_y[i][0] = vel_y[i][1];

      solid_fraction[i][dim - 1] = solid_fraction[i][dim - 2];
      vel_x[i][dim - 1] = vel_x[i][dim - 2];
      vel_y[i][dim - 1] = vel_y[i][dim - 2];
    }
  }
  if (boundary == 2) {

    for (j = 1; j < dim - 1; j++) // left and right hand side "wrap" rows
    {
      solid_fraction[0][j] = solid_fraction[dim - 2][j];
      vel_x[0][j] = vel_x[dim - 2][j];
      vel_y[0][j] = vel_y[dim - 2][j];

      solid_fraction[dim - 1][j] = solid_fraction[1][j];
      vel_x[dim - 1][j] = vel_x[1][j];
      vel_y[dim - 1][j] = vel_y[1][j];
    }

    for (i = 0; i < dim; i++) // bottom and top copy neighbour
    {
      solid_fraction[i][0] = solid_fraction[i][1];
      vel_x[i][0] = vel_x[i][1];
      vel_y[i][0] = vel_y[i][1];

      solid_fraction[i][dim - 1] = solid_fraction[i][dim - 2];
      vel_x[i][dim - 1] = vel_x[i][dim - 2];
      vel_y[i][dim - 1] = vel_y[i][dim - 2];
    }
  }
  if (boundary == 3) // use an average
  {
    // at the moment the average works best but needs to be multiplied
    // by a factor, otherwise the boundaries are too permeable
    // the main problem are the side boundaries if the fluid and reaction
    // travel upwards

    for (j = 1; j < dim - 1; j++) // check this boundary condition
    {
      solid_fraction[0][j] = 1.12 * solid_average / count;
      vel_x[0][j] = vel_x_av / count;
      vel_y[0][j] = vel_y_av / count;

      solid_fraction[dim - 1][j] = 1.12 * solid_average / count;
      vel_x[dim - 1][j] = vel_x_av / count;
      vel_y[dim - 1][j] = vel_y_av / count;
    }

    for (i = 0; i < dim; i++) // check this boundary condition
    {
      solid_fraction[i][0] = solid_average / count;
      vel_x[i][0] = vel_x_av / count;
      vel_y[i][0] = vel_y_av / count;

      solid_fraction[i][dim - 1] = solid_average / count;
      vel_x[i][dim - 1] = vel_x_av / count;
      vel_y[i][dim - 1] = vel_y_av / count;
    }
  }
}

//---------------------------------------------------------------------------------
// Caculate permeability from the porosity using Karmen Cozeny relation
// The grain size is now read in directly and changes the advection a lot,
// because of the permeability
//---------------------------------------------------------------------------------

void Fluid_Lattice::Calculate_Permeability_Matrix(double Coseny_grain_size) {
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      kappa[i][j] = 0;

      kappa[i][j] =
          (pow(Coseny_grain_size, 2.0) * pow((1 - solid_fraction[i][j]), 3.0)) /
          (45.0 * pow((solid_fraction[i][j]), 2.0));

      // if (fracture[i][j]>0.0)
      //{
      //	kappa[i][j] += fracture[i][j];
      // }

      if (kappa[i][j] >
          0.0000001) // reset if too large, probably set too small (e11)
      {
        kappa[i][j] = 0.0000001;
        cout << "reset permeability";
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CALCULATION MATRIX INVERSION
// based on implicit alternating direction method, here x and y are the same
// length, so only alpha needed
//
// alpha-set
//
// matrices set 1
//
// multiplication
//
// invert
//
// matrices set 2
//
// multiplication invers
//
// Transpose
//
// source term addition
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation. Note that pressure
// is also in this, meaning that the calucation is easily unstable because it
// varies so much if pressure varies. at the moment the matrix size is preset to
// 200 200 200 for latte simulations with 400 particles in x we are not using
// sparse matrices, this is something that should be changed, would speed up
// things this is the prefactor for each node in the grid where pf, permeability
// and porosity vary this prefactor makes the implicit method potentially
// unstable if its too large, this should be checked and corrected for in a
// future version to make the method more stable
//--------------------------------------------------------------------------------------------------------------------

int Fluid_Lattice::alpha_set(int adjust) // need only alpha
{
  int i, j;
  double max;
  max = 0.0;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      alpha[i][j] = 0;

      {
        // alpha[i][j] = (1.0 + (compressibility *
        // Pf[i][j]))*density_factor[i][j];
        alpha[i][j] = (1.0 + (compressibility * Pf[i][j]));
        alpha[i][j] *= kappa[i][j] * time_a;
        alpha[i][j] /=
            (2.0 * mu * compressibility * (1 - solid_fraction[i][j]) *
             pow(width * scale_constt, 2.0));

        alpha[i][j] = alpha[i][j] / adjust;

        if (alpha[i][j] > max)
          max = alpha[i][j];
      }
    }
  }
  cout << max << endl;

  if (max < 50) // was 50
    return 1;
  else
    return (int)max / 50;
}

//-----------------------------------------------------------------------------------
// set the matrix for the first calculation step
// this is the matrix that could be sparse potentially
//-----------------------------------------------------------------------------------

void Fluid_Lattice::matrices_set_1() {
  // creats L.H.S (al[][][])  and R.H.S (bl[][][])coeffiecient matrix
  // at first half time step of ADI method

  int k, j, i;

  for (k = 0; k < dim; k++) // matrix stack = x or y
  {
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        al[k][i][j] = 0;
        bl[k][i][j] = 0;

        if (i == j) // for diagonal elements
        {
          if (i == 0 || i == dim - 1) // set boundary conditions
          {
            al[k][i][j] = 1.0; // column set for any given column j
            bl[k][i][j] = 1.0; // row set for any given row i
          } else {
            al[k][i][j] = 1 + 2 * alpha[k][i]; // implicit
            bl[k][i][j] = 1 - 2 * alpha[i][k]; // explicit
          }
        } else if ((i == j - 1 && i > 0) || (i == j + 1 && i < (dim - 1))) {
          al[k][i][j] = -alpha[k][i]; // implicit

          bl[k][i][j] = alpha[i][k]; // explicit
        }
      }
    }
  }
}

//---------------------------------------------------------------------------
// this just copies the fluid pressure matrix in order to get the change in
// time per node for the darcy velocities that are then need for the
// advection code
//----------------------------------------------------------------------------

void Fluid_Lattice::copy_matrix() {
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      oldPf[i][j] = Pf[i][j];
    }
  }
}

//----------------------------------------------------------------------------
// multiply diagonal matrix with the fluid pressure matrix
//-----------------------------------------------------------------------------

void Fluid_Lattice::multiplication() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {

      interim_Pf[i][j] = 0;

      for (k = 0; k < dim; k++) {
        interim_Pf[i][j] += bl[j][i][k] * Pf[k][j];
      }
    }
  }
}

//----------------------------------------------------------------------------
// Matrix inversion
//----------------------------------------------------------------------------

void Fluid_Lattice::invert() {
  double factor;
  int g, i, j, k;

  for (g = 0; g < dim; g++) {
    // Creation of identity matrix
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        a_invl[g][i][j] = 0;

        a_invl[g][i][j] = (i == j) ? 1 : 0;
      }
    }

    // Development of upper-triangular matrix from the input matrix
    // along with the subsequent changes in identity matrix.
    for (i = 0; i < dim - 1; i++) {
      for (k = i + 1; k < dim; k++) {
        factor = (-al[g][k][i]) / al[g][i][i];
        if (!isinf(factor) && !isnan(factor)) {
          for (j = 0; j < dim; j++) {
            al[g][k][j] += al[g][i][j] * factor;
            a_invl[g][k][j] += a_invl[g][i][j] * factor;
          }
        }
      }
    }

    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        if (j < i)
          al[g][i][j] = 0;
      }
    }

    // Development of diagonal matrix from developed upper-triangular
    // matrix along with subsequent changes in modified identity matrix.

    for (i = dim - 1; i > 0; i--) {
      for (k = i - 1; k >= 0; k--) {
        factor = (-al[g][k][i]) / al[g][i][i];
        if (!isinf(factor) && !isnan(factor)) {
          for (j = 0; j < dim; j++) {
            al[g][k][j] += al[g][i][j] * factor;
            a_invl[g][k][j] += a_invl[g][i][j] * factor;
          }
        }
      }
    }

    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        if (j > i)
          al[g][i][j] = 0;
      }
    }

    // Conversion of given matrix into identity matrix and on contrary
    // achievement of inverse matrix from the identity matrix.

    for (i = 0; i < dim; i++) {
      if (al[g][i][i] != 0) {
        factor = al[g][i][i];

        for (j = 0; j < dim; j++) {
          al[g][i][j] /= factor;
          a_invl[g][i][j] /= factor;
        }
      }
    }
  }
}

//------------------------------------------------------------------------------------------
// Set the matrix for the second step of the calculation
//------------------------------------------------------------------------------------------

void Fluid_Lattice::matrices_set_2() {
  // creats L.H.S (al[][][])  and R.H.S (bl[][][])coeffiecient matrix
  // at second half time step of ADI method

  int k, i, j;

  for (k = 0; k < dim; k++) {
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        al[k][i][j] = 0;
        bl[k][i][j] = 0;

        if (i == j) // for diagonal elements
        {
          if (i == 0 || i == dim - 1) // set boundary conditions
          {
            al[k][i][j] = 1.0; // column set for any given column j
            bl[k][i][j] = 1.0; // row set for any given row i
          } else {
            al[k][i][j] = 1 + 2 * alpha[i][k];
            bl[k][i][j] = 1 - 2 * alpha[k][i];
          }
        } else if ((i == j - 1 && i > 0) || (i == j + 1 && i < (dim - 1))) {
          al[k][i][j] = -alpha[i][k];

          bl[k][i][j] = alpha[k][i];
        }
      }
    }
  }
}

//--------------------------------------------------------------------------------
// multiply interim matrix with inverted matrix
//--------------------------------------------------------------------------------

void Fluid_Lattice::multiplication_inv() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      Pf[i][j] = 0;

      for (k = 0; k < dim; k++)
        Pf[i][j] += a_invl[j][i][k] * interim_Pf[k][j];
    }
  }
}

//--------------------------------------------------------------------------------------
// Source term that is added to the pressure diffusion equation, uses pressure
// difference and difference in velocity of solid.
//
// check that velocity is scaled (m/second)
//--------------------------------------------------------------------------------------

void Fluid_Lattice::test_pressure_with_source(int byo) {
  int i, j, jj;
  double averageT, averageRho;

  averageT = 0;
  averageRho = 0;
  jj = 0;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      source_a[i][j] = 0;

      if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
        source_a[i][j] = ((vel_x[i][j + 1] - vel_x[i][j - 1]) /
                          (2.0 * width * scale_constt));
        source_a[i][j] += ((vel_y[i + 1][j] - vel_y[i - 1][j]) /
                           (2.0 * width * scale_constt));
        source_a[i][j] *= (1.0 / compressibility + Pf[i][j]);
        source_a[i][j] *= time_a / (2.0 * (1 - solid_fraction[i][j]));
      }
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      source_b[i][j] = 0;

      if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
        source_b[i][j] = vel_x[i][j] * ((Pf[i][j + 1] - Pf[i][j - 1]) /
                                        (2 * width * scale_constt));
        source_b[i][j] += vel_y[i][j] * ((Pf[i + 1][j] - Pf[i - 1][j]) /
                                         (2.0 * width * scale_constt));
        source_b[i][j] *= time_a / 2.0;
      }
    }
  }
  if (byo == 1) {
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        source_c[i][j] = 0;

        if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
          source_c[i][j] = solid_fraction[i][j] * 1600 + 1000;
          source_c[i][j] = source_c[i][j] * kappa[i][j] * time_a /
                           (2 * mu * (1 - solid_fraction[i][j]) *
                            pow(width * scale_constt, 2.0));
          source_c[i][j] = source_c[i][j] * compressibility *
                           (Pf[i][j + 1] - Pf[i][j - 1]) * 9.8;
        }
      }
    }
  } else if (byo == 2) {
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        source_c[i][j] = 0;

        if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
          source_c[i][j] = kappa[i][j] * time_a /
                           (2 * mu * (1 - solid_fraction[i][j]) *
                            pow(width * scale_constt, 2.0));
          source_c[i][j] = source_c[i][j] *
                           ((Pf[i][j + 1] - Pf[i][j - 1]) +
                            ((compressibility_c / compressibility) *
                             (Con[i][j + 1] - Con[i][j - 1])) +
                            ((compressibility_T / compressibility) *
                             (Temperature[i][j + 1] - Temperature[i][j - 1])));
        }
      }
    }
  } else if (byo == 3) {
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        source_c[i][j] = 0;

        if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
          source_c[i][j] = kappa[i][j] * time_a /
                           (2 * mu * (1 - solid_fraction[i][j]) *
                            pow(width * scale_constt, 2.0));
          source_c[i][j] = source_c[i][j] *
                           ((Pf[i][j + 1] - Pf[i][j - 1]) +
                            ((compressibility_c / compressibility) *
                             (Con[i][j + 1] - Con[i][j - 1])) +
                            ((compressibility_T / compressibility) *
                             (Temperature[i][j + 1] - Temperature[i][j - 1])));
        }
      }
    }
  } else if (byo == 4) {
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {

        if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
          averageT += Temperature[i][j];
          averageRho += rho_f[i][j];
          jj++;
        }
      }
    }

    averageT = averageT / jj;
    averageRho = averageRho / jj;
    cout << averageRho << endl;

    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        source_c[i][j] = 0;

        if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
          source_c[i][j] = solid_fraction[i][j] * 1600 + 1000;
          source_c[i][j] = source_c[i][j] * kappa[i][j] * time_a /
                           (2 * mu * (1 - solid_fraction[i][j]) *
                            (width * scale_constt, 2.0));
          // source_c[i][j] =
          // source_c[i][j]*compressibility*(Pf[i][j]-(Pf[i][j+1]+Pf[i][j-1])/2)*9.8
          // +
          // compressibility_c*(Con[i][j]-(Con[i][j+1]+Con[i][j-1])/2)*9.8+1000000000*compressibility_T*((Temperature[i][j+1]-Temperature[i][j-1])/2)*9.8;
          // source_c[i][j] =
          // source_c[i][j]*(compressibility*(Pf[i][j]-(Pf[i][j+1]+Pf[i][j-1])/2)*9.8
          // +
          // compressibility_c*(Con[i][j]-(Con[i][j+1]+Con[i][j-1])/2)*9.8+1000000000000*compressibility_T*((Temperature[i][j+1]-Temperature[i][j-1])/2)*9.8);
          source_c[i][j] =
              source_c[i][j] *
              ((Pf[i][j] - (Pf[i][j + 1] + Pf[i][j - 1]) / 2) * 9.8 +
               10000 * compressibility_c *
                   ((rho_f[i][j + 1] - averageRho) -
                    (rho_f[i][j - 1] - averageRho)) *
                   9.8 +
               10000 * (compressibility_T / compressibility) *
                   ((Temperature[i][j + 1] - averageT) -
                    (Temperature[i][j - 1] - averageT)) *
                   9.8);
          cout << source[i][j] << endl; // 100000000000000
        }
      }
    }
  }
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      source[i][j] = 0;

      source[i][j] = source_a[i][j] + source_b[i][j];
      if (byo > 0)
        source[i][j] = source[i][j] + source_c[i][j];
      // if (source[i][j]!=0.0)
      // cout << source[i][j] << endl;
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      interim_Pf[i][j] -= source[i][j];
      // if (interim_Pf[i][j]>1000)
      // cout << interim_Pf[i][j] << endl;
    }
  }
}

//------------------------------------------------------------------------------------
// Transpose the matrix (switch x and y)
//------------------------------------------------------------------------------------

void Fluid_Lattice::transpose_trans_P() {
  double pressure_trans[dim][dim];
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      pressure_trans[i][j] = 0;

      pressure_trans[i][j] = interim_Pf[j][i];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      interim_Pf[i][j] = 0;

      interim_Pf[i][j] = pressure_trans[i][j];
    }
  }
}

//-----------------------------------------------------------------------------------------------
// concentration calculation Advective part
// at the moment an explicit forward calculation where Ci time n+1 is Ci time n
// minus velocity vector times delta t/delta x times Ci at time n minus Ci-1 at
// time n. Needs to determine the velocity vector first to find direction of
// Ci-1 (backward to flow) need to check stability, and is supposed to induce
// negative diffusion
//
// the prefactor is adjusted internally for the worst case scenario in the code
// to make it stable however, that means every time step is potentially
// different so need to adjust that internally
//------------------------------------------------------------------------------------------------

void Fluid_Lattice::SolveAdvectionLax_Wenddroff(double time, double space,
                                                int boundary, double max) {
  double C_for_x, C_back_x, C_for_y, C_back_y;
  int i, j, k, v, iteration;
  double dx, dy, dt;
  double prefactor, maxprefactor;

  dx = space / 100;
  dy = space / 100;
  dt = time;

  // this boundary condition should probably be applied externally?

  // for (v = 0; v < dim; v++)
  //{
  //	Con[v][0]=0.5; // setting the lower boundary
  // Con[v][dim-1]=0.0;
  //}

  // determine the largest prefactor

  maxprefactor = 0.0;

  for (i = 1; i < dim - 1; i++) // x loop
  {
    for (j = 1; j < dim - 1; j++) // y loop
    {
      prefactor = fabs(vel_darcyx[i][j] * dt / dx);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
      prefactor = fabs(vel_darcyy[i][j] * dt / dy);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
    }
  }

  cout << "advection prefactor:" << maxprefactor << endl;

  // adjust prefactor to be 0.2, this should be an internal loop
  // so that the time is constant. Theoretically should be stable
  // already when prefactor is 1. So advection loop should not loop
  // more than 50 times.

  iteration = maxprefactor / 0.02;

  if (iteration < 1)
    iteration = 1;

  cout << "iteration advection:" << iteration << endl;

  // step 1 calculation of half space back and forward for x and y.

  // step 2
  for (k = 0; k < iteration; k++) {
    for (i = 0; i < dim; i++) // x loop
    {
      for (j = 1; j < dim - 1; j++) // y loop
      {
        if (i == 0) // boundary conditions are letting fluid flow in y parallel
                    // to boundary, x term constant
        {
          if (boundary == 1) {
            C_for_x =
                0.5 * (Con[i][j] + Con[i + 1][j]) -
                0.5 * vel_darcyx[i][j] * dt / dx * (Con[i + 1][j] - Con[i][j]);
            C_back_x = 0.5 * (Con[i][j] + Con[dim - 1][j]) +
                       0.5 * vel_darcyx[i][j] * dt / dx *
                           (Con[dim - 1][j] - Con[i][j]);

            // step 2 x
            dxc[i][j] =
                Con[i][j] - vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
          }
          if (boundary == 2) //
            dxc[i][j] = Con[i][j];
        } else if (i == dim - 1) // wrapping right boundary
        {
          if (boundary == 1) {
            C_for_x =
                0.5 * (Con[i][j] + Con[0][j]) -
                0.5 * vel_darcyx[i][j] * dt / dx * (Con[0][j] - Con[i][j]);
            C_back_x =
                0.5 * (Con[i][j] + Con[i - 1][j]) +
                0.5 * vel_darcyx[i][j] * dt / dx * (Con[i - 1][j] - Con[i][j]);

            // step 2 x
            dxc[i][j] =
                Con[i][j] - vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
          }
          if (boundary == 2) //
            dxc[i][j] = Con[i][j];
        } else {

          // step 1 x
          C_for_x =
              0.5 * (Con[i][j] + Con[i + 1][j]) -
              0.5 * vel_darcyx[i][j] * dt / dx * (Con[i + 1][j] - Con[i][j]);
          C_back_x =
              0.5 * (Con[i][j] + Con[i - 1][j]) +
              0.5 * vel_darcyx[i][j] * dt / dx * (Con[i - 1][j] - Con[i][j]);

          // step 2 x
          dxc[i][j] =
              Con[i][j] - vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
        }

        // step 1 y
        C_for_y =
            0.5 * (Con[i][j] + Con[i][j + 1]) -
            0.5 * vel_darcyy[i][j] * dt / dy * (Con[i][j + 1] - Con[i][j]);
        C_back_y =
            0.5 * (Con[i][j] + Con[i][j - 1]) +
            0.5 * vel_darcyy[i][j] * dt / dy * (Con[i][j - 1] - Con[i][j]);

        // step 2 y
        dyc[i][j] =
            Con[i][j] - vel_darcyy[i][j] * dt / dy * (C_for_y - C_back_y);
      }
    }
    for (i = 0; i < dim; i++) // x loop
    {
      for (j = 1; j < dim - 1; j++) // y loop
      {
        // add up x and y

        Con[i][j] = 0.5 * (dxc[i][j] + dyc[i][j]);

        // final adjustment, no overflow and no negative concentration to
        // improve stability

        if (Con[i][j] > max)
          Con[i][j] = max;
        if (Con[i][j] < 0.0)
          Con[i][j] = 0.0;
      }
    }
  }
}

//-----------------------------------------------------------------------------------------------
// concentration calculation Advective part
// at the moment an explicit forward calculation where Ci time n+1 is Ci time n
// minus velocity vector times delta t/delta x times Ci at time n minus Ci-1 at
// time n. Needs to determine the velocity vector first to find direction of
// Ci-1 (backward to flow) need to check stability, and is supposed to induce
// negative diffusion
//
// the prefactor is adjusted internally for the worst case scenario in the code
// to make it stable however, that means every time step is potentially
// different so need to adjust that internally
//------------------------------------------------------------------------------------------------

void Fluid_Lattice::SolveAdvectionLax_Wenddroff_Phreeqc(double time,
                                                        double space,
                                                        int boundary) {
  double C_for_x, C_back_x, C_for_y, C_back_y;
  int i, j, k, v, iteration;
  double dx, dy, dt;
  double prefactor, maxprefactor;

  dx = space / 100;
  dy = space / 100;
  dt = time;

  // this boundary condition should probably be applied externally?

  // for (v = 0; v < dim; v++)
  //{
  //	Con[v][0]=0.5; // setting the lower boundary
  // Con[v][dim-1]=0.0;
  //}

  // determine the largest prefactor

  maxprefactor = 0.0;

  for (i = 1; i < dim - 1; i++) // x loop
  {
    for (j = 1; j < dim - 1; j++) // y loop
    {
      prefactor = fabs(vel_darcyx[i][j] * dt / dx);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
      prefactor = fabs(vel_darcyy[i][j] * dt / dy);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
    }
  }

  cout << "advection prefactor:" << maxprefactor << endl;

  // adjust prefactor to be 0.2, this should be an internal loop
  // so that the time is constant.

  iteration = maxprefactor / 0.02;

  if (iteration < 1)
    iteration = 1;

  cout << "iteration advection:" << iteration << endl;

  // step 1 calculation of half space back and forward for x and y.

  // step 2
  for (k = 0; k < iteration; k++) {
    for (i = 0; i < dim; i++) // x loop
    {
      for (j = 1; j < dim - 1; j++) // y loop
      {
        if (i == 0) // boundary conditions are letting fluid flow in y parallel
                    // to boundary, x term constant
        {
          if (boundary == 1) {
            C_for_x = 0.5 * (Temperature[i][j] + Temperature[i + 1][j]) -
                      0.5 * vel_darcyx[i][j] * dt / dx *
                          (Temperature[i + 1][j] - Temperature[i][j]);
            C_back_x = 0.5 * (Temperature[i][j] + Temperature[dim - 1][j]) +
                       0.5 * vel_darcyx[i][j] * dt / dx *
                           (Temperature[dim - 1][j] - Temperature[i][j]);

            // step 2 x
            dxc[i][j] = Temperature[i][j] -
                        vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
          }
          if (boundary == 2) //
            dxc[i][j] = Temperature[i][j];
        } else if (i == dim - 1) // wrapping right boundary
        {
          if (boundary == 1) {
            C_for_x = 0.5 * (Temperature[i][j] + Temperature[0][j]) -
                      0.5 * vel_darcyx[i][j] * dt / dx *
                          (Temperature[0][j] - Temperature[i][j]);
            C_back_x = 0.5 * (Temperature[i][j] + Temperature[i - 1][j]) +
                       0.5 * vel_darcyx[i][j] * dt / dx *
                           (Temperature[i - 1][j] - Temperature[i][j]);

            // step 2 x
            dxc[i][j] = Temperature[i][j] -
                        vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
          }
          if (boundary == 2) //
            dxc[i][j] = Temperature[i][j];
        } else {

          // step 1 x
          C_for_x = 0.5 * (Temperature[i][j] + Temperature[i + 1][j]) -
                    0.5 * vel_darcyx[i][j] * dt / dx *
                        (Temperature[i + 1][j] - Temperature[i][j]);
          C_back_x = 0.5 * (Temperature[i][j] + Temperature[i - 1][j]) +
                     0.5 * vel_darcyx[i][j] * dt / dx *
                         (Temperature[i - 1][j] - Temperature[i][j]);

          // step 2 x
          dxc[i][j] = Temperature[i][j] -
                      vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
        }

        // step 1 y
        C_for_y = 0.5 * (Temperature[i][j] + Temperature[i][j + 1]) -
                  0.5 * vel_darcyy[i][j] * dt / dy *
                      (Temperature[i][j + 1] - Temperature[i][j]);
        C_back_y = 0.5 * (Temperature[i][j] + Temperature[i][j - 1]) +
                   0.5 * vel_darcyy[i][j] * dt / dy *
                       (Temperature[i][j - 1] - Temperature[i][j]);

        // step 2 y
        dyc[i][j] = Temperature[i][j] -
                    vel_darcyy[i][j] * dt / dy * (C_for_y - C_back_y);
      }
    }
    for (i = 0; i < dim; i++) // x loop
    {
      for (j = 1; j < dim - 1; j++) // y loop
      {
        // add up x and y

        Temperature[i][j] = 0.5 * (dxc[i][j] + dyc[i][j]);

        // cout << conc[i][j] << endl;

        // final adjustment, no overflow and no negative concentration to
        // improve stability - needs to be adjusted to Phreeqc

        if (Temperature[i][j] > 100000000)
          Temperature[i][j] = 100000000;
        if (Temperature[i][j] < 0.0)
          Temperature[i][j] = 0.0;
      }
    }
  }
}

//-----------------------------------------------------------------------------------------------
// concentration calculation Advective part
// at the moment an explicit forward calculation where Ci time n+1 is Ci time n
// minus velocity vector times delta t/delta x times Ci at time n minus Ci-1 at
// time n. Needs to determine the velocity vector first to find direction of
// Ci-1 (backward to flow) need to check stability, and is supposed to induce
// negative diffusion
//
// the prefactor is adjusted internally for the worst case scenario in the code
// to make it stable however, that means every time step is potentially
// different so need to adjust that internally
//------------------------------------------------------------------------------------------------

void Fluid_Lattice::SolveAdvectionLax_Wenddroff_Metal(double time, double space,
                                                      int boundary) {
  double C_for_x, C_back_x, C_for_y, C_back_y;
  int i, j, k, v, iteration;
  double dx, dy, dt;
  double prefactor, maxprefactor;

  dx = space / 100;
  dy = space / 100;
  dt = time;

  // this boundary condition should probably be applied externally?

  // for (v = 0; v < dim; v++)
  //{
  //	Con[v][0]=0.5; // setting the lower boundary
  // Con[v][dim-1]=0.0;
  //}

  // determine the largest prefactor

  maxprefactor = 0.0;

  for (i = 1; i < dim - 1; i++) // x loop
  {
    for (j = 1; j < dim - 1; j++) // y loop
    {
      prefactor = fabs(vel_darcyx[i][j] * dt / dx);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
      prefactor = fabs(vel_darcyy[i][j] * dt / dy);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
    }
  }

  cout << "advection prefactor:" << maxprefactor << endl;

  // adjust prefactor to be 0.2, this should be an internal loop
  // so that the time is constant.

  iteration = maxprefactor / 0.02;

  if (iteration < 1)
    iteration = 1;

  cout << "iteration advection:" << iteration << endl;

  // step 1 calculation of half space back and forward for x and y.

  // step 2
  for (k = 0; k < iteration; k++) {
    for (i = 0; i < dim; i++) // x loop
    {
      for (j = 1; j < dim - 1; j++) // y loop
      {
        if (i == 0) // boundary conditions are letting fluid flow in y parallel
                    // to boundary, x term constant
        {
          if (boundary == 1) {
            C_for_x = 0.5 * (Metal_A[i][j] + Metal_A[i + 1][j]) -
                      0.5 * vel_darcyx[i][j] * dt / dx *
                          (Metal_A[i + 1][j] - Metal_A[i][j]);
            C_back_x = 0.5 * (Metal_A[i][j] + Metal_A[dim - 1][j]) +
                       0.5 * vel_darcyx[i][j] * dt / dx *
                           (Metal_A[dim - 1][j] - Metal_A[i][j]);

            // step 2 x
            dxc[i][j] = Metal_A[i][j] -
                        vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
          }
          if (boundary == 2) //
            dxc[i][j] = Metal_A[i][j];
        } else if (i == dim - 1) // wrapping right boundary
        {
          if (boundary == 1) {
            C_for_x = 0.5 * (Metal_A[i][j] + Metal_A[0][j]) -
                      0.5 * vel_darcyx[i][j] * dt / dx *
                          (Metal_A[0][j] - Metal_A[i][j]);
            C_back_x = 0.5 * (Metal_A[i][j] + Metal_A[i - 1][j]) +
                       0.5 * vel_darcyx[i][j] * dt / dx *
                           (Metal_A[i - 1][j] - Metal_A[i][j]);

            // step 2 x
            dxc[i][j] = Metal_A[i][j] -
                        vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
          }
          if (boundary == 2) //
            dxc[i][j] = Metal_A[i][j];
        } else {

          // step 1 x
          C_for_x = 0.5 * (Metal_A[i][j] + Metal_A[i + 1][j]) -
                    0.5 * vel_darcyx[i][j] * dt / dx *
                        (Metal_A[i + 1][j] - Metal_A[i][j]);
          C_back_x = 0.5 * (Metal_A[i][j] + Metal_A[i - 1][j]) +
                     0.5 * vel_darcyx[i][j] * dt / dx *
                         (Metal_A[i - 1][j] - Metal_A[i][j]);

          // step 2 x
          dxc[i][j] =
              Metal_A[i][j] - vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
        }

        // step 1 y
        C_for_y = 0.5 * (Metal_A[i][j] + Metal_A[i][j + 1]) -
                  0.5 * vel_darcyy[i][j] * dt / dy *
                      (Metal_A[i][j + 1] - Metal_A[i][j]);
        C_back_y = 0.5 * (Metal_A[i][j] + Metal_A[i][j - 1]) +
                   0.5 * vel_darcyy[i][j] * dt / dy *
                       (Metal_A[i][j - 1] - Metal_A[i][j]);

        // step 2 y
        dyc[i][j] =
            Metal_A[i][j] - vel_darcyy[i][j] * dt / dy * (C_for_y - C_back_y);
      }
    }
    for (i = 0; i < dim; i++) // x loop
    {
      for (j = 1; j < dim - 1; j++) // y loop
      {
        // add up x and y

        Metal_A[i][j] = 0.5 * (dxc[i][j] + dyc[i][j]);

        // final adjustment, no overflow and no negative concentration to
        // improve stability

        if (Metal_A[i][j] > 0.5)
          Metal_A[i][j] = 0.5;
        if (Metal_A[i][j] < 0.0)
          Metal_A[i][j] = 0.0;
      }
    }
  }
}

//-----------------------------------------------------------------------------------------------
// concentration calculation Advective part
// at the moment an explicit forward calculation where Ci time n+1 is Ci time n
// minus velocity vector times delta t/delta x times Ci at time n minus Ci-1 at
// time n. Needs to determine the velocity vector first to find direction of
// Ci-1 (backward to flow) need to check stability, and is supposed to induce
// negative diffusion
//
// the prefactor is adjusted internally for the worst case scenario in the code
// to make it stable however, that means every time step is potentially
// different so need to adjust that internally
//------------------------------------------------------------------------------------------------

void Fluid_Lattice::SolveAdvectionLax_Wenddroff_Temperature(double time,
                                                            double space,
                                                            int boundary) {
  double C_for_x, C_back_x, C_for_y, C_back_y;
  int i, j, k, v, iteration;
  double dx, dy, dt;
  double prefactor, maxprefactor;

  dx = space / 100;
  dy = space / 100;
  dt = time;

  // this boundary condition should probably be applied externally?

  // for (v = 0; v < dim; v++)
  //{
  //	Con[v][0]=0.5; // setting the lower boundary
  // Con[v][dim-1]=0.0;
  //}

  // determine the largest prefactor

  maxprefactor = 0.0;

  for (i = 1; i < dim - 1; i++) // x loop
  {
    for (j = 1; j < dim - 1; j++) // y loop
    {
      prefactor = fabs(vel_darcyx[i][j] * dt / dx);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
      prefactor = fabs(vel_darcyy[i][j] * dt / dy);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
    }
  }

  cout << "advection prefactor:" << maxprefactor << endl;

  // adjust prefactor to be 0.2, this should be an internal loop
  // so that the time is constant.

  iteration = maxprefactor / 0.02;

  if (iteration < 1)
    iteration = 1;

  cout << "iteration advection:" << iteration << endl;

  // step 1 calculation of half space back and forward for x and y.

  // step 2
  for (k = 0; k < iteration; k++) {
    for (i = 0; i < dim; i++) // x loop
    {
      for (j = 1; j < dim - 1; j++) // y loop
      {
        if (i == 0) // boundary conditions are letting fluid flow in y parallel
                    // to boundary, x term constant
        {
          if (boundary == 1) {
            C_for_x = 0.5 * (Temperature[i][j] + Temperature[i + 1][j]) -
                      0.5 * vel_darcyx[i][j] * dt / dx *
                          (Temperature[i + 1][j] - Temperature[i][j]);
            C_back_x = 0.5 * (Temperature[i][j] + Temperature[dim - 1][j]) +
                       0.5 * vel_darcyx[i][j] * dt / dx *
                           (Temperature[dim - 1][j] - Temperature[i][j]);

            // step 2 x
            dxc[i][j] = Temperature[i][j] -
                        vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
          }
          if (boundary == 2) //
            dxc[i][j] = Temperature[i][j];
        } else if (i == dim - 1) // wrapping right boundary
        {
          if (boundary == 1) {
            C_for_x = 0.5 * (Temperature[i][j] + Temperature[0][j]) -
                      0.5 * vel_darcyx[i][j] * dt / dx *
                          (Temperature[0][j] - Temperature[i][j]);
            C_back_x = 0.5 * (Temperature[i][j] + Temperature[i - 1][j]) +
                       0.5 * vel_darcyx[i][j] * dt / dx *
                           (Temperature[i - 1][j] - Temperature[i][j]);

            // step 2 x
            dxc[i][j] = Temperature[i][j] -
                        vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
          }
          if (boundary == 2) //
            dxc[i][j] = Temperature[i][j];
        } else {

          // step 1 x
          C_for_x = 0.5 * (Temperature[i][j] + Temperature[i + 1][j]) -
                    0.5 * vel_darcyx[i][j] * dt / dx *
                        (Temperature[i + 1][j] - Temperature[i][j]);
          C_back_x = 0.5 * (Temperature[i][j] + Temperature[i - 1][j]) +
                     0.5 * vel_darcyx[i][j] * dt / dx *
                         (Temperature[i - 1][j] - Temperature[i][j]);

          // step 2 x
          dxc[i][j] = Temperature[i][j] -
                      vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
        }

        // step 1 y
        C_for_y = 0.5 * (Temperature[i][j] + Temperature[i][j + 1]) -
                  0.5 * vel_darcyy[i][j] * dt / dy *
                      (Temperature[i][j + 1] - Temperature[i][j]);
        C_back_y = 0.5 * (Temperature[i][j] + Temperature[i][j - 1]) +
                   0.5 * vel_darcyy[i][j] * dt / dy *
                       (Temperature[i][j - 1] - Temperature[i][j]);

        // step 2 y
        dyc[i][j] = Temperature[i][j] -
                    vel_darcyy[i][j] * dt / dy * (C_for_y - C_back_y);
      }
    }
    for (i = 0; i < dim; i++) // x loop
    {
      for (j = 1; j < dim - 1; j++) // y loop
      {
        // add up x and y

        Temperature[i][j] = 0.5 * (dxc[i][j] + dyc[i][j]);

        // final adjustment, no overflow and no negative concentration to
        // improve stability

        if (Temperature[i][j] > 85)
          Temperature[i][j] = 85;
        if (Temperature[i][j] < 0.0)
          Temperature[i][j] = 0.0;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------------
// this is not used at the moment, but can potentially be used to change
// concentration after a reaction
//
// sept 2018
//--------------------------------------------------------------------------------------------------------

void Fluid_Lattice::ChangeConcentration(int box_x, int box_y, double change) {
  Con[box_y][box_x] = Con[box_y][box_x] - change;
  if (Con[box_y][box_x] < 0.0)
    Con[box_y][box_x] = 0.0;
}

//------------------------------------------------------------------------------------------
// Passing back the values to the solid lattice including the fluid pressure
// gradients as force in x and y that is added to the force ballance in
// relaxation, but not to the calculated stress, so that the stress remains to
// be the solid stress only.
//
// This takes each pressure node and looks into its 4 repulsion boxes and gives
// the values to the particles in there.
//
// can potentially calculate average concetration between concentration nodes,
// but this does not seem to work well at the moment sept 2018
//
// The pressure gradient is taken from the particle repulsion box to the
// neighbouring node
//-------------------------------------------------------------------------------------------

void Fluid_Lattice::Pass_Back_Gradients(Particle **list, int av_conc) {
  int i, j, k;
  Particle *help;
  double smooth, grad_x, grad_y, conc_av_x, conc_av_y;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      //--------------------------------------------------------------------
      // 0 is lower left box, 1 lower right box, 2 upper left box
      // and 3 upper right box.
      //
      // we need special boundary conditions for the boundary boxes
      // left and right can be wrappinp, but up and down may be more
      // problematic
      //
      // at the moment right and left wrap, up and down copy next row
      //--------------------------------------------------------------------

      for (k = 0; k < 4; k++) // go to four quadrants of pressure node
      {
        if (k == 0)
          box =
              (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 + ((i - 1) * 2);
        else if (k == 1)
          box =
              (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 + ((i - 1) * 2);
        else if (k == 2)
          box =
              (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 + ((i - 1) * 2);
        else if (k == 3)
          box =
              (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 + ((i - 1) * 2);

        if (list[box]) {
          help = list[box];
          x_diff = abs(help->xpos - ((double(i) * width) + (width / 2.0)));
          y_diff = abs(help->ypos - ((double(j) * width) + (width / 2.0)));
          smooth = (1 - (x_diff / width)) * (1 - (y_diff / width));

          if (k == 0) // lower left box
          {
            if (i == 0) // left boundary, deal with x
            {
              grad_x = Pf[dim - 1][j] - Pf[i][j]; // wrap
              conc_av_x = (Con[dim - 1][j] + Con[i][j]) / 2.0;
            } else // normal
            {
              grad_x = Pf[i - 1][j] - Pf[i][j];
              conc_av_x = (Con[dim - 1][j] + Con[i][j]) / 2.0;
            }
            if (j == 0) // lower boundary, deal with y
            {
              grad_y =
                  Pf[i][j] - Pf[i][j + 1]; // use the next row and copy that
              conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
            } else // normal
            {
              grad_y = Pf[i][j - 1] - Pf[i][j];
              conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
            }
          } else if (k == 1) // lower right box
          {
            if (i == dim - 1) // wrap
            {
              grad_x = Pf[i][j] - Pf[0][j];
              conc_av_x = (Con[i][j] + Con[0][j]) / 2.0;
            } else // normal
            {
              grad_x = Pf[i][j] - Pf[i + 1][j];
              conc_av_x = (Con[i][j] + Con[i + 1][j]) / 2.0;
            }
            if (j == 0) {
              grad_y = Pf[i][j] - Pf[i][j + 1];
              conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
            } else {
              grad_y = Pf[i][j - 1] - Pf[i][j];
              conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
            }
          } else if (k == 2) {
            if (i == 0) // left boundary
            {
              grad_x = Pf[dim - 1][j] - Pf[i][j];
              conc_av_x = (Con[dim - 1][j] + Con[i][j]) / 2.0;
            } else {
              grad_x = Pf[i - 1][j] - Pf[i][j];
              conc_av_x = (Con[i - 1][j] + Con[i][j]) / 2.0;
            }
            if (j == dim - 1) {
              grad_y = Pf[i][j - 1] - Pf[i][j];
              conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
            } else {
              grad_y = Pf[i][j] - Pf[i][j + 1];
              conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
            }

          } else if (k == 3) {
            if (i == dim - 1) {
              grad_x = Pf[i][j] - Pf[0][j];
              conc_av_x = (Con[i][j] + Con[0][j]) / 2.0;
            } else {
              grad_x = Pf[i][j] - Pf[i + 1][j];
              conc_av_x = (Con[i][j] + Con[i + 1][j]) / 2.0;
            }
            if (j == dim - 1) {
              grad_y = Pf[i][j - 1] - Pf[i][j];
              conc_av_x = (Con[i][j - 1] + Con[i][j]) / 2.0;
            } else {
              grad_y = Pf[i][j] - Pf[i][j + 1];
              conc_av_x = (Con[i][j] + Con[i][j + 1]) / 2.0;
            }
          }

          // pass the pressure gradient to each particle as force in x and y

          help->F_P_y = grad_y * 3.14 * help->real_radius * help->real_radius /
                        solid_fraction[i][j];
          help->F_P_x = grad_x * 3.14 * help->real_radius * help->real_radius /
                        solid_fraction[i][j];
          help->F_P_y = help->F_P_y / (width * scale_constt);
          help->F_P_x = help->F_P_x / (width * scale_constt);

          // pass the darcy velocities in x and y. note that they are calculated
          // in a differnt function

          help->v_darcyx = vel_darcyx[i][j];
          help->v_darcyy = vel_darcyy[i][j];

          // pass aditional values for visualization purposes to the particles

          help->fluid_P = Pf[i][j];
          help->fluid_pressure_gradient = help->F_P_y;
          help->average_porosity = 1 - solid_fraction[i][j];

          // use or dont use the average concentration

          if (av_conc)
            help->conc = (conc_av_x + conc_av_y) / 2.0;
          else
            help->conc = Con_S6[i][j];

          // help->metal_conc = rho_f[i][j];
          help->metal_conc = Con_Zn[i][j];

          help->temperature = Temperature[i][j];
          help->calcite = Si[i][j];
          help->IronOxide = Si3[i][j];
          help->Goethite = Si4[i][j];
          help->Jarosite = Si5[i][j];
          help->Talc = Si6[i][j];
          help->Hamatite = Si7[i][j];

          while (help->next_inBox) // and do the whole thing for other particles
                                   // in the box
          {
            help = help->next_inBox;

            x_diff = abs(help->xpos - ((double(i) * width) + (width / 2.0)));
            y_diff = abs(help->ypos - ((double(j) * width) + (width / 2.0)));
            smooth = (1 - (x_diff / width)) * (1 - (y_diff / width));

            if (k == 0) // lower left box
            {
              if (i == 0) // left boundary, deal with x
              {
                grad_x = Pf[dim - 1][j] - Pf[i][j]; // wrap
                conc_av_x = (Con[dim - 1][j] + Con[i][j]) / 2.0;
              } else // normal
              {
                grad_x = Pf[i - 1][j] - Pf[i][j];
                conc_av_x = (Con[i - 1][j] + Con[i][j]) / 2.0;
              }
              if (j == 0) // lower boundary, deal with y
              {
                grad_y =
                    Pf[i][j] - Pf[i][j + 1]; // use the next row and copy that
                conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
              } else // normal
              {
                grad_y = Pf[i][j - 1] - Pf[i][j];
                conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
              }
            } else if (k == 1) // lower right box
            {
              if (i == dim - 1) // wrap
              {
                grad_x = Pf[i][j] - Pf[0][j];
                conc_av_x = (Con[i][j] + Con[0][j]) / 2.0;
              } else // normal
              {
                grad_x = Pf[i][j] - Pf[i + 1][j];
                conc_av_x = (Con[i][j] + Con[i + 1][j]) / 2.0;
              }
              if (j == 0) {
                grad_y = Pf[i][j] - Pf[i][j + 1];
                conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
              } else {
                grad_y = Pf[i][j - 1] - Pf[i][j];
                conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
              }
            } else if (k == 2) {
              if (i == 0) // left boundary
              {
                grad_x = Pf[dim - 1][j] - Pf[i][j];
                conc_av_x = (Con[dim - 1][j] + Con[i][j]) / 2.0;
              } else {
                grad_x = Pf[i - 1][j] - Pf[i][j];
                conc_av_x = (Con[i - 1][j] + Con[i][j]) / 2.0;
              }
              if (j == dim - 1) {
                grad_y = Pf[i][j - 1] - Pf[i][j];
                conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
              } else {
                grad_y = Pf[i][j] - Pf[i][j + 1];
                conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
              }

            } else if (k == 3) {
              if (i == dim - 1) {
                grad_x = Pf[i][j] - Pf[0][j];
                conc_av_x = (Con[i][j] + Con[0][j]) / 2.0;
              } else {
                grad_x = Pf[i][j] - Pf[i + 1][j];
                conc_av_x = (Con[i][j] + Con[i + 1][j]) / 2.0;
              }
              if (j == dim - 1) {
                grad_y = Pf[i][j - 1] - Pf[i][j];
                conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
              } else {
                grad_y = Pf[i][j] - Pf[i][j + 1];
                conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
              }
            }

            // fluid forces from pressure gradients for lattice

            help->F_P_y = grad_y * 3.14 * help->real_radius *
                          help->real_radius / solid_fraction[i][j];
            help->F_P_x = grad_x * 3.14 * help->real_radius *
                          help->real_radius / solid_fraction[i][j];
            help->F_P_y = help->F_P_y / (width * scale_constt);
            help->F_P_x = help->F_P_x / (width * scale_constt);

            // output values for visualization

            help->v_darcyx = vel_darcyx[i][j];
            help->v_darcyy = vel_darcyy[i][j];
            help->fluid_P = Pf[i][j];
            help->fluid_pressure_gradient = help->F_P_y;
            help->average_porosity = 1 - solid_fraction[i][j];

            if (av_conc)
              help->conc = (conc_av_x + conc_av_y) / 2.0;
            else
              help->conc = Con_S6[i][j];

            // help->metal_conc = rho_f[i][j];
            help->metal_conc = Con_Zn[i][j];

            help->IronOxide = Si3[i][j];
            help->Goethite = Si4[i][j];
            help->Jarosite = Si5[i][j];
            help->Talc = Si6[i][j];
            help->Hamatite = Si7[i][j];
            help->temperature = Temperature[i][j];
            help->calcite = Si[i][j];
          }
        }
      }
    }
  }
}

//--------------------------------------------------------------------------------------------------------------
// This just sets some boundaries for the fluid lattice plus potentially one for
// the sides, but that is not necessarily need because the IDE method using
// parallel flow along the boundaries already
//--------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::Set_Boundaries(double UpperBoundaryPressure,
                                   double LowerBoundaryPressure,
                                   int double_sides, double conc,
                                   double additional) {
  int i, j;

  for (i = 0; i < dim; i++) {
    if (UpperBoundaryPressure != 0.0) {
      Pf[i][dim - 1] = UpperBoundaryPressure;
      // Con[i][dim-1] = conc;
    }

    if (LowerBoundaryPressure != 0.0)

    {
      if (double_sides) {
        if (i > 40 & i < 60) {
          Pf[i][0] = LowerBoundaryPressure + additional;
        } else
          Pf[i][0] = LowerBoundaryPressure;
      } else
        Pf[i][0] = Pf[i][1];
    }
  }
}

//-------------------------------------------------------------------------------------------
// simple increase of pressure at a single node
//-------------------------------------------------------------------------------------------

void Fluid_Lattice::Fluid_Input(double pressure, int x, int y) {
  Pf[x][y] = Pf[x][y] + pressure;
  Con[x][y] += 0.5;
}

//-------------------------------------------------------------------------------------------
// this is the core for the advection code, gets the darcy velocity
// at the moment a function of changes in pressure in the pressure diffusion
// code. that means the concentration travels as a function of pressure changes
// for the advection
//
// first get delta P and then the unit vectors for old grad in x and y for the
// direction of flow.
//
// This could potentially be adjusted but works well at the moment sept 2018
//-------------------------------------------------------------------------------------------

void Fluid_Lattice::Get_Fluid_Velocity(double time) {
  int i, j;
  double grad_xold, grad_yold, unit_vector, grad_xold_in, grad_xold_out,
      grad_yold_in, grad_yold_out, u_vec_in, u_vec_out;
  double deltaPf;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      // get delta P, Pf change per time step

      deltaPf = Pf[i][j] - oldPf[i][j];

      if (i == 0) // flow parallel to left and right hand boundaries
      {
        if (j == 0) {
          grad_xold = 0.0;
          grad_yold = 1.0;
          grad_xold_in = 0.0;
          grad_yold_in = 1.0;
          grad_xold_out = 0.0;
          grad_yold_out = 1.0;

        } else if (j == dim - 1) {
          grad_xold = 0.0;
          grad_yold = 1.0;
          grad_xold_in = 0.0;
          grad_yold_in = 1.0;
          grad_xold_out = 0.0;
          grad_yold_out = 1.0;
        } else {
          grad_xold = 0.0;
          grad_yold = oldPf[i][j - 1] - oldPf[i][j + 1];
          grad_xold_in = 0.0;
          grad_yold_in = oldPf[i][j - 1] - oldPf[i][j];
          grad_xold_out = 0.0;
          grad_yold_out = oldPf[i][j] - oldPf[i][j + 1];
        }
      } else if (i == dim - 1) // right boundary wrap
      {
        if (j == 0) {
          grad_xold = 0.0;
          grad_yold = 1.0;
          grad_xold_in = 0.0;
          grad_yold_in = 1.0;
          grad_xold_out = 0.0;
          grad_yold_out = 1.0;
        } else if (j == dim - 1) {
          grad_xold = 0.0;
          grad_yold = 1.0;
          grad_xold_in = 0.0;
          grad_yold_in = 1.0;
          grad_xold_out = 0.0;
          grad_yold_out = 1.0;
        } else {
          grad_xold = 0.0;
          grad_yold = oldPf[i][j - 1] - oldPf[i][j + 1];
          grad_xold_in = 0.0;
          grad_yold_in = oldPf[i][j - 1] - oldPf[i][j];
          grad_xold_out = 0.0;
          grad_yold_out = oldPf[i][j] - oldPf[i][j + 1];
        }
      } else if (j == 0) // bottom static only upwards
      {
        grad_xold = 0.0;
        grad_yold = 1.0;
        grad_xold_in = 0.0;
        grad_yold_in = 1.0;
        grad_xold_out = 0.0;
        grad_yold_out = 1.0;
      } else if (j == dim - 1) // top static only upwards
      {
        grad_xold = 0.0;
        grad_yold = 1.0;
        grad_xold_in = 0.0;
        grad_yold_in = 1.0;
        grad_xold_out = 0.0;
        grad_yold_out = 1.0;
      } else // normal case
      {
        grad_xold = oldPf[i - 1][j] - oldPf[i + 1][j];
        grad_yold = oldPf[i][j - 1] - oldPf[i][j + 1];
        grad_xold_in = oldPf[i - 1][j] - oldPf[i][j];
        grad_yold_in = oldPf[i][j - 1] - oldPf[i][j];
        grad_xold_out = oldPf[i][j] - oldPf[i + 1][j];
        grad_yold_out = oldPf[i][j] - oldPf[i][j + 1];
      }
      // get the unit vectors for x and y

      unit_vector = sqrt((grad_xold * grad_xold) + (grad_yold * grad_yold));
      u_vec_in =
          sqrt((grad_xold_in * grad_xold_in) + (grad_yold_in * grad_yold_in));
      u_vec_out = sqrt((grad_xold_out * grad_xold_out) +
                       (grad_yold_out * grad_yold_out));

      if (unit_vector != 0.0) {
        grad_xold = grad_xold / unit_vector;
        grad_yold = grad_yold / unit_vector;
      } else {
        grad_xold = 0.0;
        grad_yold = 0.0;
      }

      if (u_vec_in != 0.0) {
        grad_xold_in = grad_xold_in / u_vec_in;
        grad_yold_in = grad_yold_in / u_vec_in;
      } else {
        grad_xold_in = 0.0;
        grad_yold_in = 0.0;
      }
      if (u_vec_out != 0.0) {
        grad_xold_out = grad_xold_out / u_vec_out;
        grad_yold_out = grad_yold_out / u_vec_out;
      } else {
        grad_xold_out = 0.0;
        grad_yold_out = 0.0;
      }
      // and finally calculate the darcy velocity in x and y

      vel_darcyx[i][j] = (((grad_xold * deltaPf) / (scale_constt * width)) *
                          (kappa[i][j] / mu)) /
                         ((1 - solid_fraction[i][j]) * time);
      vel_darcyy[i][j] = (((grad_yold * deltaPf) / (scale_constt * width)) *
                          (kappa[i][j] / mu)) /
                         ((1 - solid_fraction[i][j]) * time);

      // vel_darcyy[i][j] = -(((Pf[i][j+1] - Pf[i][j-1])/
      // (2*scale_constt*width))*
      // (kappa[i][j]/mu))/((1-solid_fraction[i][j])*time); vel_darcyx[i][j] =
      // -(((Pf[i+1][j] - Pf[i-1][j])/ (2*scale_constt*width))*
      // (kappa[i][j]/mu))/((1-solid_fraction[i][j])*time);

      vel_darcyy[i][j] =
          0.5 * (vel_darcyy[i][j] -
                 (((Pf[i][j + 1] - Pf[i][j - 1]) / (2 * scale_constt * width)) *
                  (kappa[i][j] / mu)) /
                     ((1 - solid_fraction[i][j]) * time));
      vel_darcyx[i][j] =
          0.5 * (vel_darcyx[i][j] -
                 (((Pf[i + 1][j] - Pf[i - 1][j]) / (2 * scale_constt * width)) *
                  (kappa[i][j] / mu)) /
                     ((1 - solid_fraction[i][j]) * time));

      vel_darcyx_in[i][j] =
          (((grad_xold_in * deltaPf) / (scale_constt * width)) *
           (kappa[i][j] / mu)) /
          (1 - solid_fraction[i][j]);
      vel_darcyy_in[i][j] =
          (((grad_yold_in * deltaPf) / (scale_constt * width)) *
           (kappa[i][j] / mu)) /
          (1 - solid_fraction[i][j]);

      vel_darcyy_in[i][j] +=
          -(((Pf[i][j] - Pf[i][j - 1]) / (2 * scale_constt * width)) *
            (kappa[i][j] / mu)) /
          (1 - solid_fraction[i][j]);
      vel_darcyx_in[i][j] +=
          -(((Pf[i][j] - Pf[i - 1][j]) / (2 * scale_constt * width)) *
            (kappa[i][j] / mu)) /
          (1 - solid_fraction[i][j]);

      vel_darcyx_out[i][j] =
          (((grad_xold_out * deltaPf) / (scale_constt * width)) *
           (kappa[i][j] / mu)) /
          (1 - solid_fraction[i][j]);
      vel_darcyy_out[i][j] =
          (((grad_yold_out * deltaPf) / (scale_constt * width)) *
           (kappa[i][j] / mu)) /
          (1 - solid_fraction[i][j]);

      vel_darcyy_out[i][j] +=
          -(((Pf[i][j + 1] - Pf[i][j]) / (2 * scale_constt * width)) *
            (kappa[i][j] / mu)) /
          (1 - solid_fraction[i][j]);
      vel_darcyx_out[i][j] +=
          -(((Pf[i + 1][j] - Pf[i][j]) / (2 * scale_constt * width)) *
            (kappa[i][j] / mu)) /
          (1 - solid_fraction[i][j]);

      // if (vel_darcyy[i][j] > 0.0000001)
      //{
      // cout << Pf[i][j] - oldPf[i][j] << endl;
      // cout << Pf[i][j-1] - Pf[i][j] << endl;
      // cout << vel_darcyy[i][j] << endl;
      //}

      // vel_darcyx[i][j] = grad_xold * deltaPf
      // /((1-solid_fraction[i][j])*scale_constt*width); vel_darcyy[i][j] =
      // grad_yold * deltaPf /((1-solid_fraction[i][j])*scale_constt*width);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CALCULATION MATRIX INVERSION for the concentration, simpler than pressure,
// does not include source term
//
// based on implicit alternating direction method, here x and y are the same
// length, so only alpha needed
//
// alpha-set
//
// matrices set 1
//
// multiplication
//
// invert
//
// matrices set 2
//
// multiplication invers
//
// Transpose
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation.
// at the moment the matrix size is preset to 200 200 200 for latte simulations
// with 400 particles in x. We are not using sparse matrices, this is something
// that should be changed, would speed up things. This is the prefactor for each
// node in the grid with a constant diffusion coefficient for the moment. this
// prefactor makes the implicit method potentially unstable if its too large,
// this should be checked and corrected for in a future version to make the
// method more stable.
//--------------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::SolveDiffusionImplicit(double dif_const, double time,
                                           double space) {
  alpha_set_dif(dif_const, time);

  matrices_set_1();
  multiplication_dif();
  transpose_trans_dif();
  invert();
  multiplication_inv_dif();

  matrices_set_2();
  multiplication_dif();
  transpose_trans_dif();
  invert();
  multiplication_inv_dif();
}

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation. Note that pressure
// is also in this, meaning that the calucation is easily unstable because it
// varies so much if pressure varies.
//
// this difusion constant is constant at the moment, but could be variable
// depending on porosity/permeability of nodes
//--------------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::alpha_set_dif(double dif_constant,
                                  double time_dif) // need only alpha
{
  int i, j;
  double max;

  max = 0;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      alpha[i][j] = 0;
      {
        alpha[i][j] = time_dif *
                      (dif_constant * 10 * (1.0 - solid_fraction[i][j])) /
                      pow(width * scale_constt, 2.0);
        if (alpha[i][j] > max)
          max = alpha[i][j];
      }
    }
  }

  cout << "dif prefactor:" << max << endl;
}

//----------------------------------------------------------------------------
// multiply diagonal matrix with the fluid pressure matrix
//-----------------------------------------------------------------------------

void Fluid_Lattice::multiplication_dif() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {

      ddyc[i][j] = 0; // ddyc is just used as interim here

      for (k = 0; k < dim; k++) {
        ddyc[i][j] += bl[j][i][k] * Con[k][j];
      }
    }
  }
}

//------------------------------------------------------------------------------------
// Transpose the matrix (switch x and y)
//------------------------------------------------------------------------------------

void Fluid_Lattice::transpose_trans_dif() {
  double con_trans[dim][dim];
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      con_trans[i][j] = 0;

      con_trans[i][j] = ddyc[j][i];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      ddyc[i][j] = 0;

      ddyc[i][j] = con_trans[i][j];
    }
  }
}

//--------------------------------------------------------------------------------
// multiply interim matrix with inverted matrix
//--------------------------------------------------------------------------------

void Fluid_Lattice::multiplication_inv_dif() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      Con[i][j] = 0;

      for (k = 0; k < dim; k++)
        Con[i][j] += a_invl[j][i][k] * ddyc[k][j];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CALCULATION MATRIX INVERSION for the concentration of Metals, simpler than
// pressure, does not include source term
//
// based on implicit alternating direction method, here x and y are the same
// length, so only alpha needed
//
// alpha-set
//
// matrices set 1
//
// multiplication
//
// invert
//
// matrices set 2
//
// multiplication invers
//
// Transpose
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation. Note that pressure
// is also in this, meaning that the calucation is easily unstable because it
// varies so much if pressure varies. at the moment the matrix size is preset to
// 200 200 200 for latte simulations with 400 particles in x we are not using
// sparse matrices, this is something that should be changed, would speed up
// things this is the prefactor for each node in the grid where pf, permeability
// and porosity vary this prefactor makes the implicit method potentially
// unstable if its too large, this should be checked and corrected for in a
// future version to make the method more stable
//--------------------------------------------------------------------------------------------------------------------

double Fluid_Lattice::alpha_set_time() // need only alpha
{
  int i, j;
  double max, correction;
  max = 0.0;

  for (i = 1; i < dim - 1; i++) {
    for (j = 1; j < dim - 1; j++) {
      alpha[i][j] = 0;

      {
        // alpha[i][j] = (1.0 + (compressibility *
        // Pf[i][j]))*density_factor[i][j];
        alpha[i][j] = (1.0 + (compressibility * Pf[i][j]));
        alpha[i][j] *= kappa[i][j] * time_a;
        alpha[i][j] /=
            (2.0 * mu * compressibility * (1 - solid_fraction[i][j]) *
             pow(width * scale_constt, 2.0));

        // cout << alpha[i][j] << endl;

        // if (i==20 & j == 20)
        //{
        // cout << kappa[i][j]*1000000000 << endl;
        // cout << solid_fraction[i][j] << endl;
        // cout << alpha[i][j] << endl;
        //}

        if (alpha[i][j] > max) {
          max = alpha[i][j];
          // cout << i << ";" << j << endl;
          // cout << Pf[i][j] << endl;
          // cout << kappa[i][j]*10000000000 << endl;
          // cout << solid_fraction[i][j] << endl;
        }
      }
    }
  }
  cout << max << endl;

  correction = max / 0.5;
  correction = time_a / correction;

  return correction;
}
//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation. Note that pressure
// is also in this, meaning that the calucation is easily unstable because it
// varies so much if pressure varies. at the moment the matrix size is preset to
// 200 200 200 for latte simulations with 400 particles in x we are not using
// sparse matrices, this is something that should be changed, would speed up
// things this is the prefactor for each node in the grid where pf, permeability
// and porosity vary this prefactor makes the implicit method potentially
// unstable if its too large, this should be checked and corrected for in a
// future version to make the method more stable
//--------------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::alpha_set_fluid(double nonlin_time) // need only alpha
{
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      alpha[i][j] = 0;

      {
        // alpha[i][j] = (1.0 + (compressibility *
        // Pf[i][j]))*density_factor[i][j];
        alpha[i][j] = (1.0 + (compressibility * Pf[i][j]));
        alpha[i][j] *= kappa[i][j] * nonlin_time;
        alpha[i][j] /=
            (2.0 * mu * compressibility * (1 - solid_fraction[i][j]) *
             pow(width * scale_constt, 2.0));
      }
    }
  }
}

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation.
// at the moment the matrix size is preset to 200 200 200 for latte simulations
// with 400 particles in x. We are not using sparse matrices, this is something
// that should be changed, would speed up things. This is the prefactor for each
// node in the grid with a constant diffusion coefficient for the moment. this
// prefactor makes the implicit method potentially unstable if its too large,
// this should be checked and corrected for in a future version to make the
// method more stable.
//--------------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::SolveMetalDiffusionImplicit(double dif_const, double time,
                                                double space) {
  alpha_set_met(dif_const, time);

  matrices_set_1();
  multiplication_met();
  transpose_trans_met();
  invert();
  multiplication_inv_met();

  matrices_set_2();
  multiplication_met();
  transpose_trans_met();
  invert();
  multiplication_inv_met();
}

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation.
// this difusion constant is constant at the moment, but could be variable
// depending on porosity/permeability of nodes
//--------------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::alpha_set_met(double dif_constant,
                                  double time_dif) // need only alpha
{
  int i, j;
  double max;

  max = 0;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      alpha[i][j] = 0;
      {
        alpha[i][j] = time_dif * dif_constant / pow(width * scale_constt, 2.0);
        if (alpha[i][j] > max)
          max = alpha[i][j];
      }
    }
  }

  cout << "dif prefactor:" << max << endl;
}

//----------------------------------------------------------------------------
// multiply diagonal matrix with the fluid pressure matrix
//-----------------------------------------------------------------------------

void Fluid_Lattice::multiplication_met() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {

      ddyc[i][j] = 0; // ddyc is just used as interim here

      for (k = 0; k < dim; k++) {
        ddyc[i][j] += bl[j][i][k] * Metal_A[k][j];
      }
    }
  }
}

//------------------------------------------------------------------------------------
// Transpose the matrix (switch x and y)
//------------------------------------------------------------------------------------

void Fluid_Lattice::transpose_trans_met() {
  double con_trans[dim][dim];
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      con_trans[i][j] = 0;

      con_trans[i][j] = ddyc[j][i];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      ddyc[i][j] = 0;

      ddyc[i][j] = con_trans[i][j];
    }
  }
}

//--------------------------------------------------------------------------------
// multiply interim matrix with inverted matrix
//--------------------------------------------------------------------------------

void Fluid_Lattice::multiplication_inv_met() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      Metal_A[i][j] = 0;

      for (k = 0; k < dim; k++)
        Metal_A[i][j] += a_invl[j][i][k] * ddyc[k][j];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CALCULATION MATRIX INVERSION for the Temperature, simpler than pressure, does
// not include source term
//
// based on implicit alternating direction method, here x and y are the same
// length, so only alpha needed
//
// alpha-set
//
// matrices set 1
//
// multiplication
//
// invert
//
// matrices set 2
//
// multiplication invers
//
// Transpose
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation.
// at the moment the matrix size is preset to 200 200 200 for latte simulations
// with 400 particles in x. We are not using sparse matrices, this is something
// that should be changed, would speed up things. This is the prefactor for each
// node in the grid with a constant diffusion coefficient for the moment. this
// prefactor makes the implicit method potentially unstable if its too large,
// this should be checked and corrected for in a future version to make the
// method more stable.
//--------------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::SolveTemperatureDiffusionImplicitTest(double dif_const,
                                                          double time,
                                                          double space) {
  alpha_set_temp(dif_const, time);

  matrices_set_1();
  transpose_trans_temp();
  invert();
  multiplication_temp();
  transpose_trans_temp();
  multiplication_inv_temp();

  matrices_set_2();
  invert();
  multiplication_temp();
  transpose_trans_temp();
  multiplication_inv_temp();
  transpose_trans_temp();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CALCULATION MATRIX INVERSION for the Temperature, simpler than pressure, does
// not include source term
//
// based on implicit alternating direction method, here x and y are the same
// length, so only alpha needed
//
// alpha-set
//
// matrices set 1
//
// multiplication
//
// invert
//
// matrices set 2
//
// multiplication invers
//
// Transpose
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation.
// at the moment the matrix size is preset to 200 200 200 for latte simulations
// with 400 particles in x. We are not using sparse matrices, this is something
// that should be changed, would speed up things. This is the prefactor for each
// node in the grid with a constant diffusion coefficient for the moment. this
// prefactor makes the implicit method potentially unstable if its too large,
// this should be checked and corrected for in a future version to make the
// method more stable.
//--------------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::SolveTemperatureDiffusionImplicit(double dif_const,
                                                      double time,
                                                      double space) {
  alpha_set_temp(dif_const, time);

  matrices_set_1();
  multiplication_temp();
  transpose_trans_temp();
  invert();
  multiplication_inv_temp();

  matrices_set_2();
  multiplication_temp();
  transpose_trans_temp();
  invert();
  multiplication_inv_temp();
}

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation. Note that pressure
// is also in this, meaning that the calucation is easily unstable because it
// varies so much if pressure varies.
//
// this difusion constant is constant at the moment, but could be variable
// depending on porosity/permeability of nodes
//--------------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::alpha_set_temp(double dif_constant,
                                   double time_dif) // need only alpha
{
  int i, j;
  double max;

  max = 0;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      alpha[i][j] = 0;
      {
        alpha[i][j] = time_dif *
                      (dif_constant * 100 * (1.0 - solid_fraction[i][j])) /
                      pow(width * scale_constt, 2.0);
        if (alpha[i][j] > max)
          max = alpha[i][j];
      }
    }
  }

  cout << "dif prefactor:" << max << endl;
}

//----------------------------------------------------------------------------
// multiply diagonal matrix with the fluid pressure matrix
//-----------------------------------------------------------------------------

void Fluid_Lattice::multiplication_temp() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {

      ddyc[i][j] = 0; // ddyc is just used as interim here

      for (k = 0; k < dim; k++) {
        ddyc[i][j] += bl[j][i][k] * Temperature[k][j];
      }
    }
  }
}

//------------------------------------------------------------------------------------
// Transpose the matrix (switch x and y)
//------------------------------------------------------------------------------------

void Fluid_Lattice::transpose_trans_temp() {
  double con_trans[dim][dim];
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      con_trans[i][j] = 0;

      con_trans[i][j] = ddyc[j][i];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      ddyc[i][j] = 0;

      ddyc[i][j] = con_trans[i][j];
    }
  }
}

//--------------------------------------------------------------------------------
// multiply interim matrix with inverted matrix
//--------------------------------------------------------------------------------

void Fluid_Lattice::multiplication_inv_temp() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      Temperature[i][j] = 0;

      for (k = 0; k < dim; k++)
        Temperature[i][j] += a_invl[j][i][k] * ddyc[k][j];
    }
  }
}

// this is an old version that did not work well and was not stable

void Fluid_Lattice::SolveadvDisExplicitTransient(double fluid_v,
                                                 double dif_const) {
  int i, j, v;
  double D, dt, dx, dy, dx2, dy2, conc_x, conc_y;

  D = dif_const;

  dx = width; //
  dy = width; //
  dx2 = dx * dx;
  for (v = 0; v < dim; v++) {
    Con[v][0] = 0.5; // setting the lower boundary
    Con[v][1] = 0.5;
    // Con[v][dim-1]=0;
  }
  dy2 = dy * dy;

  dt = 0.002;

  // boundary condition (background concentration is already given)

  for (i = 0; i < dim - 1; i++) // x loop
  {
    for (j = 1; j < dim - 2; j++) // y loop
    {

      ddyc[i][j] = (Con[i + 1][j] - 2 * Con[i][j] + Con[i - 1][j]) /
                   dy2; // x has to vary in y, so takes neighbours in i (=y)
      ddxc[i][j] = (Con[i][j + 1] - 2 * Con[i][j] + Con[i][j - 1]) /
                   dx2; // y has to vary in x, so takes neighbours in j (=x)

      // forward approximation takes the difference between node and node behind
      // (with respect to velocity) and adds to node

      // if (vel_darcyx[i][j]>=0.0) // x direction positive, to the right
      if (vel_darcyx[i][j] >= 0.0) {
        dxc[i][j] = (Con[i][j] - Con[i - 1][j]) / dx;
        if (dxc[i][j] > 0.0)
          dxc[i][j] = 0.0;
      } else // x direction negative, to the left
      {
        dxc[i][j] = (Con[i][j] - Con[i + 1][j]) / dx;
        if (dxc[i][j] > 0.0)
          dxc[i][j] = 0.0;
      }
      if (vel_darcyy[i][j] >= 0.0) // y direction positive, upwards
      {
        dyc[i][j] = (Con[i][j] - Con[i][j - 1]) / dy;
        if (dyc[i][j] > 0.0)
          dyc[i][j] = 0.0;
      } else // y direction positive, downwards
      {
        dyc[i][j] = (Con[i][j] - Con[i][j + 1]) / dy;
        if (dyc[i][j] > 0.0)
          dyc[i][j] = 0.0;
      }
    }
  }

  // you now have all the values to fill in the matrix and update the
  // concentration. You have to do that after the calculation, otherwise you get
  // a trend in the calculation depending on where you start because some have
  // been updated and others not. Need to update in on step.

  for (i = 0; i < dim - 1; i++) // y
  {
    for (j = 1; j < dim - 2; j++) // x
    {
      // diffusion
      // Con[i][j] = Con[i][j]+dt*D*(ddxc[i][j]+ddyc[i][j]);
      // cout << " too large y" << dt*abs(vel_darcyy[i][j])*fluid_v*dyc[i][j] <<
      // endl;
      // advection
      // if (dt*abs(vel_darcyx[i][j])*fluid_v*dxc[i][j]>0.1)
      // cout << " too large x" << dt*abs(vel_darcyx[i][j])*fluid_v*dxc[i][j] <<
      // endl; if (dt*abs(vel_darcyy[i][j])*fluid_v*dyc[i][j]<-0.0) cout << "
      // too large y" << dt*abs(vel_darcyy[i][j])*fluid_v*dyc[i][j] << endl;
      // else

      conc_x = Con[i][j] - (dt * abs(vel_darcyx[i][j]) * fluid_v * dxc[i][j]);
      conc_y = Con[i][j] - (dt * abs(vel_darcyy[i][j]) * fluid_v * dyc[i][j]);

      if (j < 4) {
        cout << "concentration x" << conc_x << endl;
        cout << "concentration y" << conc_y << endl;
      }

      Con[i][j] = conc_x + conc_y;

      // else if (vel_darcyy[i][j]< -0.00000001)
      //{
      //     Con[i][j] = Con[i][j]-(dt*vel_darcyx[i][j]*fluid_v*dxc[i][j]);
      //     Con[i][j] = Con[i][j]-(dt*vel_darcyy[i][j]*fluid_v*dyc[i][j]);
      // }

      // original

      // concent[i][j] = concent[i][j]-(dt*1000000000000*kappa[i][j]*dyc[i][j]);
      if (Con[i][j] < 0.0)
        Con[i][j] = 0.0;
      else if (Con[i][j] > 0.5)
        Con[i][j] = 0.5;
    }
  }

  for (v = 1; v < dim - 1;
       v++) // boundary condition for sides, not really wrapping... sort of
  {
    Con[0][v] = Con[1][v];

    Con[dim - 1][v] = Con[dim - 2][v];
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CALCULATION MATRIX INVERSION Phreeqc for the concentration, simpler than
// pressure, does not include source term
//
// based on implicit alternating direction method, here x and y are the same
// length, so only alpha needed
//
// alpha-set
//
// matrices set 1
//
// multiplication
//
// invert
//
// matrices set 2
//
// multiplication invers
//
// Transpose
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation.
// at the moment the matrix size is preset to 200 200 200 for latte simulations
// with 400 particles in x. We are not using sparse matrices, this is something
// that should be changed, would speed up things. This is the prefactor for each
// node in the grid with a constant diffusion coefficient for the moment. this
// prefactor makes the implicit method potentially unstable if its too large,
// this should be checked and corrected for in a future version to make the
// method more stable.
//--------------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::SolveDiffusionImplicit_Phreeqc(double dif_const,
                                                   double time, double space) {
  alpha_set_dif(dif_const, time);

  matrices_set_1();
  multiplication_dif();
  transpose_trans_dif();
  invert();
  multiplication_inv_dif();

  matrices_set_2();
  multiplication_dif();
  transpose_trans_dif();
  invert();
  multiplication_inv_dif();
}

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation. Note that pressure
// is also in this, meaning that the calucation is easily unstable because it
// varies so much if pressure varies.
//
// this difusion constant is constant at the moment, but could be variable
// depending on porosity/permeability of nodes
//--------------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::alpha_set_Phreeqc(double dif_constant,
                                      double time_dif) // need only alpha
{
  int i, j;
  double max;

  max = 0;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      alpha[i][j] = 0;
      {
        alpha[i][j] = time_dif * dif_constant / pow(width * scale_constt, 2.0);
        if (alpha[i][j] > max)
          max = alpha[i][j];
      }
    }
  }

  cout << "dif prefactor:" << max << endl;
}

//----------------------------------------------------------------------------
// multiply diagonal matrix with the fluid pressure matrix
//-----------------------------------------------------------------------------

void Fluid_Lattice::multiplication_Phreeqc() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {

      ddyc[i][j] = 0; // ddyc is just used as interim here

      for (k = 0; k < dim; k++) {
        ddyc[i][j] += bl[j][i][k] * Temperature[k][j];
      }
    }
  }
}

//------------------------------------------------------------------------------------
// Transpose the matrix (switch x and y)
//------------------------------------------------------------------------------------

void Fluid_Lattice::transpose_trans_Phreeqc() {
  double con_trans[dim][dim];
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      con_trans[i][j] = 0;

      con_trans[i][j] = ddyc[j][i];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      ddyc[i][j] = 0;

      ddyc[i][j] = con_trans[i][j];
    }
  }
}

//--------------------------------------------------------------------------------
// multiply interim matrix with inverted matrix
//--------------------------------------------------------------------------------

void Fluid_Lattice::multiplication_inv_Phreeqc() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      Temperature[i][j] = 0;

      for (k = 0; k < dim; k++)
        Temperature[i][j] += a_invl[j][i][k] * ddyc[k][j];
    }
  }
}

//---------------------------------------------------------------------------------
// Phreeqc code created by Uli 2020

//---------------------------------------------------------------------------------

void Fluid_Lattice::SetFluidPhreeqc(const char *pore_fluid,
                                    const char *inf_fluid, const char *d_base) {
  iphreeqc_obj.Version();
  iphreeqc_obj.time_step = 0;
  iphreeqc_obj.SetDatabase(d_base);
  iphreeqc_obj.SetPoreFluid(pore_fluid);
  iphreeqc_obj.SetInfFluid(inf_fluid);

  iphreeqc_obj.InitializeFluidVector();
  iphreeqc_obj.AdjustUnits(iphreeqc_obj.Conc_inf);

  // Set the pore fluid
  iphreeqc_obj.Conc_por_b = iphreeqc_obj.Conc_por;
  iphreeqc_obj.Equilibrate(iphreeqc_obj.Conc_por, "calcite");
  iphreeqc_obj.Equilibrate(iphreeqc_obj.Conc_por_b, "quartz");

  iphreeqc_obj.SetFluidMatrices(iphreeqc_obj.Conc_por, 0);
  iphreeqc_obj.AddPoreFluidLayer(iphreeqc_obj.Conc_por_b, 2, 30);
  iphreeqc_obj.AddPoreFluidLayer(iphreeqc_obj.Conc_por_b, 50, dimY);
}

void Fluid_Lattice::Transport(double timestep, double space, int adv) {
  // double conc_Zn = iphreeqc.FindConcentration("Zn");
  // cout << "\nConcentration of Zn in the infiltrating fluid: " << conc_Zn <<
  // endl;
  /**We have two types of fluid data structures, one large container that holds
  the concentrations for the entire domain and one simple vector that contains
  the concentrations of infiltrating fluid. In fluid_lattice we have an object
  adn a pointer to the phreeqc class. However, for teh parallelization we need
  to create new objects and load the databse for each of them **/

  for (PhreeqC::FLUID ::iterator it = iphreeqc_obj.Conc_inf.begin();
       it != iphreeqc_obj.Conc_inf.end(); it++) {
    if (it->first != "pH") {
      const char *type = it->first; // nmae of the element or parameter
      double constant;              // diffusion constant

      iphreeqc_obj.GetConcentration(
          type, conc); // this is the function that fills the matrix

      constant = 0.01; // 10 e -7

      // Transport
      //  first advection

      cout << it->first << conc[10][2] << endl;

      // SolveAdvectionLax_Wenddroff_Phreeqc(timestep * 1000, space, 2);

      // then diffusion

      for (int i = 0; i < dimX; i++) {
        for (int j = 0; j < dimY; j++) {
          if (type == "temp")
            Temperature[i][j] = conc[i][j];
          else
            Con[i][j] = conc[i][j];
        }
      }

      if (type == "temp") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff_Temperature(timestep, space, 2);
        SolveTemperatureDiffusionImplicitTest(constant * 1, timestep,
                                              space); // 10 e-6
      } else if (type == "S(6)") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[0]);
        // SolveAdvectionFromm(space, timestep);
        SolveDiffusionImplicit(constant * 0.01, timestep, space); // 10 e-9
      } else if (type == "Mg") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[1]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "Ca") {
        // cout << "max Ca"  << max_conc[2] << endl;
        // SolveAdvectionFromm(space, timestep);
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[2]);
        SolveDiffusionImplicit(constant * 0.0071, timestep, space);
      } else if (type == "Zn") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[3]);
        // SolveAdvectionFromm(space, timestep);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "Na") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[4]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "K") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[5]);
        SolveDiffusionImplicit(constant * 0.0177, timestep, space);
      } else if (type == "Ba") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[6]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "Si") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[7]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "Cl") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[8]);
        SolveDiffusionImplicit(constant * 0.0181, timestep, space);
      } else if (type == "Fe") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[9]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "Mn") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[10]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "Pb") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[11]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      }

      for (int i = 0; i < dimX; i++) {
        for (int j = 0; j < dimY; j++) {
          if (type == "temp")
            conc[i][j] = Temperature[i][j];
          else {
            conc[i][j] = Con[i][j];
            if (type == "S(6)")
              Con_S6[i][j] = Con[i][j];
            else if (type == "Mg")
              Con_Mg[i][j] = Con[i][j];
            else if (type == "Ca")
              Con_Ca[i][j] = Con[i][j];
            else if (type == "Zn")
              Con_Zn[i][j] = Con[i][j];
          }
        }
      }

      iphreeqc_obj.SetConcentration(
          type,
          conc); // this function passes the matric back to the fluid array
    }
  }
}

void Fluid_Lattice::GetOutput(vector<std::string> outputs) {
  cout << "phreeqC calculating fluid properties" << endl;
  boost::progress_display *show_progress =
      new boost::progress_display(dimX * dimY);
  time_t start, end;
  time(&start);

  typedef vector<std::pair<std::string, boost::array<double, (dimY * dimX)>>>
      RESULTS;
  RESULTS results;
  results.reserve(outputs.size());

  boost::array<double, (dimY * dimX)> matrix;
  matrix.fill(0);

  for (int i = 0; i < outputs.size(); i++)
    results.push_back(make_pair(outputs[i], matrix));

  // Create per-thread PhreeqC instances inside the parallel region to avoid
  // copying
#pragma omp parallel shared(results) num_threads(iphreeqc_obj.nb_threads)
  {
    // this is private for every thread (I hope...)
    PhreeqC iphreeqc_o; // construct a fresh instance (IPhreeqc is non-copyable)
    // initialize from the shared config
    iphreeqc_o.nb_threads = iphreeqc_obj.nb_threads;
    iphreeqc_o.p_database = iphreeqc_obj.p_database;
    iphreeqc_o.loader(
        iphreeqc_obj.p_database); // load the database for this object
    iphreeqc_o.Conc_matrices =
        iphreeqc_obj.Conc_matrices; // copy the fluid compositions
    iphreeqc_o.units =
        iphreeqc_obj.units; // copy units (TODO: check whats happening with the
                            // units)					!!!
    iphreeqc_o.time_step =
        iphreeqc_obj.time_step; // also copy the current time step

#pragma omp for collapse(2)
    for (int i = 0; i < dimX; i++) {
      for (int j = 0; j < dimY; j++) {
        // printf("i = %d, j= %d, threadId = %d \n", i, j,
        // omp_get_thread_num());				//for degugging
        // we loop which thread got which entry
        PhreeqC::FLUID loc_fluid =
            iphreeqc_o.GetFluid(i, j, iphreeqc_o.Conc_matrices);

        vector<std::pair<std::string, double>> local_result;
        for (int o = 0; o < outputs.size(); o++)
          local_result.push_back(make_pair(outputs[o], 0));

        iphreeqc_o.OUTPUT(iphreeqc_o.time_step, loc_fluid, local_result, false);

        for (vector<pair<std::string, double>>::const_iterator it =
                 local_result.begin();
             it != local_result.end(); it++) {
          for (RESULTS ::iterator it2 = results.begin(); it2 != results.end();
               it2++) {
            if (boost::iequals(it2->first, it->first))
              results.at(it2 - results.begin()).second[i * dimY + j] =
                  it->second;
            ;
          }
        }
        ++(*show_progress);
      }
    }
  }
  cout << " done..." << endl;
  cout << "Passing over from fluid to particle..." << endl;

  /// TODO: check how to pass over to particles ///

  for (RESULTS ::iterator it = results.begin(); it != results.end(); it++) {
    if (it->first == "density") {
      for (int x = 0; x < dimX; x++)
        for (int y = 0; y < dimY; y++) {
          rho_f[x][y] = it->second[x * dimY + y];
        }
    } else if (it->first == "barite") {
      for (int x = 0; x < dimX; x++)
        for (int y = 0; y < dimY; y++) {
          Si[x][y] = it->second[x * dimY + y];
        }
    }
  }

  // std::copy(&Si_calcite[0][0], &Si_calcite[0][0]+dimX*dimY,&Si[0][0]);
  cout << endl;
  time(&end);

  // Calculating total time taken by the program.
  double time_taken = double(end - start);
  cout << "Time taken by phreeqc: " << fixed << time_taken << setprecision(10);
  cout << " sec " << endl;
  iphreeqc_obj.time_step += 1;
}

void Fluid_Lattice::TestPhreeqc(double timestep, double space, int adv) {
  // double conc_Zn = iphreeqc.FindConcentration("Zn");
  // cout << "\nConcentration of Zn in the infiltrating fluid: " << conc_Zn <<
  // endl;
  /**We have two types of fluid data structures, one large container that holds
  the concentrations for the entire domain and one simple vector that contains
  the concentrations of infiltrating fluid. In fluid_lattice we have an object
  adn a pointer to the phreeqc class. However, for teh parallelization we need
  to create new objects and load the databse for each of them **/

  for (PhreeqC::FLUID ::iterator it = iphreeqc_obj.Conc_inf.begin();
       it != iphreeqc_obj.Conc_inf.end(); it++) {
    if (it->first != "pH") {
      const char *type = it->first; // nmae of the element or parameter
      // double inf_conc = it->second ; //concentration in filitrating fluid
      double constant; // diffusion constant

      // double conc[dimX][dimY]; 					//
      // this is the matrix we use for transport calculation
      iphreeqc_obj.GetConcentration(
          type, conc); // this is the function that fills the matrix

      constant = 0.01; // 10 e -7

      // Transport
      //  first advection

      cout << it->first << conc[10][2] << endl;

      // SolveAdvectionLax_Wenddroff_Phreeqc(timestep * 1000, space, 2);

      // then diffusion

      for (int i = 0; i < dimX; i++) {
        for (int j = 0; j < dimY; j++) {
          if (type == "temp")
            Temperature[i][j] = conc[i][j];
          else
            Con[i][j] = conc[i][j];
        }
      }

      if (type == "temp") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff_Temperature(timestep, space, 2);
        SolveTemperatureDiffusionImplicitTest(constant * 1, timestep,
                                              space); // 10 e-6
      } else if (type == "S(6)") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[0]);
        // SolveAdvectionFromm(space, timestep);
        SolveDiffusionImplicit(constant * 0.01, timestep, space); // 10 e-9
      } else if (type == "Mg") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[1]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "Ca") {
        // cout << "max Ca"  << max_conc[2] << endl;
        // SolveAdvectionFromm(space, timestep);
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[2]);
        SolveDiffusionImplicit(constant * 0.0071, timestep, space);
      } else if (type == "Zn") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[3]);
        // SolveAdvectionFromm(space, timestep);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "Na") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[4]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "K") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[5]);
        SolveDiffusionImplicit(constant * 0.0177, timestep, space);
      } else if (type == "Ba") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[6]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "Si") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[7]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "Cl") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[8]);
        SolveDiffusionImplicit(constant * 0.0181, timestep, space);
      } else if (type == "Fe") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[9]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "Mn") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[10]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      } else if (type == "Pb") {
        if (adv == 1)
          SolveAdvectionLax_Wenddroff(timestep, space, 2, max_conc[11]);
        SolveDiffusionImplicit(constant * 0.01, timestep, space);
      }

      for (int i = 0; i < dimX; i++) {
        for (int j = 0; j < dimY; j++) {
          if (type == "temp")
            conc[i][j] = Temperature[i][j];
          else {
            conc[i][j] = Con[i][j];
            if (type == "S(6)")
              Con_S6[i][j] = Con[i][j];
            else if (type == "Mg")
              Con_Mg[i][j] = Con[i][j];
            else if (type == "Ca")
              Con_Ca[i][j] = Con[i][j];
            else if (type == "Zn")
              Con_Zn[i][j] = Con[i][j];
          }
        }
      }

      iphreeqc_obj.SetConcentration(
          type,
          conc); // this function passes the matric back to the fluid array
    }
  }

  // here we do the parallelization
  cout << "Reaction calculation " << endl;
  boost::progress_display *show_progress =
      new boost::progress_display(dimX * dimY);
  time_t start, end;
  time(&start);

  double Si_calcite[dimX][dimY];
  double Si_dolomite[dimX][dimY];
  double Si_FeOH3[dimX][dimY];
  double Si_Goethite[dimX][dimY];
  double Si_Jarosite[dimX][dimY];
  double Si_Talc[dimX][dimY];
  double Si_Hamatite[dimX][dimY];

  // Create per-thread PhreeqC instances inside the parallel region to avoid
  // copying

  // at teh moment only calulation saturation state(s). We need to be careful
  // with the parallel section below. a private class instance is created after
  // the first pragma. We have to define shared variables and these CANNOT be a
  // member of a class (i.e. Si).

#pragma omp parallel shared(Si_calcite, Si_dolomite)                           \
    num_threads(iphreeqc_obj.nb_threads)
  {
    // this is private for every thread (I hope...)
    PhreeqC iphreeqc_o; // construct a fresh instance (IPhreeqc is non-copyable)
    // initialize from the shared config
    iphreeqc_o.nb_threads = iphreeqc_obj.nb_threads;
    iphreeqc_o.p_database = iphreeqc_obj.p_database;
    iphreeqc_o.loader(
        iphreeqc_obj.p_database); // load the database for this object
    iphreeqc_o.Conc_matrices =
        iphreeqc_obj.Conc_matrices;        // copy the fluid compositions
    iphreeqc_o.units = iphreeqc_obj.units; // copy units
    iphreeqc_o.time_step =
        iphreeqc_obj.time_step; // also copy the current time step

// this is the parallel loop, we collapse the nested loop because we can :-)
#pragma omp for collapse(2)

    for (int i = 1; i < dimX - 1; i++) {
      for (int j = 1; j < dimY - 1; j++) {
        double si_cal, si_Fe, si_Goe, si_Ha, si_Jar, si_Talc, si_dol;
        // printf("i = %d, j= %d, threadId = %d \n", i, j,
        // omp_get_thread_num());				//for degugging
        // we loop which thread got which entry
        ///* we first get the local fluid composition. I think we will need to
        /// adjust things here to
        /// minimize the error percentage of Anion adn Cation charge
        /// balance.*///
        PhreeqC::FLUID loc_fluid =
            iphreeqc_o.GetFluid(i, j, iphreeqc_o.Conc_matrices);

        // every 10th time step we write a file containing the saturation states
        // of the particel located at X and Y:
        int X = dimX / 2;
        int Y = 30;
        if ((iphreeqc_obj.time_step % 10) == 0) {
          if (i == X && j == Y) {
            si_cal =
                iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid, "barite", true);
            // si_Ha  = iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid,
            // "hematite", true); si_Fe  = iphreeqc_o.SI(iphreeqc_o.time_step,
            // loc_fluid, "fe(oh)3", true); si_Goe  =
            // iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid, "goethite", true);
            // si_Jar  = iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid,
            // "jarosite-k", true); si_Talc  =
            // iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid, "talc", true);
            // si_dol  = iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid,
            // "dolomite", true);
          } else
            si_cal =
                iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid, "barite", false);
          // si_Ha = iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid, "hematite",
          // false); si_Fe  = iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid,
          // "fe(oh)3", false); si_Goe  = iphreeqc_o.SI(iphreeqc_o.time_step,
          // loc_fluid, "goethite", false); si_Jar  =
          // iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid, "jarosite-k",
          // false); si_Talc  = iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid,
          // "talc", false); si_dol  = iphreeqc_o.SI(iphreeqc_o.time_step,
          // loc_fluid, "dolomite", false);
        } else {
          si_cal =
              iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid, "barite", false);
          // si_Ha = iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid, "hematite",
          // false); si_Fe  = iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid,
          // "fe(oh)3", false); si_Goe  = iphreeqc_o.SI(iphreeqc_o.time_step,
          // loc_fluid, "goethite", false); si_Jar  =
          // iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid, "jarosite-k",
          // false); si_Talc  = iphreeqc_o.SI(iphreeqc_o.time_step, loc_fluid,
          // "talc", false); si_dol  = iphreeqc_o.SI(iphreeqc_o.time_step,
          // loc_fluid, "dolomite", false);
        }

        // here is the horrible function that calulates reaction rates (full of
        // magic numbers). To see what they are look in phreeqc.cc
        // iphreeqc_o.React(iphreeqc_o.GetFluid(i, j), -13.7, 0.67, 22.6,
        // 0.162,103, 15, 1.5768e8, "Quartz", "SiO2", false);

        Si_calcite[i][j] = si_cal; // fill the matrix
        // Si_FeOH3[i][j] = si_Fe;
        // Si_Goethite[i][j] = si_Goe;
        // Si_Jarosite[dimX][dimY] = si_Jar;
        // Si_Talc[dimX][dimY] = si_Talc;
        // Si_Hamatite[i][j] = si_Ha;
        // Si_dolomite[i][j] = si_dol;
        ++(*show_progress);
      }
    }
    for (int i = 0; i < dimY; i++) {
      Si_calcite[i][0] = 0;
      Si_calcite[i][dimY - 1] = 0;
      Si_calcite[0][i] = 0;
      Si_calcite[dimY - 1][i] = 0;
      // Si_dolomite[i][0]  = 0;
      // Si_dolomite[i][dimY-1]  = 0;
      // Si_dolomite[0][i]  = 0;
      // Si_dolomite[dimY-1][i]  = 0;
      // Si_FeOH3[i][0]  = 0;
      // Si_FeOH3[i][dimY-1]  = 0;
      // Si_FeOH3[0][i]  = 0;
      // Si_FeOH3[dimY-1][i]  = 0;
      // Si_Goethite[i][0]  = 0;
      // Si_Goethite[i][dimY-1]  = 0;
      // Si_Goethite[0][i]  = 0;
      // Si_Goethite[dimY-1][i]  = 0;
      // Si_Jarosite[i][0]  = 0;
      // Si_Jarosite[i][dimY-1]  = 0;
      // Si_Jarosite[0][i]  = 0;
      // Si_Jarosite[dimY-1][i]  = 0;
      // Si_Talc[i][0]  = 0;
      // Si_Talc[i][dimY-1]  = 0;
      // Si_Talc[0][i]  = 0;
      // Si_Talc[dimY-1][i]  = 0;
      // Si_Hamatite[i][0]  = 0;
      // Si_Hamatite[i][dimY-1]  = 0;
      // Si_Hamatite[0][i]  = 0;
      // Si_Hamatite[dimY-1][i]  = 0;
    }
  }
  std::copy(&Si_calcite[0][0], &Si_calcite[0][0] + dimX * dimY, &Si[0][0]);
  // std::copy(&Si_dolomite[0][0], &Si_dolomite[0][0]+dimX*dimY,&Si2[0][0]);
  // std::copy(&Si_FeOH3[0][0], &Si_FeOH3[0][0]+dimX*dimY,&Si3[0][0]);
  // std::copy(&Si_Goethite[0][0], &Si_Goethite[0][0]+dimX*dimY,&Si4[0][0]);
  // std::copy(&Si_Jarosite[0][0], &Si_Jarosite[0][0]+dimX*dimY,&Si5[0][0]);
  // std::copy(&Si_Talc[0][0], &Si_Talc[0][0]+dimX*dimY,&Si6[0][0]);
  // std::copy(&Si_Hamatite[0][0], &Si_Hamatite[0][0]+dimX*dimY,&Si7[0][0]);
  cout << endl;

  time(&end);
  // Calculating total time taken by the program.
  double time_taken = double(end - start);
  cout << "Time taken by phreeqc: " << fixed << time_taken << setprecision(10);
  cout << " sec " << endl;
  iphreeqc_obj.time_step += 1;
}

void Fluid_Lattice::InitConcPhreeqc() {
  double max;

  /**We have two types of fluid data structures, one large container that holds
  the concentrations for the entire domain and one simple vector that contains
  the concentrations of infiltrating fluid. In fluid_lattice we have an object
  adn a pointer to the phreeqc class. However, for teh parallelization we need
  to create new objects and load the databse for each of them **/

  for (int jj = 0; jj < 12; jj++)
    max_conc[jj] = 0;

  for (PhreeqC::FLUID ::iterator it = iphreeqc_obj.Conc_inf.begin();
       it != iphreeqc_obj.Conc_inf.end(); it++) {
    if (it->first != "pH") {
      max = 0;

      const char *type = it->first; // nmae of the element or parameter
      // double inf_conc = it->second ; //concentration in filitrating fluid

      // double conc[dimX][dimY]; 					//
      // this is the matrix we use for transport calculation
      iphreeqc_obj.GetConcentration(
          type, conc); // this is the function that fills the matrix

      for (int i = 0; i < dimX; i++) {
        for (int j = 0; j < dimY; j++) {
          if (conc[i][j] > max)
            max = conc[i][j];
        }
      }
      if (type == "S(6)")
        max_conc[0] = max;
      else if (type == "Mg")
        max_conc[1] = max;
      else if (type == "Ca")
        max_conc[2] = max;
      else if (type == "Zn")
        max_conc[3] = max;
      else if (type == "Na")
        max_conc[4] = max;
      else if (type == "K")
        max_conc[5] = max;
      else if (type == "Ba")
        max_conc[6] = max;
      else if (type == "Si")
        max_conc[7] = max;
      else if (type == "Cl")
        max_conc[8] = max;
      else if (type == "Fe")
        max_conc[9] = max;
      else if (type == "Mn")
        max_conc[10] = max;
      else if (type == "Pb")
        max_conc[11] = max;
    }
  }
}

void Fluid_Lattice::SolveAdvectionFromm(double space, double time) {
  int i, j, k, iteration, kk;
  double dx, dt, maxprefactor, prefactor;

  dx = space / 100;
  dt = time;

  // determine the largest prefactor

  maxprefactor = 0.0;

  for (i = 1; i < dim - 1; i++) // x loop
  {
    for (j = 1; j < dim - 1; j++) // y loop
    {
      prefactor = fabs(vel_darcyx_in[i][j] * dt / dx);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
      prefactor = fabs(vel_darcyy_in[i][j] * dt / dx);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
      prefactor = fabs(vel_darcyx_out[i][j] * dt / dx);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
      prefactor = fabs(vel_darcyy_out[i][j] * dt / dx);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
    }
  }

  cout << "advection prefactor:" << maxprefactor << endl;

  // adjust prefactor to be 0.2, this should be an internal loop
  // so that the time is constant.

  iteration = maxprefactor / 0.02;

  if (iteration < 1)
    iteration = 1;
  else
    dt = dt / iteration;

  cout << "iteration advection:" << iteration << endl;

  for (kk = 0; kk < iteration; kk++) {

    for (i = 2; i < dim - 2; i++) // x loop
    {
      for (j = 2; j < dim - 2; j++) // y loop
      {
        // cout << "start"  << Con[i][j] << endl;
        dxc[i][j] =
            vel_darcyx_in[i][j] * (0.5 * (Con[i + 1][j] + Con[i][j]) +
                                   0.125 * (Con[i][j] - Con[i - 1][j]) -
                                   0.125 * (Con[i + 2][j] - Con[i + 1][j]));
        // cout << dxc[i][j] << endl;
        dxc[i][j] -= fabs(vel_darcyx_in[i][j]) *
                     (0.25 * (Con[i + 1][j] + Con[i][j]) -
                      0.125 * (Con[i][j] - Con[i - 1][j]) -
                      0.125 * (Con[i + 2][j] - Con[i + 1][j]));
        // cout << dxc[i][j] << endl;
        dxc[i][j] -=
            vel_darcyx_out[i][j] * (0.5 * (Con[i][j] + Con[i - 1][j]) +
                                    0.125 * (Con[i - 1][j] - Con[i - 2][j]) -
                                    0.125 * (Con[i + 1][j] - Con[i][j]));
        // cout << dxc[i][j] << endl;
        dxc[i][j] += fabs(vel_darcyx_out[i][j]) *
                     (0.25 * (Con[i][j] + Con[i - 1][j]) -
                      0.125 * (Con[i - 1][j] - Con[i - 2][j]) -
                      0.125 * (Con[i + 1][j] - Con[i][j]));
        // cout << dxc[i][j] << endl;

        dxc[i][j] +=
            vel_darcyy_in[i][j] * (0.5 * (Con[i][j + 1] + Con[i][j]) +
                                   0.125 * (Con[i][j] - Con[i][j - 1]) -
                                   0.125 * (Con[i][j + 2] - Con[i][j + 1]));
        // cout << dxc[i][j] << endl;
        dxc[i][j] -= fabs(vel_darcyy_in[i][j]) *
                     (0.25 * (Con[i][j + 1] + Con[i][j]) -
                      0.125 * (Con[i][j] - Con[i][j - 1]) -
                      0.125 * (Con[i][j + 2] - Con[i][j + 1]));
        // cout << dxc[i][j] << endl;
        dxc[i][j] -=
            vel_darcyy_out[i][j] * (0.5 * (Con[i][j] + Con[i][j - 1]) +
                                    0.125 * (Con[i][j - 1] - Con[i][j - 2]) -
                                    0.125 * (Con[i][j + 1] - Con[i][j]));
        // cout << dxc[i][j] << endl;
        dxc[i][j] += fabs(vel_darcyy_out[i][j]) *
                     (0.25 * (Con[i][j] + Con[i][j - 1]) -
                      0.125 * (Con[i][j - 1] - Con[i][j - 2]) -
                      0.125 * (Con[i][j + 1] - Con[i][j]));
        // cout << dxc[i][j] << endl;

        dxc[i][j] *= dt / dx;
        // cout << dxc[i][j] << endl;
      }
    }

    // boundaries

    for (k = 1; k < dim - 1; k++) // boundary loop
    {
      //
      i = k;
      j = 1;
      dxc[i][j] =
          vel_darcyx_in[i][j] * (0.5 * (Con[i + 1][j] + Con[i][j]) +
                                 0.125 * (Con[i][j] - Con[i - 1][j]) -
                                 0.125 * (Con[i + 2][j] - Con[i + 1][j]));
      dxc[i][j] -=
          fabs(vel_darcyx_in[i][j]) * (0.25 * (Con[i + 1][j] + Con[i][j]) -
                                       0.125 * (Con[i][j] - Con[i - 1][j]) -
                                       0.125 * (Con[i + 2][j] - Con[i + 1][j]));
      dxc[i][j] -=
          vel_darcyx_out[i][j] * (0.5 * (Con[i][j] + Con[i - 1][j]) +
                                  0.125 * (Con[i - 1][j] - Con[i - 2][j]) -
                                  0.125 * (Con[i + 1][j] - Con[i][j]));
      dxc[i][j] += fabs(vel_darcyx_out[i][j]) *
                   (0.25 * (Con[i][j] + Con[i - 1][j]) -
                    0.125 * (Con[i - 1][j] - Con[i - 2][j]) -
                    0.125 * (Con[i + 1][j] - Con[i][j]));

      dxc[i][j] +=
          vel_darcyy_in[i][j] * (0.5 * (Con[i][j + 1] + Con[i][j]) +
                                 0.125 * (Con[i][j] - Con[i][j - 1]) -
                                 0.125 * (Con[i][j + 2] - Con[i][j + 1]));
      dxc[i][j] -=
          fabs(vel_darcyy_in[i][j]) * (0.25 * (Con[i][j + 1] + Con[i][j]) -
                                       0.125 * (Con[i][j] - Con[i][j - 1]) -
                                       0.125 * (Con[i][j + 2] - Con[i][j + 1]));
      dxc[i][j] -=
          vel_darcyy_out[i][j] * (0.5 * (Con[i][j] + Con[i][j - 1]) +
                                  0.125 * (Con[i][j - 1] - Con[i][j - 2]) -
                                  0.125 * (Con[i][j + 1] - Con[i][j]));
      dxc[i][j] += fabs(vel_darcyy_out[i][j]) *
                   (0.25 * (Con[i][j] + Con[i][j - 1]) -
                    0.125 * (Con[i][j - 1] - Con[i][j - 2]) -
                    0.125 * (Con[i][j + 1] - Con[i][j]));
      dxc[i][j] *= dt / dx;

      i = k;
      j = dim - 2;
      dxc[i][j] =
          vel_darcyx_in[i][j] * (0.5 * (Con[i + 1][j] + Con[i][j]) +
                                 0.125 * (Con[i][j] - Con[i - 1][j]) -
                                 0.125 * (Con[i + 2][j] - Con[i + 1][j]));
      dxc[i][j] -=
          fabs(vel_darcyx_in[i][j]) * (0.25 * (Con[i + 1][j] + Con[i][j]) -
                                       0.125 * (Con[i][j] - Con[i - 1][j]) -
                                       0.125 * (Con[i + 2][j] - Con[i + 1][j]));
      dxc[i][j] -=
          vel_darcyx_out[i][j] * (0.5 * (Con[i][j] + Con[i - 1][j]) +
                                  0.125 * (Con[i - 1][j] - Con[i - 2][j]) -
                                  0.125 * (Con[i + 1][j] - Con[i][j]));
      dxc[i][j] += fabs(vel_darcyx_out[i][j]) *
                   (0.25 * (Con[i][j] + Con[i - 1][j]) -
                    0.125 * (Con[i - 1][j] - Con[i - 2][j]) -
                    0.125 * (Con[i + 1][j] - Con[i][j]));

      dxc[i][j] +=
          vel_darcyy_in[i][j] * (0.5 * (Con[i][j + 1] + Con[i][j]) +
                                 0.125 * (Con[i][j] - Con[i][j - 1]) -
                                 0.125 * (Con[i][j + 2] - Con[i][j + 1]));
      dxc[i][j] -=
          fabs(vel_darcyy_in[i][j]) * (0.25 * (Con[i][j + 1] + Con[i][j]) -
                                       0.125 * (Con[i][j] - Con[i][j - 1]) -
                                       0.125 * (Con[i][j + 2] - Con[i][j + 1]));
      dxc[i][j] -=
          vel_darcyy_out[i][j] * (0.5 * (Con[i][j] + Con[i][j - 1]) +
                                  0.125 * (Con[i][j - 1] - Con[i][j - 2]) -
                                  0.125 * (Con[i][j + 1] - Con[i][j]));
      dxc[i][j] += fabs(vel_darcyy_out[i][j]) *
                   (0.25 * (Con[i][j] + Con[i][j - 1]) -
                    0.125 * (Con[i][j - 1] - Con[i][j - 2]) -
                    0.125 * (Con[i][j + 1] - Con[i][j]));
      dxc[i][j] *= dt / dx;

      i = 1;
      j = k;
      dxc[i][j] =
          vel_darcyx_in[i][j] * (0.5 * (Con[i + 1][j] + Con[i][j]) +
                                 0.125 * (Con[i][j] - Con[i - 1][j]) -
                                 0.125 * (Con[i + 2][j] - Con[i + 1][j]));
      dxc[i][j] -=
          fabs(vel_darcyx_in[i][j]) * (0.25 * (Con[i + 1][j] + Con[i][j]) -
                                       0.125 * (Con[i][j] - Con[i - 1][j]) -
                                       0.125 * (Con[i + 2][j] - Con[i + 1][j]));
      dxc[i][j] -=
          vel_darcyx_out[i][j] * (0.5 * (Con[i][j] + Con[i - 1][j]) +
                                  0.125 * (Con[i - 1][j] - Con[i - 2][j]) -
                                  0.125 * (Con[i + 1][j] - Con[i][j]));
      dxc[i][j] += fabs(vel_darcyx_out[i][j]) *
                   (0.25 * (Con[i][j] + Con[i - 1][j]) -
                    0.125 * (Con[i - 1][j] - Con[i - 2][j]) -
                    0.125 * (Con[i + 1][j] - Con[i][j]));

      dxc[i][j] +=
          vel_darcyy_in[i][j] * (0.5 * (Con[i][j + 1] + Con[i][j]) +
                                 0.125 * (Con[i][j] - Con[i][j - 1]) -
                                 0.125 * (Con[i][j + 2] - Con[i][j + 1]));
      dxc[i][j] -=
          fabs(vel_darcyy_in[i][j]) * (0.25 * (Con[i][j + 1] + Con[i][j]) -
                                       0.125 * (Con[i][j] - Con[i][j - 1]) -
                                       0.125 * (Con[i][j + 2] - Con[i][j + 1]));
      dxc[i][j] -=
          vel_darcyy_out[i][j] * (0.5 * (Con[i][j] + Con[i][j - 1]) +
                                  0.125 * (Con[i][j - 1] - Con[i][j - 2]) -
                                  0.125 * (Con[i][j + 1] - Con[i][j]));
      dxc[i][j] += fabs(vel_darcyy_out[i][j]) *
                   (0.25 * (Con[i][j] + Con[i][j - 1]) -
                    0.125 * (Con[i][j - 1] - Con[i][j - 2]) -
                    0.125 * (Con[i][j + 1] - Con[i][j]));
      dxc[i][j] *= dt / dx;

      i = dim - 2;
      j = k;
      dxc[i][j] =
          vel_darcyx_in[i][j] * (0.5 * (Con[i + 1][j] + Con[i][j]) +
                                 0.125 * (Con[i][j] - Con[i - 1][j]) -
                                 0.125 * (Con[i + 2][j] - Con[i + 1][j]));
      dxc[i][j] -=
          fabs(vel_darcyx_in[i][j]) * (0.25 * (Con[i + 1][j] + Con[i][j]) -
                                       0.125 * (Con[i][j] - Con[i - 1][j]) -
                                       0.125 * (Con[i + 2][j] - Con[i + 1][j]));
      dxc[i][j] -=
          vel_darcyx_out[i][j] * (0.5 * (Con[i][j] + Con[i - 1][j]) +
                                  0.125 * (Con[i - 1][j] - Con[i - 2][j]) -
                                  0.125 * (Con[i + 1][j] - Con[i][j]));
      dxc[i][j] += fabs(vel_darcyx_out[i][j]) *
                   (0.25 * (Con[i][j] + Con[i - 1][j]) -
                    0.125 * (Con[i - 1][j] - Con[i - 2][j]) -
                    0.125 * (Con[i + 1][j] - Con[i][j]));

      dxc[i][j] +=
          vel_darcyy_in[i][j] * (0.5 * (Con[i][j + 1] + Con[i][j]) +
                                 0.125 * (Con[i][j] - Con[i][j - 1]) -
                                 0.125 * (Con[i][j + 2] - Con[i][j + 1]));
      dxc[i][j] -=
          fabs(vel_darcyy_in[i][j]) * (0.25 * (Con[i][j + 1] + Con[i][j]) -
                                       0.125 * (Con[i][j] - Con[i][j - 1]) -
                                       0.125 * (Con[i][j + 2] - Con[i][j + 1]));
      dxc[i][j] -=
          vel_darcyy_out[i][j] * (0.5 * (Con[i][j] + Con[i][j - 1]) +
                                  0.125 * (Con[i][j - 1] - Con[i][j - 2]) -
                                  0.125 * (Con[i][j + 1] - Con[i][j]));
      dxc[i][j] += fabs(vel_darcyy_out[i][j]) *
                   (0.25 * (Con[i][j] + Con[i][j - 1]) -
                    0.125 * (Con[i][j - 1] - Con[i][j - 2]) -
                    0.125 * (Con[i][j + 1] - Con[i][j]));
      dxc[i][j] *= dt / dx;

      // dxc[0][k] = dxc[1][k];
      // dxc[dim-1][k] = dxc[dim-2][k];
    }

    for (i = 1; i < dim - 1; i++) // x loop
    {
      for (j = 1; j < dim - 1; j++) // y loop
      {
        // add up x and y
        cout << "start" << Con[i][j] << endl;
        Con[i][j] = Con[i][j] - dxc[i][j];
        cout << Con[i][j] << endl;

        // final adjustment, no overflow and no negative concentration to
        // improve stability

        // if (Con[i][j] > 0.5)
        //	Con[i][j] = 0.5;
        // if (Con[i][j] < 0.0)
        //	Con[i][j] = 0.0;
      }
    }
  }
}

void Fluid_Lattice::InitHeatSheet(int Temp) {
  int i, j;

  for (i = 0; i < 200; i++) // x loop
  {
    for (j = 0; j < 100; j++) // y loop
    {
      HeatSheet[i][j] = Temp;
    }
  }
}

void Fluid_Lattice::SetHeatSheet(double pos, double pos2) {
  int i, j;

  for (i = 0; i < dim * 2; i++) // x loop
  {
    for (j = 0; j < dim; j++) // y loop
    {
      HeatSheet[i][j] = 10;

      if (i > pos) {
        if (i < pos2) {
          HeatSheet[i][j] = 1000;
        }
      }
    }
  }
}

void Fluid_Lattice::SetHeatSheetbox(double x, double x2, double y, double y2,
                                    int Temp) {
  int i, j;

  for (i = 0; i < 200; i++) // x loop
  {
    for (j = 0; j < 100; j++) // y loop
    {
      if (i > x) {
        if (i < x2) {
          if (j > y) {
            if (j < y2) {
              HeatSheet[i][j] = Temp;
            }
          }
        }
      }
    }
  }
}

void Fluid_Lattice::SolveHeatSheet(double space, double time, double dif) {
  int i, j;

  for (i = 1; i < (200) - 1; i++) {
    for (j = 1; j < dim - 1; j++) {

      // IntHeatSheet[i][j] = HeatSheet[i][j] + (time*dif/
      // pow(space,2.0))*(HeatSheet[i-1][j]-2*HeatSheet[i][j]+
      // HeatSheet[i+1][j]); IntHeatSheet[i][j] = IntHeatSheet[i][j]
      // +HeatSheet[i][j]+  (time*dif/
      // pow(space,2.0))*(HeatSheet[i][j-1]-2*HeatSheet[i][j]+
      // HeatSheet[i][j+1]);
      IntHeatSheet[i][j] =
          HeatSheet[i][j] + (0.1) * (HeatSheet[i - 1][j] - 2 * HeatSheet[i][j] +
                                     HeatSheet[i + 1][j]);
      IntHeatSheet[i][j] = IntHeatSheet[i][j] + HeatSheet[i][j] +
                           (0.1) * (HeatSheet[i][j - 1] - 2 * HeatSheet[i][j] +
                                    HeatSheet[i][j + 1]);
    }
  }

  for (i = 1; i < (200) - 1; i++) {
    for (j = 1; j < dim - 1; j++) {
      HeatSheet[i][j] = 0.5 * IntHeatSheet[i][j];
      if (HeatSheet[i][j] > 200)
        HeatSheet[i][j] = HeatSheet[i][j] - 1;
    }
  }
  for (i = 0; i < (200); i++) {
    HeatSheet[i][0] = HeatSheet[i][1];
    HeatSheet[i][dim - 1] = HeatSheet[i][dim - 2];
  }
  for (j = 1; j < dim - 1; j++) {
    HeatSheet[0][j] = HeatSheet[1][j];
    HeatSheet[199][j] = HeatSheet[198][j];
  }
}

void Fluid_Lattice::Pass_Back_SheetT(
    Particle **
        list) // again gets the replusion box as direct pointer from the lattice
{
  int i, j, k;
  Particle *help; // used as pointer to reach particles

  for (i = 0; i < dim * 2; i++) {
    for (j = 0; j < dim; j++) {
      //--------------------------------------------------------------------
      // 0 is lower left box, 1 lower right box, 2 upper left box
      // and 3 upper right box.
      //
      // we need special boundary conditions for the boundary boxes
      // left and right can be wrappinp, but up and down may be more
      // problematic
      //
      // at the moment right and left wrap, up and down copy next row
      //--------------------------------------------------------------------

      for (k = 0; k < 4; k++) // go to four quadrants of pressure node
      {
        if (k == 0)
          box =
              (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 + ((i - 1) * 2);
        else if (k == 1)
          box =
              (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 + ((i - 1) * 2);
        else if (k == 2)
          box =
              (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 + ((i - 1) * 2);
        else if (k == 3)
          box =
              (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 + ((i - 1) * 2);

        if (list[box]) {
          help = list[box]; // set the pointer

          help->Tdif = HeatSheet[i][j] - help->SheetT;

          help->SheetT = HeatSheet[i][j]; // Mg in solid

          // cout << help->SheetT << endl;

          while (help->next_inBox) // and do the whole thing for other particles
                                   // in the box
          {
            help = help->next_inBox; // and set pointer to next particle

            help->Tdif = HeatSheet[i][j] - help->SheetT;

            help->SheetT = HeatSheet[i][j]; // Mg in solid
          }
        }
      }
    }
  }
}

void Fluid_Lattice::Read_TSheet(Particle **list, int boundary) {
  int i, j, k, particle_in_box, count;
  double solid_average, T;

  // loop through the matrix, i is x, and j is y (note that Irfan did the
  // opposite in pressure lattice) Each fluid box contains 4 repulsion boxes.

  // calculate an average for the boundary conditions. These boundary conditions
  // turn out to be a major pain there are three different ones at the moment,
  // the average is number 3

  count = 0;

  for (i = 1; i < (dim * 2) - 1;
       i++) // loop in x = i and y = j except for boundaries
  {
    for (j = 1; j < dim - 1; j++) {

      particle_in_box = 0; // count sum of particles in the box (FLUID box)
      T = 0;               // solid fraction sum up

      //------------------------------------------------------------------------------------------------
      // In order to get the right solid fraction and to avoid grid effects
      // converting a triangular solid lattice into the square fluid lattice we
      // use a smoothing function and consider solid in the four repulsion boxes
      // of the fluid node plus the surrounding 12 repulsion boxes. The solid is
      // then weight using the smoothing function depending on how far the
      // particle is from the centre of the pressure node.
      //
      // this is also applied to the velocities.
      //-------------------------------------------------------------------------------------------------

      //-------------------------------------------------------------------------------------------------
      // Because of the above the outer rows of the matrix have to be treated
      // separately, otherwise we look into boxes that dont exist.
      //-------------------------------------------------------------------------------------------------

      // if((i > 0 && i < dim-1) && (j > 0 && j < dim-1))
      {
        //--------------------------------------------------------------------------------------------
        // Now define counters for all 16 repulsion boxes that we have to look
        // into to define the solid fraction for the pressure node as a function
        // of i and j. start in the lower left hand corner with k = 0, then go
        // to the right. So 0,1,2,3 are below the actual pressure node, 5,6,9,10
        // are in the actual pressure node.
        //
        // the repulsion box is a one-dimensional list starting in the lower
        // left hand corner of the lattice with 0 and running towards the right
        // with TWICE the lattice width (so for a 100 lattice up to 199). The
        // box in the second horizontal row above 0 is then 200 and so on.
        //
        // at the moment this routine is not using the boundaries because the
        // boxes there are either not completely full or are empty. There is
        // already the try to double up full boxes for the boundaries, but that
        // does not work yet, needs further debugging
        //--------------------------------------------------------------------------------------------

        for (k = 0; k < 16; k++) {
          if (k == 0) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 1 +
                  ((i - 1) * 2); // pl_size = resolution of particle lattice
            if (i == 0)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 3 + ((i - 1) * 2); // take 2
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 1 +
                    ((i - 1) * 2); // take 8
            if (j == 0 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10

          } else if (k == 1) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 + ((i - 1) * 2);
            if (i == 0)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 4 + ((i - 1) * 2); // take 3
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9
            if (j == 0 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10

          } else if (k == 2) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 3 + ((i - 1) * 2);
            if (i == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 1 + ((i - 1) * 2); // take 0
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10
            if (j == 0 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9

          } else if (k == 3) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 + ((i - 1) * 2);
            if (i == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 2 + ((i - 1) * 2); // take 1
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 4 +
                    ((i - 1) * 2); // take 11
            if (j == 0 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9

          }

          else if (k == 4) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 1 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 3 +
                    ((i - 1) * 2); // take 14
            if (j == 0 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10

          } else if (k == 5) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 4 +
                    ((i - 1) * 2); // take 7
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 2 +
                    ((i - 1) * 2); // take 13
            if (j == 0 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10
          } else if (k == 6) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 1 +
                    ((i - 1) * 2); // take 4
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 3 +
                    ((i - 1) * 2); // take 14
            if (j == 0 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9

          } else if (k == 7) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 4 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 4 +
                    ((i - 1) * 2); // take 15
            if (j == 0 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9

          } else if (k == 8) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 1 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10
            if (j == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 1 + ((i - 1) * 2); // take 0
            if (j == dim - 1 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6

          } else if (k == 9) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 4 +
                    ((i - 1) * 2); // take 11
            if (j == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 2 + ((i - 1) * 2); // take 1
            if (j == dim - 1 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6

          } else if (k == 10) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 1 +
                    ((i - 1) * 2); // take 8
            if (j == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 3 + ((i - 1) * 2); // take 2
            if (j == dim - 1 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5

          } else if (k == 11) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 4 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9
            if (j == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 4 + ((i - 1) * 2); // take 3
            if (j == dim - 1 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5

          } else if (k == 12) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 1 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 3 +
                    ((i - 1) * 2); // take 14
            if (j == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 1 +
                    ((i - 1) * 2); // take 4
            if (j == dim - 1 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6

          } else if (k == 13) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 2 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 4 +
                    ((i - 1) * 2); // take 15
            if (j == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5
            if (j == dim - 1 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6

          } else if (k == 14) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 3 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 1 +
                    ((i - 1) * 2); // take 12
            if (j == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6
            if (j == dim - 1 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5

          } else if (k == 15) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 4 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 2 +
                    ((i - 1) * 2); // take 13
            if (j == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 4 +
                    ((i - 1) * 2); // take 7
            if (j == dim - 1 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5
          }

          //-----------------------------------------------------------
          // If the box is filled deal with the particle
          //-----------------------------------------------------------

          if (list[box]) {
            help = list[box]; // help now points at particle

            T = T + help->SheetT;

            // count particles

            particle_in_box++;

            //-----------------------------------------------------------------------
            // check if there is an additional particle in the repulsion box
            // the additional particle would be attached to the help particle
            // with the next_inBox pointer.
            //-----------------------------------------------------------------------

            while (help->next_inBox) {
              // and do the same stuff as above.

              T = T + help->SheetT;

              particle_in_box++;

              help = help->next_inBox;
            }
          }
        }

        // calculate solid fraction
        if (particle_in_box > 0)
          HeatSheet[i][j] = T / particle_in_box;
      }
    }
  }

  if (boundary == 1) {
    // boundary conditions, apply after matrix filled,
    // condition 1: simply copy the next row from the matrix
    // to the boundary row on sides as well as at bottom and
    // top. Works with random background noise. Becomes
    // problematic when Elle grains are used (boundaries then
    // seem to slow things down because they are mainly oriented 90 degrees
    // to the boundary.
    // condition 2: copy neighbour row to boundary row for upper and lower
    // boundary and wrap rows for left and right boundary. Again can be
    // problematic when grain boundaries are used, but does help avoid
    // "doubling" effects during dynamic permeability development condition 3:
    // take an average of the whole lattice and apply that to boundaries. Needs
    // testing.

    for (j = 1; j < dim - 1; j++) // copy neighbour row
    {
      HeatSheet[0][j] = HeatSheet[1][j];

      HeatSheet[dim - 1][j] = HeatSheet[dim - 2][j];
    }

    for (i = 0; i < dim * 2; i++) // copy neighbour row
    {
      HeatSheet[i][0] = HeatSheet[i][1];

      HeatSheet[i][(dim * 2) - 1] = HeatSheet[i][(dim * 2) - 2];
    }
  }
}
