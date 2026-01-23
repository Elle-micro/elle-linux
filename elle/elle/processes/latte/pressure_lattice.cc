#include <iostream>
#include <stdio.h>

#include "pressure_lattice.h"

using namespace std;

// ---------------------------------------------------------------
// Constructor of Pressure_Lattice class
// ---------------------------------------------------------------

Pressure_Lattice::Pressure_Lattice(Particle *first)
    : // -----------------------------------------------------------------------------
      // default variable values. Variables are declared in Pressure_lattice
      // header
      // -----------------------------------------------------------------------------

      x_wall(1.0), y_wall(1.0),

      Seal_hor(true),  // horirontal seal present or not? true for counting
                       // broken bonds in zones
      LowSeal_y(0.45), // lower y_limit of hor_seal
      UppSeal_y(0.55), // upper y_limit of hor_seal

      ModelHalve(true), ModHalf_y(0.5) {
  // f = first;
  // first->radius					//radius of particle.

  time_a = 86400; // one day is 86400!!

  /**********************************************
   * stability/accuracy: alpha[i][j], beta[i][j] < 0.5
   * 5e-10 fur d und dd_dasmal, 1e-09 fur d_dasmal (both with phi*r^2)
   * 5e-10 fur c und cc_dasmal, 1e-09 for c_dasmal (both with phi*r^2)
   *
   * Time_interval		Diffusion		Remained	Remarks
   * 1e-09				99.7%			0.34%
   * unstable 8e-09				53.3%			37.7%
   * unstable 1e-10
   * 44%				56.56%		stable 2e-10
   * 66%				33.41%		stable 3e-10
   * 80%				20%			stable
   * *********************************************/

  // scale_constt = 2000.0; //size of box in meters
  CozenyGrainSize = 0.000001; // in meters - used to calculate permeability
                              // (kappa) from porosity

  rad = 0.02;
  area_circle = PI * pow(rad, 2);

  mu = 1.0e-03;
  // mu = 0.890e-03;   //0.984×10-3 , 1.004*10-3concern PhD disertation "Coupled
  // Flow Discrete Element"

  // The dynamic viscosity of water is 8.90 × 10−4 Pa·s or 8.90 × 10−3 dyn·s/cm2
  // or 0.890 cP at about 25 °C. Water has a viscosity of 0.0091 poise at 25 °C,
  // or 1 centipoise at 20 °C. As a function of temperature T (K): (Pa·s) = A ×
  // 10B/(T−C) where A=2.414 × 10−5 Pa·s ; B = 247.8 K ; and C = 140 K.

  // dynamic shear viscosity of water = 1.0e-3 Pa.s at room tempertaure,
  // whereas volume-viscosity(bulk viscosity) is additional to the usual
  // dynamic-viscosity. Volume-viscosity becomes important only for such effects
  // where fluid compressibility is essential and important. This parameter
  // appears in Navier-Stockes equation, if it is written for compressible
  // fluid. The only known volume viscosity value of simple Newtonian fluids
  // come from old Litovitz and Davis review. They reported volume viscosity of
  // water at 15 C° = 3.09 centipois. 1 P = 1 g/cm.s or 1 Pa.s = 1 Kg/m.s = 10
  // P. water has 0.0091 pois or 91e-5 Pa.s at 25 C°

  compressibility = 4.5e-10; // water compressibility = 4.5e10 m.m/N at 25 C°
  rho_fluid = 1.0e+03;       // water mass density (1000 kg/m³)
  rho_rock = 2.5e+03;        // kg/m³
  g_accl = 9.8;              // gravitational acceleration 9.8 m/s²
}

void Pressure_Lattice::Make_Melt(double change_mu,
                                 double change_compressibility) {
  mu = change_mu;
  compressibility = change_compressibility;
}

void Pressure_Lattice::Change_Time(float change) { time_a = time_a * change; }

void Pressure_Lattice::d_pressure_node_pos_n_area(int resolution) {
  scale_constt = resolution; //(x dimension in meters)

  delta_x = x_wall / (dim - 3);
  delta_y = y_wall / (dim - 3);

  area_node = delta_x * delta_y;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      Pressure_x[i][j] = 0;
      Pressure_y[i][j] = 0;

      Pressure_x[i][j] = (j - 1) * delta_x;
      Pressure_y[i][j] = (i - 1) * delta_y; // j = x direction, i = y direction
    }
  }
}

void Pressure_Lattice::press_nodes(int press_y, int press_x) {
  press_pos_x = Pressure_x[press_y][press_x];
  press_pos_y = Pressure_y[press_y][press_x];

  // cout << press_pos_x << " ::: "  << press_pos_x << endl;
}

void Pressure_Lattice::pressure_background(double init_pressure_val,
                                           int z_distance) {

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      pressure_dist[i][j] = init_pressure_val;
    }
  }
}

void Pressure_Lattice::background_press() {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      press_background[i][j] = 0;

      press_background[i][j] = pressure_dist[i][j];
    }
  }
}

void Pressure_Lattice::source_pressure() {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      sourced_pressure[i][j] = 0;

      sourced_pressure[i][j] = pressure_dist[i][j];
    }
  }
}

void Pressure_Lattice::ohne_background_press() {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++)
      pressure_dist[i][j] -= press_background[i][j];
  }
}

void Pressure_Lattice::mit_background_press() {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++)
      pressure_dist[i][j] += press_background[i][j];
  }
}

/*********************************************************
 * This function is called from Particle_fluid_parameters()
 * in lattice.cc which is called from experiment.cc in
 * order to set the fluid background
 *
 * Should calculate the porosities for rho_s[i][j] and
 *  PArea[i][j]
 *
 * *******************************************************/

void Pressure_Lattice::d_solid_frac_n_vel(Particle **list) // used now
{
  double rho;
  double smoothx, smoothy;
  int res;

  res = (dim - 3) * 2;

  Particles = 0;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      par_rad = 0;
      par_area = 0;

      particle_in_box = 0;
      rho = 0;

      vel_x[i][j] = 0;
      vel_y[i][j] = 0;

      PRadius[i][j] = 0;
      PArea[i][j] = 0;

      rho_s[i][j] = 0;

      if ((i > 1 && i < dim - 2) && (j > 1 && j < dim - 2)) {
        for (k = 0; k < 16; k++) // note that fluid lattice is dim + 3 - empty
                                 // rows around lattice
        {
          if (k == 0)
            box = ((i - 2) * (4 * res)) + 0 +
                  ((j - 2) * 2); // res = resolution of particle lattice
          else if (k == 1)
            box = ((i - 2) * (4 * res)) + 1 + ((j - 2) * 2);
          else if (k == 2)
            box = ((i - 2) * (4 * res)) + 2 + ((j - 2) * 2);
          else if (k == 3)
            box = ((i - 2) * (4 * res)) + 3 + ((j - 2) * 2);

          else if (k == 4)
            box = ((i - 2) * (4 * res)) + (2 * res) + ((j - 2) * 2);
          else if (k == 5)
            box = ((i - 2) * (4 * res)) + (2 * res) + 1 + ((j - 2) * 2);
          else if (k == 6)
            box = ((i - 2) * (4 * res)) + (2 * res) + 2 + ((j - 2) * 2);
          else if (k == 7)
            box = ((i - 2) * (4 * res)) + (2 * res) + 3 + ((j - 2) * 2);

          else if (k == 8)
            box = ((i - 2) * (4 * res)) + (4 * res) + ((j - 2) * 2);
          else if (k == 9)
            box = ((i - 2) * (4 * res)) + (4 * res) + 1 + ((j - 2) * 2);
          else if (k == 10)
            box = ((i - 2) * (4 * res)) + (4 * res) + 2 + ((j - 2) * 2);
          else if (k == 11)
            box = ((i - 2) * (4 * res)) + (4 * res) + 3 + ((j - 2) * 2);

          else if (k == 12)
            box = ((i - 2) * (4 * res)) + (6 * res) + ((j - 2) * 2);
          else if (k == 13)
            box = ((i - 2) * (4 * res)) + (6 * res) + 1 + ((j - 2) * 2);
          else if (k == 14)
            box = ((i - 2) * (4 * res)) + (6 * res) + 2 + ((j - 2) * 2);
          else
            box = ((i - 2) * (4 * res)) + (6 * res) + 3 + ((j - 2) * 2);

          if (list[box]) {
            help = list[box];

            x_diff = help->xpos - Pressure_x[i][j];
            y_diff = help->ypos - Pressure_y[i][j];

            // cout << y_diff << endl;

            if (x_diff < 0)
              x_diff = x_diff * -1;
            if (y_diff < 0)
              y_diff = y_diff * -1;

            if (y_diff < (help->radius / 2.0))
              y_diff = help->radius / 2.0;
            // x_diff = sqrt((x_diff*x_diff)+(y_diff*y_diff));

            if (x_diff > delta_x)
              x_diff = delta_x;
            if (y_diff > delta_y)
              y_diff = delta_y;

            smooth_func =
                (1.0 - (x_diff / delta_x)) * (1.0 - (y_diff / delta_y));

            // smooth_func = 1-(x_diff/delta_x);

            if (help->draw_break) // if there is a fracture solid fraction is
                                  // reduced
            {
              rho = rho + smooth_func * (help->area_par_fluid * 0.8);
            } else {
              rho = rho + smooth_func * help->area_par_fluid;
            }
            // rho = rho + help->area_par_fluid;

            vel_x[i][j] =
                vel_x[i][j] +
                (help->velx * smooth_func); // velocity (flux) per node
            vel_y[i][j] = vel_y[i][j] + (help->vely * smooth_func);

            par_rad += help->rad_par_fluid;
            par_area += help->area_par_fluid;

            particle_in_box++;

            while (help->next_inBox) {
              x_diff = help->next_inBox->xpos - Pressure_x[i][j];
              y_diff = help->next_inBox->ypos - Pressure_y[i][j];

              if (x_diff < 0)
                x_diff = x_diff * -1;
              if (y_diff < 0)
                y_diff = y_diff * -1;

              if (y_diff < (help->radius / 2.0))
                y_diff = help->radius / 2.0;

              if (x_diff > delta_x)
                x_diff = delta_x;
              if (y_diff > delta_y)
                y_diff = delta_y;

              smooth_func =
                  (1.0 - (x_diff / delta_x)) * (1.0 - (y_diff / delta_y));

              //~ rho = rho + smooth_func;

              if (help->draw_break) {
                rho = rho +
                      smooth_func * (help->next_inBox->area_par_fluid * 0.9);
              } else {
                rho = rho + smooth_func * help->next_inBox->area_par_fluid;
              }

              // rho = rho + smooth_func*help->next_inBox->area_par_fluid;

              // rho = rho + help->next_inBox->area_par_fluid;

              vel_x[i][j] =
                  vel_x[i][j] + (help->next_inBox->velx * smooth_func);
              vel_y[i][j] =
                  vel_y[i][j] + (help->next_inBox->vely * smooth_func);

              par_rad += help->next_inBox->rad_par_fluid;
              par_area += help->next_inBox->area_par_fluid;

              particle_in_box++;

              help = help->next_inBox;
            }
          }
        }

        particle[i][j] = particle_in_box;

        par_rad /=
            particle_in_box; // average particle radius w.r.t no. of particles
        par_area /=
            particle_in_box; // average particle area w.r.t no. of particles

        PRadius[i][j] = par_rad;
        // cout << par_rad << endl;
        PArea[i][j] = par_area;
        // cout << rho << endl;
        rho_s[i][j] = rho / area_node;
        // rho_s[i][j] = rho;

        // cout << rho_s[i][j] << endl;

        // rho_s[i][j] *= (2.0/3.0); daniel jul 2013
        rho_s[i][j] *= 0.8;

        if (rho_s[i][j] >= 1.0)
          rho_s[i][j] = 0.999;

        if (rho_s[i][j] <= 0.01)
          rho_s[i][j] = 0.01;

        // cout << "solid fraction" << rho_s[i][j]<< endl;
        //~ rho_s[i][j] = (rho*PArea[i][j])/area_node;
        // rho_s[i][j] = (area_particle * particle[i][j]) / (delta_x * delta_y);

        Particles += particle_in_box;

        vel_x[i][j] = vel_x[i][j] / rho; // average velocity(flux) per node
        vel_y[i][j] = vel_y[i][j] / rho;
      }
      // if(j == dim-1)
      // cout << endl << i << endl;
    }
  }

  for (i = 0; i < dim; i++) // check this boundary condition
  {
    rho_s[i][0] = rho_s[i][2];
    vel_x[i][0] = vel_x[i][2];
    vel_y[i][0] = vel_y[i][2];
    PRadius[i][0] = PRadius[i][2];
    PArea[i][0] = PArea[i][2];

    rho_s[i][1] = rho_s[i][3];
    vel_x[i][1] = vel_x[i][3];
    vel_y[i][1] = vel_y[i][3];
    PRadius[i][1] = PRadius[i][3];
    PArea[i][1] = PArea[i][3];

    rho_s[i][dim - 1] = rho_s[i][dim - 3];
    vel_x[i][dim - 1] = vel_x[i][dim - 3];
    vel_y[i][dim - 1] = vel_y[i][dim - 3];
    PRadius[i][dim - 1] = PRadius[i][dim - 3];
    PArea[i][dim - 1] = PArea[i][dim - 3];

    rho_s[i][dim - 2] = rho_s[i][dim - 4];
    vel_x[i][dim - 2] = vel_x[i][dim - 4];
    vel_y[i][dim - 2] = vel_y[i][dim - 4];
    PRadius[i][dim - 2] = PRadius[i][dim - 4];
    PArea[i][dim - 2] = PArea[i][dim - 4];
  }

  for (j = 0; j < dim; j++) // check this boundary condition
  {
    rho_s[0][j] = rho_s[2][j];
    vel_x[0][j] = vel_x[2][j];
    vel_y[0][j] = vel_y[2][j];
    PRadius[0][j] = PRadius[2][j];
    PArea[0][j] = PArea[2][j];

    rho_s[1][j] = rho_s[3][j];
    vel_x[1][j] = vel_x[3][j];
    vel_y[1][j] = vel_y[3][j];
    PRadius[1][j] = PRadius[3][j];
    PArea[1][j] = PArea[3][j];

    rho_s[dim - 1][j] = rho_s[dim - 3][j];
    vel_x[dim - 1][j] = vel_x[dim - 3][j];
    vel_y[dim - 1][j] = vel_y[dim - 3][j];
    PRadius[dim - 1][j] = PRadius[dim - 3][j];
    PArea[dim - 1][j] = PArea[dim - 3][j];

    rho_s[dim - 2][j] = rho_s[dim - 4][j];
    vel_x[dim - 2][j] = vel_x[dim - 4][j];
    vel_y[dim - 2][j] = vel_y[dim - 4][j];
    PRadius[dim - 2][j] = PRadius[dim - 4][j];
    PArea[dim - 2][j] = PArea[dim - 4][j];
  }
}

void Pressure_Lattice::numberDensity_porosity_kappa() {

  /**************************************************
   * solid fraction = (no. of particles * volume of particle)/volume of cell
   * f/(n*V_p) = 1/ V_cell
   * number density -> n/v_cell = f/V_p
   * ***********************************************/
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      rho[i][j] = 0;

      rho[i][j] = rho_s[i][j] / PArea[i][j];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      porosity[i][j] = 0;

      // rho_s[i][j] *= (2.0/3.0);			//3D fraction of 2D
      // rho_s value.

      porosity[i][j] = 1.0 - rho_s[i][j];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      kappa[i][j] = 0;

      kappa[i][j] = (pow(CozenyGrainSize, 2.0) * pow(porosity[i][j], 3.0));
      kappa[i][j] /= 45.0 * pow((1 - porosity[i][j]), 2.0);

      // cout << kappa[i][j]<< endl;

      if (kappa[i][j] > 0.000000000001) {
        kappa[i][j] = 0.000000000001;

        // cout << "reset kappa" << porosity[i][j] << endl; // just a help to
        // make the code more stable
      }
    }
  }
}

void Pressure_Lattice::alpha_beta_set(double time) {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      alpha[i][j] = 0;
      beta[i][j] = 0;

      {
        alpha[i][j] = (1.0 + (compressibility * pressure_dist[i][j]));
        alpha[i][j] *= kappa[i][j] * time_a;
        alpha[i][j] /= (2.0 * mu * compressibility * porosity[i][j] *
                        pow(delta_x * scale_constt, 2.0));

        beta[i][j] = (1.0 + (compressibility * pressure_dist[i][j]));
        beta[i][j] *= kappa[i][j] * time_a;
        beta[i][j] /= (2.0 * mu * compressibility * porosity[i][j] *
                       pow(delta_y * scale_constt, 2.0));
      }
    }
  }
}

void Pressure_Lattice::matrices_set_1() {
  // from 1st half-time step of ADI method
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      a[i][j] = 0;
      b[i][j] = 0;

      if (i == j) // for diagonal elements
      {
        if (i == 0 || i == dim - 1) // set boundary conditions
        {
          a[i][j] = 1.0; // column set for any given column j

          b[i][j] = 1.0; // row set for any given row i
        } else {
          a[i][j] = 1.0 + 2.0 * alpha[i][j];

          b[i][j] = 1.0 - 2.0 * beta[i][j];
        }
      } else if ((i == j - 1 && i > 0) || (i == j + 1 && i < (dim - 1))) {
        a[i][j] = -alpha[i][j];

        b[i][j] = beta[i][j];
      }
    }
  }
}

void Pressure_Lattice::matrices_set_2() {
  // from 2nd half-time step of ADI method
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      a[i][j] = 0;
      b[i][j] = 0;

      if (i == j) // for diagonal elements
      {
        if (i == 0 || i == dim - 1) // set boundary conditions
        {
          a[i][j] = 1.0; // column set for any given column j

          b[i][j] = 1.0; // row set for any given row i
        } else {
          a[i][j] = 1.0 + 2.0 * beta[i][j];

          b[i][j] = 1.0 - 2.0 * alpha[i][j];
        }
      } else if ((i == j - 1 && i > 0) || (i == j + 1 && i < (dim - 1))) {
        a[i][j] = -beta[i][j];

        b[i][j] = alpha[i][j];
      }
    }
  }
}

void Pressure_Lattice::matrix_a_1st() {
  // from 1st half-time step of ADI method
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      a[i][j] = 0;

      if (i == j) // for diagonal elements
      {
        if (i == 0 || i == dim - 1) // set boundary conditions
          a[i][j] = 1.0;            // column set for any given column j
        else
          a[i][j] = 1.0 + 2.0 * alpha[i][j];
      } else if ((i == j - 1 && i > 0) || (i == j + 1 && i < (dim - 1)))
        a[i][j] = -alpha[i][j];
    }
  }

  // for(i = 0; i < dim; i++)
  //{
  // for(j = 0; j < dim; j++)
  //{
  // cout << a[i][j] << " ";

  // if(j == dim-1)
  // cout << endl << i << endl;
  //}
  //}
}

void Pressure_Lattice::matrix_a_2_step() {
  // from 2nd half-time step of ADI method
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      a[i][j] = 0;

      if (i == j) // for diagonal elements
      {
        if (i == 0 || i == dim - 1) // set boundary conditions
          a[i][j] = 1.0;            // column set for any given column j
        else
          a[i][j] = 1 + 2 * beta[i][j];
      } else if ((i == j - 1 && i > 0) || (i == j + 1 && i < (dim - 1)))
        a[i][j] = -beta[i][j];
    }
  }
}

void Pressure_Lattice::invert() {
  double factor;

  // Creation of identity matrix
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      a_inv[i][j] = 0;

      a_inv[i][j] = (i == j) ? 1 : 0;
    }
  }

  // Development of upper-triangular matrix from the input matrix
  // along with the subsequent changes in identity matrix.
  for (i = 0; i < dim - 1; i++) {
    for (k = i + 1; k < dim; k++) {
      factor = (-a[k][i]) / a[i][i];
      if (!isinf(factor) && !isnan(factor)) {
        for (j = 0; j < dim; j++) {
          a[k][j] += a[i][j] * factor;
          a_inv[k][j] += a_inv[i][j] * factor;
        }
      }
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      if (j < i)
        a[i][j] = 0;
    }
  }

  // Development of diagonal matrix from developed upper-triangular
  // matrix along with subsequent changes in modified identity matrix.

  for (i = dim - 1; i > 0; i--) {
    for (k = i - 1; k >= 0; k--) {
      factor = (-a[k][i]) / a[i][i];
      if (!isinf(factor) && !isnan(factor)) {
        for (j = 0; j < dim; j++) {
          a[k][j] += a[i][j] * factor;
          a_inv[k][j] += a_inv[i][j] * factor;
        }
      }
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      if (j > i)
        a[i][j] = 0;
    }
  }

  // Conversion of given matrix into identity matrix and on contrary
  // achievement of inverse matrix from the identity matrix.

  for (i = 0; i < dim; i++) {
    if (a[i][i] != 0) {
      factor = a[i][i];

      for (j = 0; j < dim; j++) {
        a[i][j] /= factor;
        a_inv[i][j] /= factor;
      }
    }
  }
}

void Pressure_Lattice::multiplication_1st() {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      interim_pressure_dist[i][j] = 0;

      for (k = 0; k < dim; k++)
        interim_pressure_dist[i][j] += b[i][k] * pressure_dist[k][j];
    }
  }
}

void Pressure_Lattice::multiplication_x() {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      interim_pressure_dist[i][j] = pressure_dist[i][j];

      if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
        interim_pressure_dist[i][j] =
            (1 - 2 * beta[i][j]) * pressure_dist[i][j];
        interim_pressure_dist[i][j] += 2 * beta[i][j] * pressure_dist[i + 1][j];
        interim_pressure_dist[i][j] += 2 * beta[i][j] * pressure_dist[i - 1][j];
      }
    }
  }
}

void Pressure_Lattice::multiplication_y() {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      interim_pressure_dist[i][j] = pressure_dist[i][j];

      if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
        interim_pressure_dist[i][j] =
            (1 - 2 * alpha[j][i]) * pressure_dist[i][j];
        interim_pressure_dist[i][j] +=
            2 * alpha[j][i] * pressure_dist[i + 1][j];
        interim_pressure_dist[i][j] +=
            2 * alpha[j][i] * pressure_dist[i - 1][j];
      }
    }
  }
}

void Pressure_Lattice::multiplication_inv() {

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      pressure_dist[i][j] = 0;

      for (k = 0; k < dim; k++)
        pressure_dist[i][j] += a_inv[i][k] * interim_pressure_dist[k][j];
    }
  }
}

void Pressure_Lattice::transpose_P() {
  double pressure_trans[dim][dim];

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      pressure_trans[i][j] = 0;

      pressure_trans[i][j] = pressure_dist[j][i];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      pressure_dist[i][j] = 0;

      pressure_dist[i][j] = pressure_trans[i][j];
    }
  }

  // for(i = 0; i < dim; i++)
  //{
  // for(j = 0; j < dim; j++)
  //{
  // cout << pressure_dist[i][j] << " ";

  // if(j == dim-1)
  // cout << endl << i << endl;
  //}
  //}
}

void Pressure_Lattice::transpose_trans_P() {
  double pressure_trans[dim][dim];

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      pressure_trans[i][j] = 0;

      pressure_trans[i][j] = interim_pressure_dist[j][i];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      interim_pressure_dist[i][j] = 0;

      interim_pressure_dist[i][j] = pressure_trans[i][j];
    }
  }
}

void Pressure_Lattice::pressure_with_source(double time) {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      source_a[i][j] = 0;

      if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
        source_a[i][j] = ((vel_x[i][j + 1] - vel_x[i][j - 1]) /
                          (2.0 * delta_x * scale_constt));
        source_a[i][j] += ((vel_y[i + 1][j] - vel_y[i - 1][j]) /
                           (2.0 * delta_y * scale_constt));
        source_a[i][j] *= (1.0 / compressibility + sourced_pressure[i][j]);
        source_a[i][j] *= time_a / porosity[i][j];
      }
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      source_b[i][j] = 0;

      if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
        source_b[i][j] =
            vel_x[i][j] *
            ((sourced_pressure[i][j + 1] - sourced_pressure[i][j - 1]) /
             (2.0 * delta_x * scale_constt));
        source_b[i][j] +=
            vel_y[i][j] *
            ((sourced_pressure[i + 1][j] - sourced_pressure[i - 1][j]) /
             (2.0 * delta_y * scale_constt));
        source_b[i][j] *= time_a;
      }
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      source[i][j] = 0;

      source[i][j] = source_a[i][j] + source_b[i][j];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      pressure_dist[i][j] -= source[i][j];
    }
  }
}

void Pressure_Lattice::pressure_incre_source(
    int press_y, int press_x, double additional_press) //, int max_press)
{

  pressure_dist[press_y][press_x] += additional_press;
}

void Pressure_Lattice::Pressure_incre_twoXsources(int y_press, int x_press_diff,
                                                  double press_additional) {
  for (i = 1; i < dim - 1; i++) {
    for (j = 1; j < dim - 1; j++) {
      if (i == y_press &&
          (j == (dim / 2) - x_press_diff || j == (dim / 2) + x_press_diff)) {
        pressure_dist[i][j] += press_additional;
      }
    }
  }
}

void Pressure_Lattice::Pressure_incre_twoYsources(int x_press, int y_press_diff,
                                                  double press_additional) {
  for (i = 1; i < dim - 1; i++) {
    for (j = 1; j < dim - 1; j++) {
      if ((i == (dim / 2) - y_press_diff || i == (dim / 2) + y_press_diff) &&
          j == x_press) {
        pressure_dist[i][j] += press_additional;
        // cout << i << " " << j << endl;
      }
    }
  }
}

void Pressure_Lattice::Pressure_source_n_sink(int press_y_a, int press_x_a,
                                              int press_y_b, int press_x_b,
                                              double additional_press,
                                              double reduced_press) {
  pressure_dist[press_y_a][press_x_a] += additional_press;

  pressure_dist[press_y_b][press_x_b] -= reduced_press;
}

void Pressure_Lattice::pressure_point_sink(int press_y, int press_x,
                                           double reduced_press) {
  pressure_dist[press_y][press_x] -= reduced_press;
}

void Pressure_Lattice::pressure_sink(float y_press_min, float y_press_max,
                                     float x_press_min, float x_press_max,
                                     double reducing_press) {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      if (i >= dim * y_press_min && i <= dim * y_press_max) {
        if (j >= dim * x_press_min && j <= dim * x_press_max) {
          pressure_dist[i][j] = reducing_press;
        }
      }
    }
  }
}

void Pressure_Lattice::pressure_sink_cont(float y_press_min, float y_press_max,
                                          float x_press_min, float x_press_max,
                                          double reducing_press) {
  int y_min, y_max, x_min, x_max, med_y, med_x;
  float diff;

  y_min = dim * y_press_min;
  y_max = dim * y_press_max;
  x_min = dim * x_press_min;
  x_max = dim * x_press_max;

  med_y = (y_min + y_max) / 2;
  med_x = (x_min + x_max) / 2;

  diff = (y_press_max - y_press_min) * dim;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      if (i >= y_min && i <= y_max) {
        if (j >= x_min && j <= x_max) {
          // if(pressure_dist[med_y][med_x] > -1e+07)
          {
            if (pressure_dist[i][j] > -5 * (2 * reducing_press)) {
              for (k = 0; k < (diff / 2); k++) {
                if ((i > y_min + k && i < y_max - k) &&
                    (j > x_min + k && j < x_max - k))
                  pressure_dist[i][j] -= reducing_press / (diff / 2);
              }
            } else
              pressure_dist[i][j] += 0;
          }
          // else
          // pressure_dist[med_y][med_x] = -1e+07; //-1 * reducing_press;
        }
      }
    }
  }
}

void Pressure_Lattice::pressure_particles(Particle **list) {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
        for (k = 0; k < 4; k++) {
          if (k == 0)
            box = ((i - 1) * 400) + ((j - 1) * 2);
          else if (k == 1)
            box = ((i - 1) * 400) + 1 + ((j - 1) * 2);
          else if (k == 2)
            box = ((i - 1) * 400) + 200 + ((j - 1) * 2);
          else
            box = ((i - 1) * 400) + 201 + ((j - 1) * 2);

          if (list[box]) {
            help = list[box];

            x_diff = help->xpos - Pressure_x[i][j];
            y_diff = help->ypos - Pressure_y[i][j];
            if (x_diff < 0)
              x_diff = x_diff * -1;
            if (y_diff < 0)
              y_diff = y_diff * -1;
            smooth_func = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

            x_diff = help->xpos - Pressure_x[i - 1][j];
            y_diff = help->ypos - Pressure_y[i - 1][j];
            if (x_diff < 0)
              x_diff = x_diff * -1;
            if (y_diff < 0)
              y_diff = y_diff * -1;
            if (i - 1 > 0) {
              if (x_diff < delta_x && y_diff < delta_y)
                //{
                smooth_i_low =
                    (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

              // if(x_diff == 0 && y_diff == 0)
              // deta_i_low = 0;
              // else
              // deta_i_low = atan (y_diff/x_diff) * 180/PI;
              // sin_deta_i_low = sin (deta_i_low * PI/180);
              // cos_deta_i_low = cos (deta_i_low * PI/180);

              // dl_i_low = sqrt(pow(x_diff,2)+pow(y_diff,2));
              //}
              else
                //{
                smooth_i_low = 0;
              // sin_deta_i_low = 0;
              // cos_deta_i_low = 0;
              // dl_i_low = 0;
              //}
            } else
              //{
              smooth_i_low = 0;
            // sin_deta_i_low = 0;
            // cos_deta_i_low = 0;
            // dl_i_low = 0;
            //}

            x_diff = help->xpos - Pressure_x[i + 1][j];
            y_diff = help->ypos - Pressure_y[i + 1][j];
            // if(x_diff == 0 && y_diff == 0)
            // deta_i_up = 0;
            // else
            // deta_i_up = atan2 (y_diff,x_diff) * 180/PI;
            if (x_diff < 0)
              x_diff = x_diff * -1;
            if (y_diff < 0)
              y_diff = y_diff * -1;
            if (i + 1 < dim - 1) {
              if (x_diff < delta_x && y_diff < delta_y)
                //{
                smooth_i_up =
                    (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

              // if(x_diff == 0 && y_diff == 0)
              // deta_i_up = 0;
              // else
              // deta_i_up = atan (y_diff/x_diff) * 180/PI;
              // sin_deta_i_up = sin (deta_i_up * PI/180);
              // cos_deta_i_up = cos (deta_i_up * PI/180);

              // dl_i_up = sqrt(pow(x_diff,2)+pow(y_diff,2));
              //}
              else
                //{
                smooth_i_up = 0;
              // sin_deta_i_up = 0;
              // cos_deta_i_up = 0;
              // dl_i_up = 0;
              //}
            } else
              //{
              smooth_i_up = 0;
            // sin_deta_i_up = 0;
            // cos_deta_i_up = 0;
            // dl_i_up = 0;
            //}

            x_diff = help->xpos - Pressure_x[i][j - 1];
            y_diff = help->ypos - Pressure_y[i][j - 1];
            // if(x_diff == 0 && y_diff == 0)
            // deta_j_back = 0;
            // else
            // deta_j_back = atan2 (y_diff,x_diff) * 180/PI;
            if (x_diff < 0)
              x_diff = x_diff * -1;
            if (y_diff < 0)
              y_diff = y_diff * -1;
            if (j - 1 > 0) {
              if (x_diff < delta_x && y_diff < delta_y)
                //{
                smooth_j_back =
                    (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

              // if(x_diff == 0 && y_diff == 0)
              // deta_j_back = 0;
              // else
              // deta_j_back = atan (y_diff/x_diff) * 180/PI;
              // cos_deta_j_back = cos (deta_j_back * PI/180);
              // sin_deta_j_back = sin (deta_j_back * PI/180);

              // dl_j_back = sqrt(pow(x_diff,2)+pow(y_diff,2));
              //}
              else
                //{
                smooth_j_back = 0;
              // cos_deta_j_back = 0;
              // sin_deta_j_back = 0;
              // dl_j_back = 0;
              //}
            } else
              //{
              smooth_j_back = 0;
            // cos_deta_j_back = 0;
            // sin_deta_j_back = 0;
            // dl_j_back = 0;
            //}

            x_diff = help->xpos - Pressure_x[i][j + 1];
            y_diff = help->ypos - Pressure_y[i][j + 1];
            // if(x_diff == 0 && y_diff == 0)
            // deta_j_fwd = 0;
            // else
            // deta_j_fwd = atan2 (y_diff,x_diff) * 180/PI;
            if (x_diff < 0)
              x_diff = x_diff * -1;
            if (y_diff < 0)
              y_diff = y_diff * -1;
            if (j + 1 < dim - 1) {
              if (x_diff < delta_x && y_diff < delta_y)
                //{
                smooth_j_fwd =
                    (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

              // if(x_diff == 0 && y_diff == 0)
              // deta_j_fwd = 0;
              // else
              // deta_j_fwd = atan (y_diff/x_diff) * 180/PI;
              // cos_deta_j_fwd = cos (deta_j_fwd * PI/180);
              // sin_deta_j_fwd = sin (deta_j_fwd * PI/180);

              // dl_j_fwd = sqrt(pow(x_diff,2)+pow(y_diff,2));
              //}
              else
                //{
                smooth_j_fwd = 0;
              // cos_deta_j_fwd = 0;
              // sin_deta_j_fwd = 0;
              // dl_j_fwd = 0;
              //}
            } else
              //{
              smooth_j_fwd = 0;
            // cos_deta_j_fwd = 0;
            // sin_deta_j_fwd = 0;
            // dl_j_fwd = 0;
            //}

            if (smooth_i_low == 0 || rho[i - 1][j] == 0) {
              F_P_i_low = 0;
              F_P_j_low = 0;

              // F_P_ij_low = 0;

              // F_P_iy_low = 0;
              // F_P_ix_low = 0;
            } else {
              // delta_P_i_low = (pressure_dist[i-1][j] * sin_deta_i_low) -
              // (pressure_dist[i][j] * sin_deta); F_P_i_low = delta_P_i_low *
              // (smooth_i_low/rho[i-1][j]);

              // delta_P_ij_low = (pressure_dist[i-1][j] * cos_deta_i_low) -
              // (pressure_dist[i][j] * cos_deta); F_P_ij_low = delta_P_ij_low *
              // (smooth_i_low/rho[i-1][j]);

              // F_P_iy_low = pressure_dist[i-1][j] * sin_deta_i_low; // *
              // (smooth_i_low/rho[i-1][j]); F_P_ix_low = pressure_dist[i-1][j]
              // * cos_deta_i_low; // * (smooth_i_low/rho[i-1][j]);

              // delta_P_i_low = (pressure_dist[i][j]/dl) -
              // (pressure_dist[i-1][j]/dl_i_low); if(delta_P_i_low < 0)
              // F_P_i_low = -1 * delta_P_i_low * (smooth_i_low/rho[i-1][j]);
              // else
              // F_P_i_low = delta_P_i_low * (smooth_i_low/rho[i-1][j]);

              // F_P_i_low = -1
              // *((pressure_dist[i][j]-pressure_dist[i-1][j])/delta_y) *
              // (smooth_i_low/rho[i][j]);

              if (i - 2 > 0) {
                delta_P_i_low =
                    -1 * ((pressure_dist[i][j] - pressure_dist[i - 2][j]) /
                          (2 * delta_y));
                F_P_i_low = delta_P_i_low * (smooth_i_low / rho[i - 1][j]);
              } else
                F_P_i_low = 0;

              if (i - 1 > 0 && j - 1 > 0 && j + 1 < dim - 1) {
                delta_P_j_low = -1 * ((pressure_dist[i - 1][j + 1] -
                                       pressure_dist[i - 1][j - 1]) /
                                      (2 * delta_x));
                F_P_j_low = delta_P_j_low * (smooth_i_low / rho[i - 1][j]);
              } else
                F_P_j_low = 0;
              // if(i-1 > 1)
              // delta_P_i_low = -1 * ((pressure_dist[i][j] -
              // pressure_dist[i-2][j])/(2 * delta_y));
              // else
              // delta_P_i_low = 0;
              // delta_P_i_low = -1 * (pressure_dist[i][j]/(2 * delta_y));

              // F_P_i_low = smooth_i_low * (delta_P_i_low / rho[i-1][j]);
            }

            if (smooth_i_up == 0 || rho[i + 1][j] == 0) {
              F_P_i_up = 0;
              F_P_j_up = 0;
              // F_P_ij_up = 0;

              // F_P_iy_up = 0;
              // F_P_ix_up = 0;
            } else {
              // delta_P_i_up = (pressure_dist[i][j] * sin_deta) -
              // (pressure_dist[i+1][j] * sin_deta_i_up); F_P_i_up =
              // delta_P_i_up * (smooth_i_up/rho[i+1][j]);

              // delta_P_ij_up = (pressure_dist[i][j] * cos_deta) -
              // (pressure_dist[i+1][j] * cos_deta_i_up); F_P_ij_up =
              // delta_P_ij_up * (smooth_i_up/rho[i+1][j]);

              // F_P_iy_up = pressure_dist[i+1][j] * sin_deta_i_up; // *
              // (smooth_i_up/rho[i+1][j]); F_P_ix_up = pressure_dist[i+1][j] *
              // cos_deta_i_up; // * (smooth_i_up/rho[i+1][j]);

              // delta_P_i_up = (pressure_dist[i+1][j]/dl_i_up) -
              // (pressure_dist[i][j]/dl); if(delta_P_i_up < 0) F_P_i_up = -1 *
              // delta_P_i_up * (smooth_i_up/rho[i+1][j]);
              // else
              // F_P_i_up = delta_P_i_up * (smooth_i_up/rho[i+1][j]);

              // F_P_i_up = -1 *
              // ((pressure_dist[i+1][j]-pressure_dist[i][j])/delta_y) *
              // (smooth_i_up/rho[i][j]);

              if (i + 2 < dim - 1) {
                delta_P_i_up =
                    -1 * ((pressure_dist[i + 2][j] - pressure_dist[i][j]) /
                          (2 * delta_y));
                F_P_i_up = delta_P_i_up * (smooth_i_up / rho[i + 1][j]);
              } else
                F_P_i_up = 0;

              if (j - 1 > 0 && i + 1 < dim - 1 && j + 1 < dim - 1) {
                delta_P_j_up = -1 * ((pressure_dist[i + 1][j + 1] -
                                      pressure_dist[i + 1][j - 1]) /
                                     (2 * delta_x));
                F_P_j_up = delta_P_j_up * (smooth_i_up / rho[i + 1][j]);
              } else
                F_P_j_up = 0;

              // if(i+1 < dim-2)
              // delta_P_i_up = -1 * ((pressure_dist[i+2][j] -
              // pressure_dist[i][j])/(2 * delta_y));
              // else
              // delta_P_i_up = 0;
              ////delta_P_i_up = pressure_dist[i][j]/(2 * delta_y);

              // F_P_i_up = smooth_i_up * (delta_P_i_up / rho[i+1][j]);
            }

            if (smooth_j_back == 0 || rho[i][j - 1] == 0) {
              F_P_i_back = 0;
              F_P_j_back = 0;
              // F_P_ij_back = 0;

              // F_P_jx_back = 0;
              // F_P_jy_back = 0;
            } else {
              // delta_P_j_back = (pressure_dist[i][j-1] * cos_deta_j_back) -
              // (pressure_dist[i][j] * cos_deta); F_P_j_back = delta_P_j_back *
              // (smooth_j_back/rho[i][j-1]);

              // delta_P_ij_back = (pressure_dist[i][j-1] * sin_deta_j_back) -
              // (pressure_dist[i][j] * sin_deta); F_P_ij_back = delta_P_ij_back
              // * (smooth_j_back/rho[i][j-1]);

              // F_P_jx_back = pressure_dist[i][j-1] * cos_deta_j_back; // *
              // (smooth_j_back/rho[i][j-1]); F_P_jy_back =
              // pressure_dist[i][j-1] * sin_deta_j_back; // *
              // (smooth_j_back/rho[i][j-1]);

              // delta_P_j_back = (pressure_dist[i][j]/dl) -
              // (pressure_dist[i][j-1]/dl_j_back); if(delta_P_j_back < 0)
              // F_P_j_back = -1 * delta_P_j_back * (smooth_j_back/rho[i][j-1]);
              // else
              // F_P_j_back = delta_P_j_back * (smooth_j_back/rho[i][j-1]);

              // F_P_j_back = -1 *
              // ((pressure_dist[i][j]-pressure_dist[i][j-1])/delta_x) *
              // (smooth_j_back/rho[i][j]);

              if (i - 1 > 0 && j - 1 > 0 && i + 1 < dim - 1) {
                delta_P_i_back = -1 * ((pressure_dist[i + 1][j - 1] -
                                        pressure_dist[i - 1][j - 1]) /
                                       (2 * delta_y));
                F_P_i_back = delta_P_i_back * (smooth_j_back / rho[i][j - 1]);
              } else
                F_P_i_back = 0;

              if (j - 2 > 0) {
                delta_P_j_back =
                    -1 * ((pressure_dist[i][j] - pressure_dist[i][j - 2]) /
                          (2 * delta_x));
                F_P_j_back = delta_P_j_back * (smooth_j_back / rho[i][j - 1]);
              } else
                F_P_j_back = 0;

              // if(j-1 > 1)
              // delta_P_j_back = -1 * ((pressure_dist[i][j] -
              // pressure_dist[i][j-2])/(2 * delta_x));
              // else
              // delta_P_j_back = 0;
              ////delta_P_j_back = -1 * (pressure_dist[i][j]/(2 * delta_x));

              // F_P_j_back = smooth_j_back * (delta_P_j_back / rho[i][j-1]);
            }

            if (smooth_j_fwd == 0 || rho[i][j + 1] == 0) {
              F_P_i_fwd = 0;
              F_P_j_fwd = 0;
              // F_P_ij_fwd = 0;

              // F_P_jx_fwd = 0;
              // F_P_jy_fwd = 0;
            } else {
              // delta_P_j_fwd = (pressure_dist[i][j] * cos_deta) -
              // (pressure_dist[i][j+1] * cos_deta_j_fwd); F_P_j_fwd =
              // delta_P_j_fwd * (smooth_j_fwd/rho[i][j+1]);

              // delta_P_ij_fwd = (pressure_dist[i][j] * sin_deta) -
              // (pressure_dist[i][j+1] * sin_deta_j_fwd); F_P_ij_fwd =
              // delta_P_ij_fwd * (smooth_j_fwd/rho[i][j+1]);

              // F_P_jx_fwd = pressure_dist[i][j+1] * cos_deta_j_fwd; // *
              // (smooth_j_fwd/rho[i][j+1]); F_P_jy_fwd = pressure_dist[i][j+1]
              // * sin_deta_j_fwd; // * (smooth_j_fwd/rho[i][j+1]);

              // delta_P_j_fwd = (pressure_dist[i][j+1]/dl_j_fwd) -
              // (pressure_dist[i][j]/dl); if(delta_P_j_fwd < 0) F_P_j_fwd = -1
              // * delta_P_j_fwd * (smooth_j_fwd/rho[i][j+1]);
              // else
              // F_P_j_fwd = delta_P_j_fwd * (smooth_j_fwd/rho[i][j+1]);

              // F_P_j_fwd = -1 *
              // ((pressure_dist[i][j+1]-pressure_dist[i][j])/delta_x) *
              // (smooth_j_fwd/rho[i][j]);

              if (i - 1 > 0 && i + 1 < dim - 1 && j + 1 < dim - 1) {
                delta_P_i_fwd = -1 * ((pressure_dist[i + 1][j + 1] -
                                       pressure_dist[i - 1][j + 1]) /
                                      (2 * delta_y));
                F_P_i_fwd = delta_P_i_fwd * (smooth_j_fwd / rho[i][j + 1]);
              } else
                F_P_i_fwd = 0;

              if (j + 2 < dim - 1) {
                delta_P_j_fwd =
                    -1 * ((pressure_dist[i][j + 2] - pressure_dist[i][j]) /
                          (2 * delta_x));
                F_P_j_fwd = delta_P_j_fwd * (smooth_j_fwd / rho[i][j + 1]);
              } else
                F_P_j_fwd = 0;

              // if(j+1 < dim-2)
              // delta_P_j_fwd = -1 * ((pressure_dist[i][j+2] -
              // pressure_dist[i][j])/(2 * delta_x));
              // else
              // delta_P_j_fwd = 0;
              // delta_P_j_fwd = pressure_dist[i][j]/(2 * delta_x);

              // F_P_j_fwd = smooth_j_fwd * (delta_P_j_fwd / rho[i][j+1]);
            }

            if (j - 1 > 0 && j + 1 < dim - 1)
              F_P_ij_x = -1 *
                         ((pressure_dist[i][j + 1] - pressure_dist[i][j - 1]) /
                          (2 * delta_x)) *
                         (smooth_func / rho[i][j]);
            else
              F_P_ij_x = 0;

            if (i - 1 > 0 && i + 1 < dim - 1)
              F_P_ij_y = -1 *
                         ((pressure_dist[i + 1][j] - pressure_dist[i - 1][j]) /
                          (2 * delta_y)) *
                         (smooth_func / rho[i][j]);
            else
              F_P_ij_y = 0;

            // F_P_ij_x = -1 * (smooth_func *
            // (pressure_dist[i][j+1]-pressure_dist[i][j-1]))/(rho[i][j]*2*delta_x);
            // F_P_ij_y = -1 * (smooth_func *
            // (pressure_dist[i+1][j]-pressure_dist[i-1][j]))/(rho[i][j]*2*delta_y);

            help->F_P_x = F_P_ij_x;
            help->F_P_y = F_P_ij_y;

            // deleted next in box pointer below for both lines Daniel June 2013

            help->fluid_pressure_gradient =
                sqrt((F_P_ij_x * F_P_ij_x) + (F_P_ij_y * F_P_ij_y));
            help->average_porosity = porosity[i][j];

            // cout << F_P_j_back << " ";

            // help->F_P_x = F_P_i_low + F_P_i_up + F_P_i_back + F_P_i_fwd;
            // help->F_P_y = F_P_j_low + F_P_j_up + F_P_j_back + F_P_j_fwd;

            help->temperature = pressure_dist[i][j];

            // help->F_P_x = (F_P_j_back + F_P_j_fwd)/2;
            // help->F_P_y = (F_P_i_low + F_P_i_up)/2;

            // help->F_P_x = (F_P_j_back + F_P_j_fwd + F_P_ij_low + F_P_ij_up);
            // help->F_P_y = (F_P_i_low + F_P_i_up + F_P_ij_back + F_P_ij_fwd);

            // help->F_P_x = (F_P_ijx + F_P_jx_back + F_P_jx_fwd + F_P_ix_low +
            // F_P_ix_up); help->F_P_y = (F_P_ijy + F_P_jy_back + F_P_jy_fwd +
            // F_P_iy_low + F_P_iy_up);

            // help->sheet_young = F_P_x + F_P_y;

            // cout << help->F_P_y << " ";

            while (help->next_inBox) {
              x_diff = help->next_inBox->xpos - Pressure_x[i][j];
              y_diff = help->next_inBox->ypos - Pressure_y[i][j];
              // if(x_diff == 0 && y_diff == 0)
              // deta = 0;
              // else
              // deta = atan2 (y_diff,x_diff) * 180/PI;
              if (x_diff < 0)
                x_diff = x_diff * -1;
              if (y_diff < 0)
                y_diff = y_diff * -1;
              smooth_func = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));
              // if(x_diff == 0 && y_diff == 0)
              // deta = 0;
              // else
              // deta = atan (y_diff/x_diff) * 180/PI;
              // sin_deta = sin (deta * PI/180);
              // cos_deta = cos (deta * PI/180);
              // dl = sqrt(pow(x_diff,2)+pow(y_diff,2));

              x_diff = help->xpos - Pressure_x[i - 1][j];
              y_diff = help->ypos - Pressure_y[i - 1][j];
              // if(x_diff == 0 && y_diff == 0)
              // deta_i_low = 0;
              // else
              // deta_i_low = atan2 (y_diff,x_diff) * 180/PI;
              if (x_diff < 0)
                x_diff = x_diff * -1;
              if (y_diff < 0)
                y_diff = y_diff * -1;
              if (i - 1 > 0) {
                if (x_diff < delta_x && y_diff < delta_y)
                  //{
                  smooth_i_low =
                      (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

                // if(x_diff == 0 && y_diff == 0)
                // deta_i_low = 0;
                // else
                // deta_i_low = atan (y_diff/x_diff) * 180/PI;
                // sin_deta_i_low = sin (deta_i_low * PI/180);
                // cos_deta_i_low = cos (deta_i_low * PI/180);

                // dl_i_low = sqrt(pow(x_diff,2)+pow(y_diff,2));
                //}
                else
                  //{
                  smooth_i_low = 0;
                // sin_deta_i_low = 0;
                // cos_deta_i_low = 0;
                // dl_i_low = 0;
                //}
              } else
                //{
                smooth_i_low = 0;
              // sin_deta_i_low = 0;
              // cos_deta_i_low = 0;
              // dl_i_low = 0;
              //}

              x_diff = help->xpos - Pressure_x[i + 1][j];
              y_diff = help->ypos - Pressure_y[i + 1][j];
              // if(x_diff == 0 && y_diff == 0)
              // deta_i_up = 0;
              // else
              // deta_i_up = atan2 (y_diff,x_diff) * 180/PI;
              if (x_diff < 0)
                x_diff = x_diff * -1;
              if (y_diff < 0)
                y_diff = y_diff * -1;
              if (i + 1 < dim - 1) {
                if (x_diff < delta_x && y_diff < delta_y)
                  //{
                  smooth_i_up =
                      (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

                // if(x_diff == 0 && y_diff == 0)
                // deta_i_up = 0;
                // else
                // deta_i_up = atan (y_diff/x_diff) * 180/PI;
                // sin_deta_i_up = sin (deta_i_up * PI/180);
                // cos_deta_i_up = cos (deta_i_up * PI/180);

                // dl_i_up = sqrt(pow(x_diff,2)+pow(y_diff,2));
                //}
                else
                  //{
                  smooth_i_up = 0;
                // sin_deta_i_up = 0;
                // cos_deta_i_up = 0;
                // dl_i_up = 0;
                //}
              } else
                //{
                smooth_i_up = 0;
              // sin_deta_i_up = 0;
              // cos_deta_i_up = 0;
              // dl_i_up = 0;
              //}

              x_diff = help->xpos - Pressure_x[i][j - 1];
              y_diff = help->ypos - Pressure_y[i][j - 1];
              // if(x_diff == 0 && y_diff == 0)
              // deta_j_back = 0;
              // else
              // deta_j_back = atan2 (y_diff,x_diff) * 180/PI;
              if (x_diff < 0)
                x_diff = x_diff * -1;
              if (y_diff < 0)
                y_diff = y_diff * -1;
              if (j - 1 > 0) {
                if (x_diff < delta_x && y_diff < delta_y)
                  //{
                  smooth_j_back =
                      (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

                // if(x_diff == 0 && y_diff == 0)
                // deta_j_back = 0;
                // else
                // deta_j_back = atan (y_diff/x_diff) * 180/PI;
                // cos_deta_j_back = cos (deta_j_back * PI/180);
                // sin_deta_j_back = sin (deta_j_back * PI/180);

                // dl_j_back = sqrt(pow(x_diff,2)+pow(y_diff,2));
                //}
                else
                  //{
                  smooth_j_back = 0;
                // cos_deta_j_back = 0;
                // sin_deta_j_back = 0;
                // dl_j_back = 0;
                //}
              } else
                //{
                smooth_j_back = 0;
              // cos_deta_j_back = 0;
              // sin_deta_j_back = 0;
              // dl_j_back = 0;
              //}

              x_diff = help->xpos - Pressure_x[i][j + 1];
              y_diff = help->ypos - Pressure_y[i][j + 1];
              // if(x_diff == 0 && y_diff == 0)
              // deta_j_fwd = 0;
              // else
              // deta_j_fwd = atan2 (y_diff,x_diff) * 180/PI;
              if (x_diff < 0)
                x_diff = x_diff * -1;
              if (y_diff < 0)
                y_diff = y_diff * -1;
              if (j + 1 < dim - 1) {
                if (x_diff < delta_x && y_diff < delta_y)
                  //{
                  smooth_j_fwd =
                      (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

                // if(x_diff == 0 && y_diff == 0)
                // deta_j_fwd = 0;
                // else
                // deta_j_fwd = atan (y_diff/x_diff) * 180/PI;
                // cos_deta_j_fwd = cos (deta_j_fwd * PI/180);
                // sin_deta_j_fwd = sin (deta_j_fwd * PI/180);

                // dl_j_fwd = sqrt(pow(x_diff,2)+pow(y_diff,2));
                //}
                else
                  //{
                  smooth_j_fwd = 0;
                // cos_deta_j_fwd = 0;
                // sin_deta_j_fwd = 0;
                // dl_j_fwd = 0;
                //}
              } else
                //{
                smooth_j_fwd = 0;
              // cos_deta_j_fwd = 0;
              // sin_deta_j_fwd = 0;
              // dl_j_fwd = 0;
              //}

              if (smooth_i_low == 0 || rho[i - 1][j] == 0) {
                F_P_i_low = 0;
                F_P_j_low = 0;

                // F_P_ij_low = 0;

                // F_P_iy_low = 0;
                // F_P_ix_low = 0;
              } else {
                // delta_P_i_low = (pressure_dist[i-1][j] * sin_deta_i_low) -
                // (pressure_dist[i][j] * sin_deta); F_P_i_low = delta_P_i_low *
                // (smooth_i_low/rho[i-1][j]);

                // delta_P_ij_low = (pressure_dist[i-1][j] * cos_deta_i_low) -
                // (pressure_dist[i][j] * cos_deta); F_P_ij_low = delta_P_ij_low
                // * (smooth_i_low/rho[i-1][j]);

                // F_P_iy_low = pressure_dist[i-1][j] * sin_deta_i_low; // *
                // (smooth_i_low/rho[i-1][j]); F_P_ix_low =
                // pressure_dist[i-1][j] * cos_deta_i_low; // *
                // (smooth_i_low/rho[i-1][j]);

                // delta_P_i_low = (pressure_dist[i][j]/dl) -
                // (pressure_dist[i-1][j]/dl_i_low); if(delta_P_i_low < 0)
                // F_P_i_low = -1 * delta_P_i_low * (smooth_i_low/rho[i-1][j]);
                // else
                // F_P_i_low = delta_P_i_low * (smooth_i_low/rho[i-1][j]);

                // F_P_i_low = -1
                // *((pressure_dist[i][j]-pressure_dist[i-1][j])/delta_y) *
                // (smooth_i_low/rho[i][j]);

                if (i - 2 > 0) {
                  delta_P_i_low =
                      -1 * ((pressure_dist[i][j] - pressure_dist[i - 2][j]) /
                            (2 * delta_y));
                  F_P_i_low = delta_P_i_low * (smooth_i_low / rho[i - 1][j]);
                } else
                  F_P_i_low = 0;

                if (i - 1 > 0 && j - 1 > 0 && j + 1 < dim - 1) {
                  delta_P_j_low = -1 * ((pressure_dist[i - 1][j + 1] -
                                         pressure_dist[i - 1][j - 1]) /
                                        (2 * delta_x));
                  F_P_j_low = delta_P_j_low * (smooth_i_low / rho[i - 1][j]);
                } else
                  F_P_j_low = 0;
                // if(i-1 > 1)
                // delta_P_i_low = -1 * ((pressure_dist[i][j] -
                // pressure_dist[i-2][j])/(2 * delta_y));
                // else
                // delta_P_i_low = 0;
                // delta_P_i_low = -1 * (pressure_dist[i][j]/(2 * delta_y));

                // F_P_i_low = smooth_i_low * (delta_P_i_low / rho[i-1][j]);
              }

              if (smooth_i_up == 0 || rho[i + 1][j] == 0) {
                F_P_i_up = 0;
                F_P_j_up = 0;
                // F_P_ij_up = 0;

                // F_P_iy_up = 0;
                // F_P_ix_up = 0;
              } else {
                // delta_P_i_up = (pressure_dist[i][j] * sin_deta) -
                // (pressure_dist[i+1][j] * sin_deta_i_up); F_P_i_up =
                // delta_P_i_up * (smooth_i_up/rho[i+1][j]);

                // delta_P_ij_up = (pressure_dist[i][j] * cos_deta) -
                // (pressure_dist[i+1][j] * cos_deta_i_up); F_P_ij_up =
                // delta_P_ij_up * (smooth_i_up/rho[i+1][j]);

                // F_P_iy_up = pressure_dist[i+1][j] * sin_deta_i_up; // *
                // (smooth_i_up/rho[i+1][j]); F_P_ix_up = pressure_dist[i+1][j]
                // * cos_deta_i_up; // * (smooth_i_up/rho[i+1][j]);

                // delta_P_i_up = (pressure_dist[i+1][j]/dl_i_up) -
                // (pressure_dist[i][j]/dl); if(delta_P_i_up < 0) F_P_i_up = -1
                // * delta_P_i_up * (smooth_i_up/rho[i+1][j]);
                // else
                // F_P_i_up = delta_P_i_up * (smooth_i_up/rho[i+1][j]);

                // F_P_i_up = -1 *
                // ((pressure_dist[i+1][j]-pressure_dist[i][j])/delta_y) *
                // (smooth_i_up/rho[i][j]);

                if (i + 2 < dim - 1) {
                  delta_P_i_up =
                      -1 * ((pressure_dist[i + 2][j] - pressure_dist[i][j]) /
                            (2 * delta_y));
                  F_P_i_up = delta_P_i_up * (smooth_i_up / rho[i + 1][j]);
                } else
                  F_P_i_up = 0;

                if (j - 1 > 0 && i + 1 < dim - 1 && j + 1 < dim - 1) {
                  delta_P_j_up = -1 * ((pressure_dist[i + 1][j + 1] -
                                        pressure_dist[i + 1][j - 1]) /
                                       (2 * delta_x));
                  F_P_j_up = delta_P_j_up * (smooth_i_up / rho[i + 1][j]);
                } else
                  F_P_j_up = 0;

                // if(i+1 < dim-2)
                // delta_P_i_up = -1 * ((pressure_dist[i+2][j] -
                // pressure_dist[i][j])/(2 * delta_y));
                // else
                // delta_P_i_up = 0;
                ////delta_P_i_up = pressure_dist[i][j]/(2 * delta_y);

                // F_P_i_up = smooth_i_up * (delta_P_i_up / rho[i+1][j]);
              }

              if (smooth_j_back == 0 || rho[i][j - 1] == 0) {
                F_P_i_back = 0;
                F_P_j_back = 0;
                // F_P_ij_back = 0;

                // F_P_jx_back = 0;
                // F_P_jy_back = 0;
              } else {
                // delta_P_j_back = (pressure_dist[i][j-1] * cos_deta_j_back) -
                // (pressure_dist[i][j] * cos_deta); F_P_j_back = delta_P_j_back
                // * (smooth_j_back/rho[i][j-1]);

                // delta_P_ij_back = (pressure_dist[i][j-1] * sin_deta_j_back) -
                // (pressure_dist[i][j] * sin_deta); F_P_ij_back =
                // delta_P_ij_back * (smooth_j_back/rho[i][j-1]);

                // F_P_jx_back = pressure_dist[i][j-1] * cos_deta_j_back; // *
                // (smooth_j_back/rho[i][j-1]); F_P_jy_back =
                // pressure_dist[i][j-1] * sin_deta_j_back; // *
                // (smooth_j_back/rho[i][j-1]);

                // delta_P_j_back = (pressure_dist[i][j]/dl) -
                // (pressure_dist[i][j-1]/dl_j_back); if(delta_P_j_back < 0)
                // F_P_j_back = -1 * delta_P_j_back *
                // (smooth_j_back/rho[i][j-1]);
                // else
                // F_P_j_back = delta_P_j_back * (smooth_j_back/rho[i][j-1]);

                // F_P_j_back = -1 *
                // ((pressure_dist[i][j]-pressure_dist[i][j-1])/delta_x) *
                // (smooth_j_back/rho[i][j]);

                if (i - 1 > 0 && j - 1 > 0 && i + 1 < dim - 1) {
                  delta_P_i_back = -1 * ((pressure_dist[i + 1][j - 1] -
                                          pressure_dist[i - 1][j - 1]) /
                                         (2 * delta_y));
                  F_P_i_back = delta_P_i_back * (smooth_j_back / rho[i][j - 1]);
                } else
                  F_P_i_back = 0;

                if (j - 2 > 0) {
                  delta_P_j_back =
                      -1 * ((pressure_dist[i][j] - pressure_dist[i][j - 2]) /
                            (2 * delta_x));
                  F_P_j_back = delta_P_j_back * (smooth_j_back / rho[i][j - 1]);
                } else
                  F_P_j_back = 0;

                // if(j-1 > 1)
                // delta_P_j_back = -1 * ((pressure_dist[i][j] -
                // pressure_dist[i][j-2])/(2 * delta_x));
                // else
                // delta_P_j_back = 0;
                ////delta_P_j_back = -1 * (pressure_dist[i][j]/(2 * delta_x));

                // F_P_j_back = smooth_j_back * (delta_P_j_back / rho[i][j-1]);
              }

              if (smooth_j_fwd == 0 || rho[i][j + 1] == 0) {
                F_P_i_fwd = 0;
                F_P_j_fwd = 0;
                // F_P_ij_fwd = 0;

                // F_P_jx_fwd = 0;
                // F_P_jy_fwd = 0;
              } else {
                // delta_P_j_fwd = (pressure_dist[i][j] * cos_deta) -
                // (pressure_dist[i][j+1] * cos_deta_j_fwd); F_P_j_fwd =
                // delta_P_j_fwd * (smooth_j_fwd/rho[i][j+1]);

                // delta_P_ij_fwd = (pressure_dist[i][j] * sin_deta) -
                // (pressure_dist[i][j+1] * sin_deta_j_fwd); F_P_ij_fwd =
                // delta_P_ij_fwd * (smooth_j_fwd/rho[i][j+1]);

                // F_P_jx_fwd = pressure_dist[i][j+1] * cos_deta_j_fwd; // *
                // (smooth_j_fwd/rho[i][j+1]); F_P_jy_fwd =
                // pressure_dist[i][j+1] * sin_deta_j_fwd; // *
                // (smooth_j_fwd/rho[i][j+1]);

                // delta_P_j_fwd = (pressure_dist[i][j+1]/dl_j_fwd) -
                // (pressure_dist[i][j]/dl); if(delta_P_j_fwd < 0) F_P_j_fwd =
                // -1 * delta_P_j_fwd * (smooth_j_fwd/rho[i][j+1]);
                // else
                // F_P_j_fwd = delta_P_j_fwd * (smooth_j_fwd/rho[i][j+1]);

                // F_P_j_fwd = -1 *
                // ((pressure_dist[i][j+1]-pressure_dist[i][j])/delta_x) *
                // (smooth_j_fwd/rho[i][j]);

                if (i - 1 > 0 && i + 1 < dim - 1 && j + 1 < dim - 1) {
                  delta_P_i_fwd = -1 * ((pressure_dist[i + 1][j + 1] -
                                         pressure_dist[i - 1][j + 1]) /
                                        (2 * delta_y));
                  F_P_i_fwd = delta_P_i_fwd * (smooth_j_fwd / rho[i][j + 1]);
                } else
                  F_P_i_fwd = 0;

                if (j + 2 < dim - 1) {
                  delta_P_j_fwd =
                      -1 * ((pressure_dist[i][j + 2] - pressure_dist[i][j]) /
                            (2 * delta_x));
                  F_P_j_fwd = delta_P_j_fwd * (smooth_j_fwd / rho[i][j + 1]);
                } else
                  F_P_j_fwd = 0;

                // if(j+1 < dim-2)
                // delta_P_j_fwd = -1 * ((pressure_dist[i][j+2] -
                // pressure_dist[i][j])/(2 * delta_x));
                // else
                // delta_P_j_fwd = 0;
                // delta_P_j_fwd = pressure_dist[i][j]/(2 * delta_x);

                // F_P_j_fwd = smooth_j_fwd * (delta_P_j_fwd / rho[i][j+1]);
              }

              if (j - 1 > 0 && j + 1 < dim - 1)
                F_P_ij_x =
                    -1 *
                    ((pressure_dist[i][j + 1] - pressure_dist[i][j - 1]) /
                     (2 * delta_x)) *
                    (smooth_func / rho[i][j]);
              else
                F_P_ij_x = 0;

              if (i - 1 > 0 && i + 1 < dim - 1)
                F_P_ij_y =
                    -1 *
                    ((pressure_dist[i + 1][j] - pressure_dist[i - 1][j]) /
                     (2 * delta_y)) *
                    (smooth_func / rho[i][j]);
              else
                F_P_ij_y = 0;

              // F_P_ij_x = -1 * (smooth_func *
              // (pressure_dist[i][j+1]-pressure_dist[i][j-1]))/(rho[i][j]*2*delta_x);
              // F_P_ij_y = -1 * (smooth_func *
              // (pressure_dist[i+1][j]-pressure_dist[i-1][j]))/(rho[i][j]*2*delta_y);

              help->next_inBox->F_P_x = F_P_ij_x;
              help->next_inBox->F_P_y = F_P_ij_y;

              help->next_inBox->fluid_pressure_gradient =
                  sqrt((F_P_ij_x * F_P_ij_x) + (F_P_ij_y * F_P_ij_y));
              help->next_inBox->average_porosity = rho[i][j];

              // cout << F_P_j_back << " ";

              // help->next_inBox->F_P_x = F_P_i_low + F_P_i_up + F_P_i_back +
              // F_P_i_fwd; help->next_inBox->F_P_y = F_P_j_low + F_P_j_up +
              // F_P_j_back + F_P_j_fwd;

              help->next_inBox->temperature = pressure_dist[i][j];

              // help->next_inBox->F_P_x = (F_P_j_back + F_P_j_fwd)/2;
              // help->next_inBox->F_P_y = (F_P_i_low + F_P_i_up)/2;

              // help->next_inBox->F_P_x = (F_P_j_back + F_P_j_fwd + F_P_ij_low
              // + F_P_ij_up); help->next_inBox->F_P_y = (F_P_i_low + F_P_i_up +
              // F_P_ij_back + F_P_ij_fwd);

              // help->next_inBox->F_P_x = (F_P_ijx + F_P_jx_back + F_P_jx_fwd +
              // F_P_ix_low + F_P_ix_up); help->next_inBox->F_P_y = (F_P_ijy +
              // F_P_jy_back + F_P_jy_fwd + F_P_iy_low + F_P_iy_up);

              // help->next_inBox->sheet_young = F_P_x + F_P_y;

              // cout << help->next_inBox->F_P_y << " ";

              help = help->next_inBox;
            }
          }
        }
      }

      // cout << pressure_dist[i][j] << " ";
      // if(j==dim-2)
      // cout << endl << i << endl;
    }
  }
}

void Pressure_Lattice::d_pressure_particles(Particle **list) {

  /******************************************************
   *  No extension (keep the dimensions 53*53)
   *
   * each inner node covers 16 repBox in its surroundings.
   *
   * Edged nodes already considered as boundary and dont cover
   * any rep_Box individually. here they just used to contribute
   * in pressure gradient calculation to the inner nodes, in the way
   *
   * delta_P_x[i][j] = -1 * (P_x[i][j+1] - P_x[i][j-1]) / (2*delta_x)
   * delta_P_y[i][j] = -1 * (P_y[i+1][j] - P_y[i-1][j]) / (2*delta_y)
   *
   ******************************************************/

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
        if (i == 1 && j == 1) {
          d_1x4_Boxes(list, 1, i, j);

        }

        else if (i == 1 && (j > 1 && j < dim - 2)) {
          d_2x4_Boxes(list, 1, i, j);

        }

        else if (i == 1 && j == dim - 2) {
          d_1x4_Boxes(list, 2, i, j);

        }

        else if ((i > 1 && i < dim - 2) && j == 1) {
          d_2x4_Boxes(list, 2, i, j);

        }

        else if ((i > 1 && i < dim - 2) && j == dim - 2) {
          d_2x4_Boxes(list, 3, i, j);

        }

        else if (i == dim - 2 && j == 1) {
          d_1x4_Boxes(list, 3, i, j);

        }

        else if (i == dim - 2 && (j > 1 && j < dim - 2)) {
          d_2x4_Boxes(list, 4, i, j);

        }

        else if (i == dim - 2 && j == dim - 2) {
          d_1x4_Boxes(list, 4, i, j);

        }

        else {
          d_4x4_Boxes(list, 1, i, j);
          d_4x4_Boxes(list, 2, i, j);
          d_4x4_Boxes(list, 3, i, j);
          d_4x4_Boxes(list, 4, i, j);
        }
      }
    }
  }
}

void Pressure_Lattice::d_1x4_Boxes(Particle **list, int position, int l,
                                   int m) {
  switch (position) {
  case 1: {
    //~ cout << "i " << i << " "<< " j " << j << endl;
    //~ cout << "l " << l << " "<< " m " << m << endl;
    for (k = 0; k < 4; k++) {
      d_UppRight_4d_Boxes(k, l, m);
      d_d_Boxes(list, l, m);
    }
    //~ cout << endl ;
  } break;

  case 2: {
    //~ cout << "i " << i << " "<< " j " << j << endl;
    for (k = 0; k < 4; k++) {
      d_UppLeft_4c_Boxes(k, l, m);
      d_c_Boxes(list, l, m);
    }
    //~ cout << endl ;
  } break;

  case 3: {
    //~ cout << "i " << i << " "<< " j " << j << endl;
    for (k = 0; k < 4; k++) {
      d_LowRight_4b_Boxes(k, l, m);
      d_b_Boxes(list, l, m);
    }
    //~ cout << endl ;
  } break;

  case 4: {
    //~ cout << "i " << i << " "<< " j " << j << endl;
    for (k = 0; k < 4; k++) {
      d_LowLeft_4a_Boxes(k, l, m);
      d_a_Boxes(list, l, m);
    }
    //~ cout << endl ;
  } break;

  default:
    cout << "No box in d_1x4_Boxes" << endl;
  }
}

void Pressure_Lattice::d_2x4_Boxes(Particle **list, int position, int l,
                                   int m) {
  switch (position) {
  case 1: {
    for (k = 0; k < 4; k++) {
      d_UppLeft_4c_Boxes(k, l, m);
      d_c_Boxes(list, l, m);
    }

    for (k = 0; k < 4; k++) {
      d_UppRight_4d_Boxes(k, l, m);
      d_d_Boxes(list, l, m);
    }
  } break;

  case 2: {
    for (k = 0; k < 4; k++) {
      d_LowRight_4b_Boxes(k, l, m);
      d_b_Boxes(list, l, m);
    }

    for (k = 0; k < 4; k++) {
      d_UppRight_4d_Boxes(k, l, m);
      d_d_Boxes(list, l, m);
    }
  } break;

  case 3: {
    for (k = 0; k < 4; k++) {
      d_LowLeft_4a_Boxes(k, l, m);

      d_a_Boxes(list, l, m);
    }

    for (k = 0; k < 4; k++) {
      d_UppLeft_4c_Boxes(k, l, m);
      d_c_Boxes(list, l, m);
    }
  } break;

  case 4: {
    for (k = 0; k < 4; k++) {
      d_LowLeft_4a_Boxes(k, l, m);
      d_a_Boxes(list, l, m);
    }

    for (k = 0; k < 4; k++) {
      d_LowRight_4b_Boxes(k, l, m);
      d_b_Boxes(list, l, m);
    }
  } break;

  default:
    cout << "No box in d_2x4_Boxes" << endl;
  }
}

void Pressure_Lattice::d_4x4_Boxes(Particle **list, int position, int l,
                                   int m) {
  switch (position) {
  case 1: {
    for (k = 0; k < 4; k++) {
      d_LowLeft_4a_Boxes(k, l, m);
      d_a_Boxes(list, l, m);
    }
  } break;

  case 2: {
    for (k = 0; k < 4; k++) {
      d_LowRight_4b_Boxes(k, l, m);
      d_b_Boxes(list, l, m);
    }
  } break;

  case 3: {
    for (k = 0; k < 4; k++) {
      d_UppLeft_4c_Boxes(k, l, m);
      d_c_Boxes(list, l, m);
    }
  } break;

  case 4: {
    for (k = 0; k < 4; k++) {
      d_UppRight_4d_Boxes(k, l, m);
      d_d_Boxes(list, l, m);
    }
  } break;

  default:
    cout << "No box in d_4x4_Boxes" << endl;
  }
}

void Pressure_Lattice::d_LowLeft_4a_Boxes(int pos_count, int q, int r) {
  int res;

  res = (dim - 3) * 2;

  switch (pos_count) {
  case 0:
    box = ((q - 1) * (4 * res)) - (4 * res) - 2 + ((r - 1) * 2);
    break;
  case 1:
    box = ((q - 1) * (4 * res)) - (4 * res) - 1 + ((r - 1) * 2);
    break;
  case 2:
    box = ((q - 1) * (4 * res)) - (2 * res) - 2 + ((r - 1) * 2);
    break;
  case 3:
    box = ((q - 1) * (4 * res)) - (2 * res) - 1 + ((r - 1) * 2);
    break;

  default:
    cout << "No box in d_LowLeft_4a_Boxes " << endl;
  }
}

void Pressure_Lattice::d_LowRight_4b_Boxes(int pos_count, int q, int r) {
  int res;

  res = (dim - 3) * 2;
  switch (pos_count) {
  case 0:
    box = ((q - 1) * (4 * res)) - (4 * res) + ((r - 1) * 2);
    break;
  case 1:
    box = ((q - 1) * (4 * res)) - (4 * res) + 1 + ((r - 1) * 2);
    break;
  case 2:
    box = ((q - 1) * (4 * res)) - (2 * res) + ((r - 1) * 2);
    break;
  case 3:
    box = ((q - 1) * (4 * res)) - (2 * res) + 1 + ((r - 1) * 2);
    break;

  default:
    cout << "No box in d_LowRight_4b_Boxes " << endl;
  }
}

void Pressure_Lattice::d_UppLeft_4c_Boxes(int pos_count, int q, int r) {
  int res;

  res = (dim - 3) * 2;

  switch (pos_count) {
  case 0:
    box = ((q - 1) * (4 * res)) - 2 + ((r - 1) * 2);
    break;
  case 1:
    box = ((q - 1) * (4 * res)) - 1 + ((r - 1) * 2);
    break;
  case 2:
    box = ((q - 1) * (4 * res)) + (2 * res) - 2 + ((r - 1) * 2);
    break;
  case 3:
    box = ((q - 1) * (4 * res)) + (2 * res) - 1 + ((r - 1) * 2);
    break;

  default:
    cout << "No box in d_UppLeft_4c_Boxes " << endl;
  }
}

void Pressure_Lattice::d_UppRight_4d_Boxes(int pos_count, int q, int r) {
  int res;

  res = (dim - 3) * 2;

  switch (pos_count) {
  case 0:
    box = ((q - 1) * (4 * res)) + ((r - 1) * 2);
    break;
  case 1:
    box = ((q - 1) * (4 * res)) + 1 + ((r - 1) * 2);
    break;
  case 2:
    box = ((q - 1) * (4 * res)) + (2 * res) + ((r - 1) * 2);
    break;
  case 3:
    box = ((q - 1) * (4 * res)) + (2 * res) + 1 + ((r - 1) * 2);
    break;

  default:
    cout << "No box in d_UppRight_4d_Boxes " << endl;
  }
}

void Pressure_Lattice::d_a_Boxes(Particle **list, int n, int o) {
  int res;

  res = (dim - 3) * 2;

  if (list[box]) {
    help = list[box];

    x_diff = abs(help->xpos - Pressure_x[n - 1][o - 1]);
    y_diff = abs(help->ypos - Pressure_y[n - 1][o - 1]);
    smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n - 1][o]);
    y_diff = abs(help->ypos - Pressure_y[n - 1][o]);
    smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n][o - 1]);
    y_diff = abs(help->ypos - Pressure_y[n][o - 1]);
    smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n][o]);
    y_diff = abs(help->ypos - Pressure_y[n][o]);
    smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    delta_P_x_1st =
        -1 * (pressure_dist[n - 1][o] -
              pressure_dist[n - 1][o - 2]); // removed divided by scale constant
    delta_P_x_2nd =
        -1 * (pressure_dist[n - 1][o + 1] - pressure_dist[n - 1][o - 1]);
    delta_P_x_3rd = -1 * (pressure_dist[n][o] - pressure_dist[n][o - 2]);
    delta_P_x_4th = -1 * (pressure_dist[n][o + 1] - pressure_dist[n][o - 1]);

    delta_P_y_1st =
        -1 * (pressure_dist[n][o - 1] - pressure_dist[n - 2][o - 1]);
    delta_P_y_2nd = -1 * (pressure_dist[n][o] - pressure_dist[n - 2][o]);
    delta_P_y_3rd =
        -1 * (pressure_dist[n + 1][o - 1] - pressure_dist[n - 1][o - 1]);
    delta_P_y_4th = -1 * (pressure_dist[n + 1][o] - pressure_dist[n - 1][o]);

    // F_P_x_1st = (delta_P_x_1st * smooth_1st) / rho[n-1][o-1];
    // F_P_x_2nd = (delta_P_x_2nd * smooth_2nd) / rho[n-1][o];
    // F_P_x_3rd = (delta_P_x_3rd * smooth_3rd) / rho[n][o-1];
    // F_P_x_4th = (delta_P_x_4th * smooth_4th) / rho[n][o];

    // F_P_y_1st = (delta_P_y_1st * smooth_1st) / rho[n-1][o-1];
    // F_P_y_2nd = (delta_P_y_2nd * smooth_2nd) / rho[n-1][o];
    // F_P_y_3rd = (delta_P_y_3rd * smooth_3rd) / rho[n][o-1];
    // F_P_y_4th = (delta_P_y_4th * smooth_4th) / rho[n][o];

    F_P_x_1st = (delta_P_x_1st * smooth_1st);
    F_P_x_2nd = (delta_P_x_2nd * smooth_2nd);
    F_P_x_3rd = (delta_P_x_3rd * smooth_3rd);
    F_P_x_4th = (delta_P_x_4th * smooth_4th);

    F_P_y_1st = (delta_P_y_1st * smooth_1st);
    F_P_y_2nd = (delta_P_y_2nd * smooth_2nd);
    F_P_y_3rd = (delta_P_y_3rd * smooth_3rd);
    F_P_y_4th = (delta_P_y_4th * smooth_4th);

    help->F_P_x = F_P_x_1st + F_P_x_2nd + F_P_x_3rd + F_P_x_4th;
    help->F_P_y = F_P_y_1st + F_P_y_2nd + F_P_y_3rd + F_P_y_4th;

    help->F_P_y = help->F_P_y * 3.14 * help->real_radius * help->real_radius /
                  rho_s[n][o];
    help->F_P_x = help->F_P_x * 3.14 * help->real_radius * help->real_radius /
                  rho_s[n][o];

    if (box == ((i - 1) * (4 * res)) - (2 * res) - 1 + ((j - 1) * 2)) {
      help->temperature = pressure_dist[n][o];
      help->fluid_pressure_gradient =
          sqrt((help->F_P_x * help->F_P_x) + (help->F_P_y * help->F_P_y));
      help->average_porosity = porosity[n][o];
      help->average_permeability = kappa[n][o];
    }

    v_darcyx_1st =
        (help->average_permeability / help->average_porosity) * delta_P_x_1st;
    v_darcyx_2nd =
        (help->average_permeability / help->average_porosity) * delta_P_x_2nd;
    v_darcyx_3rd =
        (help->average_permeability / help->average_porosity) * delta_P_x_3rd;
    v_darcyx_4th =
        (help->average_permeability / help->average_porosity) * delta_P_x_4th;
    v_darcyy_1st =
        (help->average_permeability / help->average_porosity) * delta_P_y_1st;
    v_darcyy_2nd =
        (help->average_permeability / help->average_porosity) * delta_P_y_2nd;
    v_darcyy_3rd =
        (help->average_permeability / help->average_porosity) * delta_P_y_3rd;
    v_darcyy_4th =
        (help->average_permeability / help->average_porosity) * delta_P_y_4th;

    // computed at the particule by taking average and looking velocity of the
    // fluid 11/12/13
    help->v_darcyx = (v_darcyx_1st * smooth_1st + v_darcyx_2nd * smooth_2nd +
                      v_darcyx_3rd * smooth_3rd + v_darcyx_4th * smooth_4th) /
                     (4 * mu);
    help->v_darcyy = (v_darcyy_1st * smooth_1st + v_darcyy_2nd * smooth_2nd +
                      v_darcyy_3rd * smooth_3rd + v_darcyy_4th * smooth_4th) /
                     (4 * mu);

    // help->v_darcyx = help->velx - (help->v_darcyx/mu);
    // help->v_darcyy = help->vely - (help->v_darcyy/mu);

    // cout << help->v_darcyx << endl;
    // cout << help->v_darcyy << endl;

    while (help->next_inBox) {
      x_diff = abs(help->next_inBox->xpos - Pressure_x[n - 1][o - 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n - 1][o - 1]);
      smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n - 1][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n - 1][o]);
      smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o - 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o - 1]);
      smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o]);
      smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      delta_P_x_1st =
          -1 * (pressure_dist[n - 1][o] - pressure_dist[n - 1][o - 2]);
      delta_P_x_2nd =
          -1 * (pressure_dist[n - 1][o + 1] - pressure_dist[n - 1][o - 1]);
      delta_P_x_3rd = -1 * (pressure_dist[n][o] - pressure_dist[n][o - 2]);
      delta_P_x_4th = -1 * (pressure_dist[n][o + 1] - pressure_dist[n][o - 1]);

      delta_P_y_1st =
          -1 * (pressure_dist[n][o - 1] - pressure_dist[n - 2][o - 1]);
      delta_P_y_2nd = -1 * (pressure_dist[n][o] - pressure_dist[n - 2][o]);
      delta_P_y_3rd =
          -1 * (pressure_dist[n + 1][o - 1] - pressure_dist[n - 1][o - 1]);
      delta_P_y_4th = -1 * (pressure_dist[n + 1][o] - pressure_dist[n - 1][o]);

      // F_P_x_1st = (delta_P_x_1st * smooth_1st) / rho[n-1][o-1];
      // F_P_x_2nd = (delta_P_x_2nd * smooth_2nd) / rho[n-1][o];
      // F_P_x_3rd = (delta_P_x_3rd * smooth_3rd) / rho[n][o-1];
      // F_P_x_4th = (delta_P_x_4th * smooth_4th) / rho[n][o];

      // F_P_y_1st = (delta_P_y_1st * smooth_1st) / rho[n-1][o-1];
      // F_P_y_2nd = (delta_P_y_2nd * smooth_2nd) / rho[n-1][o];
      // F_P_y_3rd = (delta_P_y_3rd * smooth_3rd) / rho[n][o-1];
      // F_P_y_4th = (delta_P_y_4th * smooth_4th) / rho[n][o];

      F_P_x_1st = (delta_P_x_1st * smooth_1st);
      F_P_x_2nd = (delta_P_x_2nd * smooth_2nd);
      F_P_x_3rd = (delta_P_x_3rd * smooth_3rd);
      F_P_x_4th = (delta_P_x_4th * smooth_4th);

      F_P_y_1st = (delta_P_y_1st * smooth_1st);
      F_P_y_2nd = (delta_P_y_2nd * smooth_2nd);
      F_P_y_3rd = (delta_P_y_3rd * smooth_3rd);
      F_P_y_4th = (delta_P_y_4th * smooth_4th);

      help->next_inBox->F_P_x = F_P_x_1st + F_P_x_2nd + F_P_x_3rd + F_P_x_4th;
      help->next_inBox->F_P_y = F_P_y_1st + F_P_y_2nd + F_P_y_3rd + F_P_y_4th;

      help->next_inBox->F_P_y = help->next_inBox->F_P_y * 3.14 *
                                help->next_inBox->real_radius *
                                help->next_inBox->real_radius / rho_s[n][o];
      help->next_inBox->F_P_x = help->next_inBox->F_P_x * 3.14 *
                                help->next_inBox->real_radius *
                                help->next_inBox->real_radius / rho_s[n][o];

      if (box == ((i - 1) * (4 * res)) - (2 * res) - 1 + ((j - 1) * 2)) {
        help->next_inBox->temperature = pressure_dist[n][o];
        help->next_inBox->fluid_pressure_gradient =
            sqrt((help->next_inBox->F_P_x * help->next_inBox->F_P_x) +
                 (help->next_inBox->F_P_y * help->next_inBox->F_P_y));
        help->next_inBox->average_porosity = porosity[n][o];
        help->next_inBox->average_permeability = kappa[n][o];
      }

      // added the 10/12/13 darcy velocity of the fluid
      // computed at each node of the box karine
      v_darcyx_1st = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_1st;
      v_darcyx_2nd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_2nd;
      v_darcyx_3rd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_3rd;
      v_darcyx_4th = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_4th;
      v_darcyy_1st = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_1st;
      v_darcyy_2nd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_2nd;
      v_darcyy_3rd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_3rd;
      v_darcyy_4th = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_4th;

      // computed at the particule by taking average and looking velocity of the
      // fluid 11/12/13 karine
      help->next_inBox->v_darcyx =
          (v_darcyx_1st * smooth_1st + v_darcyx_2nd * smooth_2nd +
           v_darcyx_3rd * smooth_3rd + v_darcyx_4th * smooth_4th) /
          (4 * mu);
      help->next_inBox->v_darcyy =
          (v_darcyy_1st * smooth_1st + v_darcyy_2nd * smooth_2nd +
           v_darcyy_3rd * smooth_3rd + v_darcyy_4th * smooth_4th) /
          (4 * mu);

      // help->next_inBox->v_darcyx = help->next_inBox->velx -
      // (help->next_inBox->v_darcyx/mu); help->next_inBox->v_darcyy =
      // help->next_inBox->vely - (help->next_inBox->v_darcyy/mu);

      help = help->next_inBox;
    }
  }
}

void Pressure_Lattice::d_b_Boxes(Particle **list, int n, int o) {

  int res;

  res = (dim - 3) * 2;

  if (list[box]) {
    help = list[box];

    x_diff = abs(help->xpos - Pressure_x[n - 1][o]);
    y_diff = abs(help->ypos - Pressure_y[n - 1][o]);
    smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n - 1][o + 1]);
    y_diff = abs(help->ypos - Pressure_y[n - 1][o + 1]);
    smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n][o]);
    y_diff = abs(help->ypos - Pressure_y[n][o]);
    smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n][o + 1]);
    y_diff = abs(help->ypos - Pressure_y[n][o + 1]);
    smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    delta_P_x_1st =
        -1 * (pressure_dist[n - 1][o + 1] - pressure_dist[n - 1][o - 1]);
    delta_P_x_2nd =
        -1 * (pressure_dist[n - 1][o + 2] - pressure_dist[n - 1][o]);
    delta_P_x_3rd = -1 * (pressure_dist[n][o + 1] - pressure_dist[n][o - 1]);
    delta_P_x_4th = -1 * (pressure_dist[n][o + 2] - pressure_dist[n][o]);

    delta_P_y_1st = -1 * (pressure_dist[n][o] - pressure_dist[n - 2][o]);
    delta_P_y_2nd =
        -1 * (pressure_dist[n][o + 1] - pressure_dist[n - 2][o + 1]);
    delta_P_y_3rd = -1 * (pressure_dist[n + 1][o] - pressure_dist[n - 1][o]);
    delta_P_y_4th =
        -1 * (pressure_dist[n + 1][o + 1] - pressure_dist[n - 1][o + 1]);

    // F_P_x_1st = (delta_P_x_1st * smooth_1st) / rho[n-1][o];
    // F_P_x_2nd = (delta_P_x_2nd * smooth_2nd) / rho[n-1][o+1];
    // F_P_x_3rd = (delta_P_x_3rd * smooth_3rd) / rho[n][o];
    // F_P_x_4th = (delta_P_x_4th * smooth_4th) / rho[n][o+1];

    // F_P_y_1st = (delta_P_y_1st * smooth_1st) / rho[n-1][o];
    // F_P_y_2nd = (delta_P_y_2nd * smooth_2nd) / rho[n-1][o+1];
    // F_P_y_3rd = (delta_P_y_3rd * smooth_3rd) / rho[n][o];
    // F_P_y_4th = (delta_P_y_4th * smooth_4th) / rho[n][o+1];

    F_P_x_1st = (delta_P_x_1st * smooth_1st);
    F_P_x_2nd = (delta_P_x_2nd * smooth_2nd);
    F_P_x_3rd = (delta_P_x_3rd * smooth_3rd);
    F_P_x_4th = (delta_P_x_4th * smooth_4th);

    F_P_y_1st = (delta_P_y_1st * smooth_1st);
    F_P_y_2nd = (delta_P_y_2nd * smooth_2nd);
    F_P_y_3rd = (delta_P_y_3rd * smooth_3rd);
    F_P_y_4th = (delta_P_y_4th * smooth_4th);

    help->F_P_x = F_P_x_1st + F_P_x_2nd + F_P_x_3rd + F_P_x_4th;
    help->F_P_y = F_P_y_1st + F_P_y_2nd + F_P_y_3rd + F_P_y_4th;

    help->F_P_y = help->F_P_y * 3.14 * help->real_radius * help->real_radius /
                  rho_s[n][o];
    help->F_P_x = help->F_P_x * 3.14 * help->real_radius * help->real_radius /
                  rho_s[n][o];

    if (box == ((i - 1) * (4 * res)) - (2 * res) + ((j - 1) * 2)) {
      help->temperature = pressure_dist[n][o];
      help->fluid_pressure_gradient =
          sqrt((help->F_P_x * help->F_P_x) + (help->F_P_y * help->F_P_y));
      help->average_porosity = porosity[n][o];
      help->average_permeability = kappa[n][o];
    }

    v_darcyx_1st =
        (help->average_permeability / help->average_porosity) * delta_P_x_1st;
    v_darcyx_2nd =
        (help->average_permeability / help->average_porosity) * delta_P_x_2nd;
    v_darcyx_3rd =
        (help->average_permeability / help->average_porosity) * delta_P_x_3rd;
    v_darcyx_4th =
        (help->average_permeability / help->average_porosity) * delta_P_x_4th;
    v_darcyy_1st =
        (help->average_permeability / help->average_porosity) * delta_P_y_1st;
    v_darcyy_2nd =
        (help->average_permeability / help->average_porosity) * delta_P_y_2nd;
    v_darcyy_3rd =
        (help->average_permeability / help->average_porosity) * delta_P_y_3rd;
    v_darcyy_4th =
        (help->average_permeability / help->average_porosity) * delta_P_y_4th;

    // computed at the particule by taking average and looking velocity of the
    // fluid 11/12/13
    help->v_darcyx = (v_darcyx_1st * smooth_1st + v_darcyx_2nd * smooth_2nd +
                      v_darcyx_3rd * smooth_3rd + v_darcyx_4th * smooth_4th) /
                     (4 * mu);
    help->v_darcyy = (v_darcyy_1st * smooth_1st + v_darcyy_2nd * smooth_2nd +
                      v_darcyy_3rd * smooth_3rd + v_darcyy_4th * smooth_4th) /
                     (4 * mu);

    // help->v_darcyx = help->velx - (help->v_darcyx/mu);
    // help->v_darcyy = help->vely - (help->v_darcyy/mu);

    while (help->next_inBox) {
      x_diff = abs(help->next_inBox->xpos - Pressure_x[n - 1][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n - 1][o]);
      smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n - 1][o + 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n - 1][o + 1]);
      smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o]);
      smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o + 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o + 1]);
      smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      delta_P_x_1st =
          -1 * (pressure_dist[n - 1][o + 1] - pressure_dist[n - 1][o - 1]);
      delta_P_x_2nd =
          -1 * (pressure_dist[n - 1][o + 2] - pressure_dist[n - 1][o]);
      delta_P_x_3rd = -1 * (pressure_dist[n][o + 1] - pressure_dist[n][o - 1]);
      delta_P_x_4th = -1 * (pressure_dist[n][o + 2] - pressure_dist[n][o]);

      delta_P_y_1st = -1 * (pressure_dist[n][o] - pressure_dist[n - 2][o]);
      delta_P_y_2nd =
          -1 * (pressure_dist[n][o + 1] - pressure_dist[n - 2][o + 1]);
      delta_P_y_3rd = -1 * (pressure_dist[n + 1][o] - pressure_dist[n - 1][o]);
      delta_P_y_4th =
          -1 * (pressure_dist[n + 1][o + 1] - pressure_dist[n - 1][o + 1]);

      // F_P_x_1st = (delta_P_x_1st * smooth_1st) / rho[n-1][o];
      // F_P_x_2nd = (delta_P_x_2nd * smooth_2nd) / rho[n-1][o+1];
      // F_P_x_3rd = (delta_P_x_3rd * smooth_3rd) / rho[n][o];
      // F_P_x_4th = (delta_P_x_4th * smooth_4th) / rho[n][o+1];

      // F_P_y_1st = (delta_P_y_1st * smooth_1st) / rho[n-1][o];
      // F_P_y_2nd = (delta_P_y_2nd * smooth_2nd) / rho[n-1][o+1];
      // F_P_y_3rd = (delta_P_y_3rd * smooth_3rd) / rho[n][o];
      // F_P_y_4th = (delta_P_y_4th * smooth_4th) / rho[n][o+1];

      F_P_x_1st = (delta_P_x_1st * smooth_1st);
      F_P_x_2nd = (delta_P_x_2nd * smooth_2nd);
      F_P_x_3rd = (delta_P_x_3rd * smooth_3rd);
      F_P_x_4th = (delta_P_x_4th * smooth_4th);

      F_P_y_1st = (delta_P_y_1st * smooth_1st);
      F_P_y_2nd = (delta_P_y_2nd * smooth_2nd);
      F_P_y_3rd = (delta_P_y_3rd * smooth_3rd);
      F_P_y_4th = (delta_P_y_4th * smooth_4th);

      help->next_inBox->F_P_x = F_P_x_1st + F_P_x_2nd + F_P_x_3rd + F_P_x_4th;
      help->next_inBox->F_P_y = F_P_y_1st + F_P_y_2nd + F_P_y_3rd + F_P_y_4th;

      help->next_inBox->F_P_y = help->next_inBox->F_P_y * 3.14 *
                                help->next_inBox->real_radius *
                                help->next_inBox->real_radius / rho_s[n][o];
      help->next_inBox->F_P_x = help->next_inBox->F_P_x * 3.14 *
                                help->next_inBox->real_radius *
                                help->next_inBox->real_radius / rho_s[n][o];

      if (box == ((i - 1) * (4 * res)) - (2 * res) + ((j - 1) * 2)) {
        help->next_inBox->temperature = pressure_dist[n][o];
        help->next_inBox->fluid_pressure_gradient =
            sqrt((help->next_inBox->F_P_x * help->next_inBox->F_P_x) +
                 (help->next_inBox->F_P_y * help->next_inBox->F_P_y));
        help->next_inBox->average_porosity = porosity[n][o];
        help->next_inBox->average_permeability = kappa[n][o];
      }

      // added the 10/12/13 darcy velocity of the fluid
      // computed at each node of the box karine
      v_darcyx_1st = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_1st;
      v_darcyx_2nd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_2nd;
      v_darcyx_3rd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_3rd;
      v_darcyx_4th = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_4th;
      v_darcyy_1st = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_1st;
      v_darcyy_2nd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_2nd;
      v_darcyy_3rd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_3rd;
      v_darcyy_4th = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_4th;

      // computed at the particule by taking average and looking velocity of the
      // fluid 11/12/13 karine
      help->next_inBox->v_darcyx =
          (v_darcyx_1st * smooth_1st + v_darcyx_2nd * smooth_2nd +
           v_darcyx_3rd * smooth_3rd + v_darcyx_4th * smooth_4th) /
          (4 * mu);
      help->next_inBox->v_darcyy =
          (v_darcyy_1st * smooth_1st + v_darcyy_2nd * smooth_2nd +
           v_darcyy_3rd * smooth_3rd + v_darcyy_4th * smooth_4th) /
          (4 * mu);

      // help->next_inBox->v_darcyx = help->next_inBox->velx -
      // (help->next_inBox->v_darcyx/mu); help->next_inBox->v_darcyy =
      // help->next_inBox->vely - (help->next_inBox->v_darcyy/mu);

      help = help->next_inBox;
    }
  }
}

void Pressure_Lattice::d_c_Boxes(Particle **list, int n, int o) {
  int res;

  res = (dim - 3) * 2;

  if (list[box]) {
    help = list[box];

    x_diff = abs(help->xpos - Pressure_x[n][o - 1]);
    y_diff = abs(help->ypos - Pressure_y[n][o - 1]);
    smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n][o]);
    y_diff = abs(help->ypos - Pressure_y[n][o]);
    smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n + 1][o - 1]);
    y_diff = abs(help->ypos - Pressure_y[n + 1][o - 1]);
    smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n + 1][o]);
    y_diff = abs(help->ypos - Pressure_y[n + 1][o]);
    smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    delta_P_x_1st = -1 * (pressure_dist[n][o] - pressure_dist[n][o - 2]);
    delta_P_x_2nd = -1 * (pressure_dist[n][o + 1] - pressure_dist[n][o - 1]);
    delta_P_x_3rd =
        -1 * (pressure_dist[n + 1][o] - pressure_dist[n + 1][o - 2]);
    delta_P_x_4th =
        -1 * (pressure_dist[n + 1][o + 1] - pressure_dist[n + 1][o - 1]);

    delta_P_y_1st =
        -1 * (pressure_dist[n + 1][o - 1] - pressure_dist[n - 1][o - 1]);
    delta_P_y_2nd = -1 * (pressure_dist[n + 1][o] - pressure_dist[n - 1][o]);
    delta_P_y_3rd =
        -1 * (pressure_dist[n + 2][o - 1] - pressure_dist[n][o - 1]);
    delta_P_y_4th = -1 * (pressure_dist[n + 2][o] - pressure_dist[n][o]);

    // F_P_x_1st = (delta_P_x_1st * smooth_1st) / rho[n][o-1];
    // F_P_x_2nd = (delta_P_x_2nd * smooth_2nd) / rho[n][o];
    // F_P_x_3rd = (delta_P_x_3rd * smooth_3rd) / rho[n+1][o-1];
    // F_P_x_4th = (delta_P_x_4th * smooth_4th) / rho[n+1][o];

    // F_P_y_1st = (delta_P_y_1st * smooth_1st) / rho[n][o-1];
    // F_P_y_2nd = (delta_P_y_2nd * smooth_2nd) / rho[n][o];
    // F_P_y_3rd = (delta_P_y_3rd * smooth_3rd) / rho[n+1][o-1];
    // F_P_y_4th = (delta_P_y_4th * smooth_4th) / rho[n+1][o];

    F_P_x_1st = (delta_P_x_1st * smooth_1st);
    F_P_x_2nd = (delta_P_x_2nd * smooth_2nd);
    F_P_x_3rd = (delta_P_x_3rd * smooth_3rd);
    F_P_x_4th = (delta_P_x_4th * smooth_4th);

    F_P_y_1st = (delta_P_y_1st * smooth_1st);
    F_P_y_2nd = (delta_P_y_2nd * smooth_2nd);
    F_P_y_3rd = (delta_P_y_3rd * smooth_3rd);
    F_P_y_4th = (delta_P_y_4th * smooth_4th);

    help->F_P_x = F_P_x_1st + F_P_x_2nd + F_P_x_3rd + F_P_x_4th;
    help->F_P_y = F_P_y_1st + F_P_y_2nd + F_P_y_3rd + F_P_y_4th;

    help->F_P_y = help->F_P_y * 3.14 * help->real_radius * help->real_radius /
                  rho_s[n][o];
    help->F_P_x = help->F_P_x * 3.14 * help->real_radius * help->real_radius /
                  rho_s[n][o];

    if (box == ((i - 1) * (4 * res)) - 1 + ((j - 1) * 2)) {
      help->temperature = pressure_dist[n][o];
      help->fluid_pressure_gradient =
          sqrt((help->F_P_x * help->F_P_x) + (help->F_P_y * help->F_P_y));
      help->average_porosity = porosity[n][o];
      help->average_permeability = kappa[n][o];
    }

    v_darcyx_1st =
        (help->average_permeability / help->average_porosity) * delta_P_x_1st;
    v_darcyx_2nd =
        (help->average_permeability / help->average_porosity) * delta_P_x_2nd;
    v_darcyx_3rd =
        (help->average_permeability / help->average_porosity) * delta_P_x_3rd;
    v_darcyx_4th =
        (help->average_permeability / help->average_porosity) * delta_P_x_4th;
    v_darcyy_1st =
        (help->average_permeability / help->average_porosity) * delta_P_y_1st;
    v_darcyy_2nd =
        (help->average_permeability / help->average_porosity) * delta_P_y_2nd;
    v_darcyy_3rd =
        (help->average_permeability / help->average_porosity) * delta_P_y_3rd;
    v_darcyy_4th =
        (help->average_permeability / help->average_porosity) * delta_P_y_4th;

    // computed at the particule by taking average and looking velocity of the
    // fluid 11/12/13
    help->v_darcyx = (v_darcyx_1st * smooth_1st + v_darcyx_2nd * smooth_2nd +
                      v_darcyx_3rd * smooth_3rd + v_darcyx_4th * smooth_4th) /
                     (4 * mu);
    help->v_darcyy = (v_darcyy_1st * smooth_1st + v_darcyy_2nd * smooth_2nd +
                      v_darcyy_3rd * smooth_3rd + v_darcyy_4th * smooth_4th) /
                     (4 * mu);

    // help->v_darcyx = help->velx - (help->v_darcyx/mu);
    // help->v_darcyy = help->vely - (help->v_darcyy/mu);

    while (help->next_inBox) {
      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o - 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o - 1]);
      smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o]);
      smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n + 1][o - 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n + 1][o - 1]);
      smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n + 1][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n + 1][o]);
      smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      delta_P_x_1st = -1 * (pressure_dist[n][o] - pressure_dist[n][o - 2]);
      delta_P_x_2nd = -1 * (pressure_dist[n][o + 1] - pressure_dist[n][o - 1]);
      delta_P_x_3rd =
          -1 * (pressure_dist[n + 1][o] - pressure_dist[n + 1][o - 2]);
      delta_P_x_4th =
          -1 * (pressure_dist[n + 1][o + 1] - pressure_dist[n + 1][o - 1]);

      delta_P_y_1st =
          -1 * (pressure_dist[n + 1][o - 1] - pressure_dist[n - 1][o - 1]);
      delta_P_y_2nd = -1 * (pressure_dist[n + 1][o] - pressure_dist[n - 1][o]);
      delta_P_y_3rd =
          -1 * (pressure_dist[n + 2][o - 1] - pressure_dist[n][o - 1]);
      delta_P_y_4th = -1 * (pressure_dist[n + 2][o] - pressure_dist[n][o]);

      // F_P_x_1st = (delta_P_x_1st * smooth_1st) / rho[n][o-1];
      // F_P_x_2nd = (delta_P_x_2nd * smooth_2nd) / rho[n][o];
      // F_P_x_3rd = (delta_P_x_3rd * smooth_3rd) / rho[n+1][o-1];
      // F_P_x_4th = (delta_P_x_4th * smooth_4th) / rho[n+1][o];

      // F_P_y_1st = (delta_P_y_1st * smooth_1st) / rho[n][o-1];
      // F_P_y_2nd = (delta_P_y_2nd * smooth_2nd) / rho[n][o];
      // F_P_y_3rd = (delta_P_y_3rd * smooth_3rd) / rho[n+1][o-1];
      // F_P_y_4th = (delta_P_y_4th * smooth_4th) / rho[n+1][o];

      F_P_x_1st = (delta_P_x_1st * smooth_1st);
      F_P_x_2nd = (delta_P_x_2nd * smooth_2nd);
      F_P_x_3rd = (delta_P_x_3rd * smooth_3rd);
      F_P_x_4th = (delta_P_x_4th * smooth_4th);

      F_P_y_1st = (delta_P_y_1st * smooth_1st);
      F_P_y_2nd = (delta_P_y_2nd * smooth_2nd);
      F_P_y_3rd = (delta_P_y_3rd * smooth_3rd);
      F_P_y_4th = (delta_P_y_4th * smooth_4th);

      help->next_inBox->F_P_x = F_P_x_1st + F_P_x_2nd + F_P_x_3rd + F_P_x_4th;
      help->next_inBox->F_P_y = F_P_y_1st + F_P_y_2nd + F_P_y_3rd + F_P_y_4th;

      help->next_inBox->F_P_y = help->next_inBox->F_P_y * 3.14 *
                                help->next_inBox->real_radius *
                                help->next_inBox->real_radius / rho_s[n][o];
      help->next_inBox->F_P_x = help->next_inBox->F_P_x * 3.14 *
                                help->next_inBox->real_radius *
                                help->next_inBox->real_radius / rho_s[n][o];

      if (box == ((i - 1) * (4 * res)) - 1 + ((j - 1) * 2)) {
        help->next_inBox->temperature = pressure_dist[n][o];
        help->next_inBox->fluid_pressure_gradient =
            sqrt((help->next_inBox->F_P_x * help->next_inBox->F_P_x) +
                 (help->next_inBox->F_P_y * help->next_inBox->F_P_y));
        help->next_inBox->average_porosity = porosity[n][o];
        help->next_inBox->average_permeability = kappa[n][o];
      }

      // added the 10/12/13 darcy velocity of the fluid
      // computed at each node of the box karine
      v_darcyx_1st = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_1st;
      v_darcyx_2nd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_2nd;
      v_darcyx_3rd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_3rd;
      v_darcyx_4th = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_4th;
      v_darcyy_1st = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_1st;
      v_darcyy_2nd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_2nd;
      v_darcyy_3rd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_3rd;
      v_darcyy_4th = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_4th;

      // computed at the particule by taking average and looking velocity of the
      // fluid 11/12/13 karine
      help->next_inBox->v_darcyx =
          (v_darcyx_1st * smooth_1st + v_darcyx_2nd * smooth_2nd +
           v_darcyx_3rd * smooth_3rd + v_darcyx_4th * smooth_4th) /
          (4 * mu);
      help->next_inBox->v_darcyy =
          (v_darcyy_1st * smooth_1st + v_darcyy_2nd * smooth_2nd +
           v_darcyy_3rd * smooth_3rd + v_darcyy_4th * smooth_4th) /
          (4 * mu);

      // help->next_inBox->v_darcyx = help->next_inBox->velx -
      // (help->next_inBox->v_darcyx/mu); help->next_inBox->v_darcyy =
      // help->next_inBox->vely - (help->next_inBox->v_darcyy/mu);

      help = help->next_inBox;
    }
  }
}

void Pressure_Lattice::d_d_Boxes(Particle **list, int n, int o) {
  //~ cout << "p " << p << " "<< " k " << k << endl;

  int res;

  res = (dim - 3) * 2;

  if (list[box]) {
    help = list[box];
    //~ cout << box << " " << "xpos " << help->xpos << " " ;

    x_diff = abs(help->xpos - Pressure_x[n][o]);
    y_diff = abs(help->ypos - Pressure_y[n][o]);
    smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n][o + 1]);
    y_diff = abs(help->ypos - Pressure_y[n][o + 1]);
    smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n + 1][o]);
    y_diff = abs(help->ypos - Pressure_y[n + 1][o]);
    smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n + 1][o + 1]);
    y_diff = abs(help->ypos - Pressure_y[n + 1][o + 1]);
    smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    delta_P_x_1st = -1 * (pressure_dist[n][o + 1] - pressure_dist[n][o - 1]);
    delta_P_x_2nd = -1 * (pressure_dist[n][o + 2] - pressure_dist[n][o]);
    delta_P_x_3rd =
        -1 * (pressure_dist[n + 1][o + 1] - pressure_dist[n + 1][o - 1]);
    delta_P_x_4th =
        -1 * (pressure_dist[n + 1][o + 2] - pressure_dist[n + 1][o]);

    delta_P_y_1st = -1 * (pressure_dist[n + 1][o] - pressure_dist[n - 1][o]);
    delta_P_y_2nd =
        -1 * (pressure_dist[n + 1][o + 1] - pressure_dist[n - 1][o + 1]);
    delta_P_y_3rd = -1 * (pressure_dist[n + 2][o] - pressure_dist[n][o]);
    delta_P_y_4th =
        -1 * (pressure_dist[n + 2][o + 1] - pressure_dist[n][o + 1]);

    // F_P_x_1st = (delta_P_x_1st * smooth_1st) / rho[n][o];
    // F_P_x_2nd = (delta_P_x_2nd * smooth_2nd) / rho[n][o+1];
    // F_P_x_3rd = (delta_P_x_3rd * smooth_3rd) / rho[n+1][o];
    // F_P_x_4th = (delta_P_x_4th * smooth_4th) / rho[n+1][o+1];

    // F_P_y_1st = (delta_P_y_1st * smooth_1st) / rho[n][o];
    // F_P_y_2nd = (delta_P_y_2nd * smooth_2nd) / rho[n][o+1];
    // F_P_y_3rd = (delta_P_y_3rd * smooth_3rd) / rho[n+1][o];
    // F_P_y_4th = (delta_P_y_4th * smooth_4th) / rho[n+1][o+1];

    F_P_x_1st = (delta_P_x_1st * smooth_1st);
    F_P_x_2nd = (delta_P_x_2nd * smooth_2nd);
    F_P_x_3rd = (delta_P_x_3rd * smooth_3rd);
    F_P_x_4th = (delta_P_x_4th * smooth_4th);

    F_P_y_1st = (delta_P_y_1st * smooth_1st);
    F_P_y_2nd = (delta_P_y_2nd * smooth_2nd);
    F_P_y_3rd = (delta_P_y_3rd * smooth_3rd);
    F_P_y_4th = (delta_P_y_4th * smooth_4th);

    help->F_P_x = F_P_x_1st + F_P_x_2nd + F_P_x_3rd + F_P_x_4th;
    help->F_P_y = F_P_y_1st + F_P_y_2nd + F_P_y_3rd + F_P_y_4th;

    help->F_P_y = help->F_P_y * 3.14 * help->real_radius * help->real_radius /
                  rho_s[n][o];
    help->F_P_x = help->F_P_x * 3.14 * help->real_radius * help->real_radius /
                  rho_s[n][o];

    if (box == ((i - 1) * (4 * res)) + ((j - 1) * 2)) {
      help->temperature = pressure_dist[n][o];
      help->fluid_pressure_gradient =
          sqrt((help->F_P_x * help->F_P_x) + (help->F_P_y * help->F_P_y));
      help->average_porosity = porosity[n][o];
      help->average_permeability = kappa[n][o];
      // cout << "in " << endl;
    }

    v_darcyx_1st =
        (help->average_permeability / help->average_porosity) * delta_P_x_1st;
    v_darcyx_2nd =
        (help->average_permeability / help->average_porosity) * delta_P_x_2nd;
    v_darcyx_3rd =
        (help->average_permeability / help->average_porosity) * delta_P_x_3rd;
    v_darcyx_4th =
        (help->average_permeability / help->average_porosity) * delta_P_x_4th;
    v_darcyy_1st =
        (help->average_permeability / help->average_porosity) * delta_P_y_1st;
    v_darcyy_2nd =
        (help->average_permeability / help->average_porosity) * delta_P_y_2nd;
    v_darcyy_3rd =
        (help->average_permeability / help->average_porosity) * delta_P_y_3rd;
    v_darcyy_4th =
        (help->average_permeability / help->average_porosity) * delta_P_y_4th;

    // computed at the particule by taking average and looking velocity of the
    // fluid 11/12/13
    help->v_darcyx = (v_darcyx_1st * smooth_1st + v_darcyx_2nd * smooth_2nd +
                      v_darcyx_3rd * smooth_3rd + v_darcyx_4th * smooth_4th) /
                     (4 * mu);
    help->v_darcyy = (v_darcyy_1st * smooth_1st + v_darcyy_2nd * smooth_2nd +
                      v_darcyy_3rd * smooth_3rd + v_darcyy_4th * smooth_4th) /
                     (4 * mu);

    // help->v_darcyx = help->velx - (help->v_darcyx/mu);
    // help->v_darcyy = help->vely - (help->v_darcyy/mu);

    while (help->next_inBox) {
      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o]);
      smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o + 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o + 1]);
      smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n + 1][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n + 1][o]);
      smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n + 1][o + 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n + 1][o + 1]);
      smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      delta_P_x_1st = -1 * (pressure_dist[n][o + 1] - pressure_dist[n][o - 1]);
      delta_P_x_2nd = -1 * (pressure_dist[n][o + 2] - pressure_dist[n][o]);
      delta_P_x_3rd =
          -1 * (pressure_dist[n + 1][o + 1] - pressure_dist[n + 1][o - 1]);
      delta_P_x_4th =
          -1 * (pressure_dist[n + 1][o + 2] - pressure_dist[n + 1][o]);

      delta_P_y_1st = -1 * (pressure_dist[n + 1][o] - pressure_dist[n - 1][o]);
      delta_P_y_2nd =
          -1 * (pressure_dist[n + 1][o + 1] - pressure_dist[n - 1][o + 1]);
      delta_P_y_3rd = -1 * (pressure_dist[n + 2][o] - pressure_dist[n][o]);
      delta_P_y_4th =
          -1 * (pressure_dist[n + 2][o + 1] - pressure_dist[n][o + 1]);

      // F_P_x_1st = (delta_P_x_1st * smooth_1st) / rho[n][o];
      // F_P_x_2nd = (delta_P_x_2nd * smooth_2nd) / rho[n][o+1];
      // F_P_x_3rd = (delta_P_x_3rd * smooth_3rd) / rho[n+1][o];
      // F_P_x_4th = (delta_P_x_4th * smooth_4th) / rho[n+1][o+1];

      // F_P_y_1st = (delta_P_y_1st * smooth_1st) / rho[n][o];
      // F_P_y_2nd = (delta_P_y_2nd * smooth_2nd) / rho[n][o+1];
      // F_P_y_3rd = (delta_P_y_3rd * smooth_3rd) / rho[n+1][o];
      // F_P_y_4th = (delta_P_y_4th * smooth_4th) / rho[n+1][o+1];

      F_P_x_1st = (delta_P_x_1st * smooth_1st);
      F_P_x_2nd = (delta_P_x_2nd * smooth_2nd);
      F_P_x_3rd = (delta_P_x_3rd * smooth_3rd);
      F_P_x_4th = (delta_P_x_4th * smooth_4th);

      F_P_y_1st = (delta_P_y_1st * smooth_1st);
      F_P_y_2nd = (delta_P_y_2nd * smooth_2nd);
      F_P_y_3rd = (delta_P_y_3rd * smooth_3rd);
      F_P_y_4th = (delta_P_y_4th * smooth_4th);

      help->next_inBox->F_P_x = F_P_x_1st + F_P_x_2nd + F_P_x_3rd + F_P_x_4th;
      help->next_inBox->F_P_y = F_P_y_1st + F_P_y_2nd + F_P_y_3rd + F_P_y_4th;

      help->next_inBox->F_P_y = help->next_inBox->F_P_y * 3.14 *
                                help->next_inBox->real_radius *
                                help->next_inBox->real_radius / rho_s[n][o];
      help->next_inBox->F_P_x = help->next_inBox->F_P_x * 3.14 *
                                help->next_inBox->real_radius *
                                help->next_inBox->real_radius / rho_s[n][o];

      if (box == ((i - 1) * (4 * res)) + ((j - 1) * 2)) {
        help->next_inBox->temperature = pressure_dist[n][o];
        help->next_inBox->fluid_pressure_gradient =
            sqrt((help->next_inBox->F_P_x * help->next_inBox->F_P_x) +
                 (help->next_inBox->F_P_y * help->next_inBox->F_P_y));
        help->next_inBox->average_porosity = porosity[n][o];
        help->next_inBox->average_permeability = kappa[n][o];
      }

      // added the 10/12/13 darcy velocity of the fluid
      // computed at each node of the box karine
      v_darcyx_1st = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_1st;
      v_darcyx_2nd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_2nd;
      v_darcyx_3rd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_3rd;
      v_darcyx_4th = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_x_4th;
      v_darcyy_1st = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_1st;
      v_darcyy_2nd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_2nd;
      v_darcyy_3rd = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_3rd;
      v_darcyy_4th = (help->next_inBox->average_permeability /
                      help->next_inBox->average_porosity) *
                     delta_P_y_4th;

      // computed at the particule by taking average and looking velocity of the
      // fluid 11/12/13 karine
      help->next_inBox->v_darcyx =
          (v_darcyx_1st * smooth_1st + v_darcyx_2nd * smooth_2nd +
           v_darcyx_3rd * smooth_3rd + v_darcyx_4th * smooth_4th) /
          (4 * mu);
      help->next_inBox->v_darcyy =
          (v_darcyy_1st * smooth_1st + v_darcyy_2nd * smooth_2nd +
           v_darcyy_3rd * smooth_3rd + v_darcyy_4th * smooth_4th) /
          (4 * mu);

      // help->next_inBox->v_darcyx = help->next_inBox->velx -
      // (help->next_inBox->v_darcyx); help->next_inBox->v_darcyy =
      // help->next_inBox->vely - (help->next_inBox->v_darcyy);

      help = help->next_inBox;
    }
  }
}

void Pressure_Lattice::boundary_set_hydrostatic(int box_size, int depth,
                                                double increase) {

  for (j = 0; j < dim; j++) {

    pressure_dist[0][j] = ((depth + box_size) * 1000 * 9.81) + increase;

    pressure_dist[dim - 1][j] = depth * 1000 * 9.81;
  }

  for (i = 0; i < dim; i++) {
    pressure_dist[i][0] = pressure_dist[i][1];
    pressure_dist[i][dim - 1] = pressure_dist[i][dim - 2];
  }
}
void Pressure_Lattice::boundary_set_hydrostatic_gradient(int box_size,
                                                         int depth) {

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {

      pressure_dist[i][j] =
          (depth + (box_size * (float((dim - 1) - i) / (dim - 1)))) * 1000 *
          9.81;

      // pressure_dist[i][j] =  (box_size+depth)*1000*9.81;
    }
  }

  // for(i = 0; i < dim; i++)
  //{
  //	pressure_dist[i][0] = pressure_dist[i][1];
  // pressure_dist[i][dim-1] = pressure_dist[i][dim-2];

  //}
}
void Pressure_Lattice::boundary_fixed(int boundry_sourc,
                                      double source_val_boundry,
                                      int z_distance) {
  for (i = 0; i < dim; i++) {
    pressure_dist[i][0] = pressure_dist[i][1];
    pressure_dist[i][dim - 1] = pressure_dist[i][dim - 2];
  }

  for (j = 0; j < dim; j++) {
    if (boundry_sourc)
      pressure_dist[0][j] = source_val_boundry;
    else
      pressure_dist[0][j] = pressure_dist[1][j];

    pressure_dist[dim - 1][j] = pressure_dist[dim - 2][j];
  }
}

void Pressure_Lattice::rand_source_PressIncrement(int press_y, int press_x,
                                                  double pressure_val) {
  pressure_dist[press_y][press_x] += pressure_val;
}

void Pressure_Lattice::source_hor_PressIncrement(int press_y, int press_min_x,
                                                 int press_max_x,
                                                 double pressure_val) {
  for (i = press_min_x; i <= press_max_x; i++)
    pressure_dist[press_y][i] += pressure_val;
}

void Pressure_Lattice::source_ver_PressIncrement(int press_x, int press_min_y,
                                                 int press_max_y,
                                                 double pressure_val) {
  for (i = press_min_y; i <= press_max_y; i++)
    pressure_dist[i][press_x] += pressure_val;
}

void Pressure_Lattice::source_area_PressIncrement(int press_min_y,
                                                  int press_max_y,
                                                  int press_min_x,
                                                  int press_max_x,
                                                  double pressure_val) {
  for (i = press_min_y; i <= press_max_y; i++) {
    for (j = press_min_x; j <= press_max_x; j++)
      pressure_dist[i][j] += pressure_val;
  }
}

void Pressure_Lattice::test_pressure_with_source(double time) {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      source_a[i][j] = 0;

      if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
        source_a[i][j] = ((vel_x[i][j + 1] - vel_x[i][j - 1]) /
                          (2 * delta_x * scale_constt));
        source_a[i][j] += ((vel_y[i + 1][j] - vel_y[i - 1][j]) /
                           (2.0 * delta_y * scale_constt));
        source_a[i][j] *= (1.0 / compressibility + sourced_pressure[i][j]);
        source_a[i][j] *= time_a / (2.0 * porosity[i][j]);
      }
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      source_b[i][j] = 0;

      if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
        source_b[i][j] =
            vel_x[i][j] *
            ((sourced_pressure[i][j + 1] - sourced_pressure[i][j - 1]) /
             (2 * delta_x * scale_constt));
        source_b[i][j] +=
            vel_y[i][j] *
            ((sourced_pressure[i + 1][j] - sourced_pressure[i - 1][j]) /
             (2.0 * delta_y * scale_constt));
        source_b[i][j] *= time_a / 2.0;
      }
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      source[i][j] = 0;

      source[i][j] = source_a[i][j] + source_b[i][j];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++)
      interim_pressure_dist[i][j] -= source[i][j];
  }
}

void Pressure_Lattice::dan_matrices_set_1() {
  // creats L.H.S (al[][][])  and R.H.S (bl[][][])coeffiecient matrix
  // at first half time step of ADI method

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
            al[k][i][j] = 1 + 2 * alpha[k][i];
            bl[k][i][j] = 1 - 2 * beta[i][k];
          }
        } else if ((i == j - 1 && i > 0) || (i == j + 1 && i < (dim - 1))) {
          al[k][i][j] = -alpha[k][i];

          bl[k][i][j] = beta[i][k];
        }
      }
    }
  }
}

void Pressure_Lattice::dan_multiplication_1st() {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      interim_pressure_dist[i][j] = 0;

      for (k = 0; k < dim; k++) {
        interim_pressure_dist[i][j] += bl[j][i][k] * pressure_dist[k][j];
      }
    }
  }
}

void Pressure_Lattice::dan_invert() {
  double factor;
  int g;

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

void Pressure_Lattice::dan_multiplication_inv() {

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      pressure_dist[i][j] = 0;

      for (k = 0; k < dim; k++)
        pressure_dist[i][j] += a_invl[j][i][k] * interim_pressure_dist[k][j];

      // cout << pressure_dist[i][j] << " ";
      // if (j == dim-1)
      // cout << endl << j << endl;
    }
  }
}

void Pressure_Lattice::dan_matrices_set_2() {
  // creats L.H.S (al[][][])  and R.H.S (bl[][][])coeffiecient matrix
  // at second half time step of ADI method

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
            al[k][i][j] = 1 + 2 * beta[i][k];
            bl[k][i][j] = 1 - 2 * alpha[k][i];
          }
        } else if ((i == j - 1 && i > 0) || (i == j + 1 && i < (dim - 1))) {
          al[k][i][j] = -beta[i][k];

          bl[k][i][j] = alpha[k][i];
        }
      }
    }
  }
}

void Pressure_Lattice::grad_pressure(float y_limit) {
  seal_y_limit = y_limit * dim;

  for (i = seal_y_limit; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      pressure_dist[i][j] = 0;
    }
  }
}

/**********************************
 * functions for the Advection Diffusion
 *
 * *******************************/

void Pressure_Lattice::setconcent(double concentration) {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      concent[i][j] = concentration;
    }
  }
}

void Pressure_Lattice::concentration_particles(Particle **list) {
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
        if (i == 1 && j == 1) {
          conc_1x4_Boxes(list, 1, i, j);

        }

        else if (i == 1 && (j > 1 && j < dim - 2)) {
          conc_2x4_Boxes(list, 1, i, j);

        }

        else if (i == 1 && j == dim - 2) {
          conc_1x4_Boxes(list, 2, i, j);

        }

        else if ((i > 1 && i < dim - 2) && j == 1) {
          conc_2x4_Boxes(list, 2, i, j);

        }

        else if ((i > 1 && i < dim - 2) && j == dim - 2) {
          conc_2x4_Boxes(list, 3, i, j);

        }

        else if (i == dim - 2 && j == 1) {
          conc_1x4_Boxes(list, 3, i, j);

        }

        else if (i == dim - 2 && (j > 1 && j < dim - 2)) {
          conc_2x4_Boxes(list, 4, i, j);

        }

        else if (i == dim - 2 && j == dim - 2) {
          conc_1x4_Boxes(list, 4, i, j);

        }

        else {
          conc_4x4_Boxes(list, 1, i, j);
          conc_4x4_Boxes(list, 2, i, j);
          conc_4x4_Boxes(list, 3, i, j);
          conc_4x4_Boxes(list, 4, i, j);
        }
      }
    }
  }
}

void Pressure_Lattice::conc_1x4_Boxes(Particle **list, int position, int l,
                                      int m) {
  switch (position) {
  case 1: {
    //~ cout << "i " << i << " "<< " j " << j << endl;
    //~ cout << "l " << l << " "<< " m " << m << endl;
    for (k = 0; k < 4; k++) {
      conc_UppRight_4d_Boxes(k, l, m);
      conc_d_Boxes(list, l, m);
    }
    //~ cout << endl ;
  } break;

  case 2: {
    //~ cout << "i " << i << " "<< " j " << j << endl;
    for (k = 0; k < 4; k++) {
      conc_UppLeft_4c_Boxes(k, l, m);
      conc_c_Boxes(list, l, m);
    }
    //~ cout << endl ;
  } break;

  case 3: {
    //~ cout << "i " << i << " "<< " j " << j << endl;
    for (k = 0; k < 4; k++) {
      conc_LowRight_4b_Boxes(k, l, m);
      conc_b_Boxes(list, l, m);
    }
    //~ cout << endl ;
  } break;

  case 4: {
    //~ cout << "i " << i << " "<< " j " << j << endl;
    for (k = 0; k < 4; k++) {
      conc_LowLeft_4a_Boxes(k, l, m);
      conc_a_Boxes(list, l, m);
    }
    //~ cout << endl ;
  } break;

  default:
    cout << "No box in d_1x4_Boxes" << endl;
  }
}

void Pressure_Lattice::conc_2x4_Boxes(Particle **list, int position, int l,
                                      int m) {
  switch (position) {
  case 1: {
    for (k = 0; k < 4; k++) {
      conc_UppLeft_4c_Boxes(k, l, m);
      conc_c_Boxes(list, l, m);
    }

    for (k = 0; k < 4; k++) {
      conc_UppRight_4d_Boxes(k, l, m);
      conc_d_Boxes(list, l, m);
    }
  } break;

  case 2: {
    for (k = 0; k < 4; k++) {
      conc_LowRight_4b_Boxes(k, l, m);
      conc_b_Boxes(list, l, m);
    }

    for (k = 0; k < 4; k++) {
      conc_UppRight_4d_Boxes(k, l, m);
      conc_d_Boxes(list, l, m);
    }
  } break;

  case 3: {
    for (k = 0; k < 4; k++) {
      conc_LowLeft_4a_Boxes(k, l, m);

      conc_a_Boxes(list, l, m);
    }

    for (k = 0; k < 4; k++) {
      conc_UppLeft_4c_Boxes(k, l, m);
      conc_c_Boxes(list, l, m);
    }
  } break;

  case 4: {
    for (k = 0; k < 4; k++) {
      conc_LowLeft_4a_Boxes(k, l, m);
      conc_a_Boxes(list, l, m);
    }

    for (k = 0; k < 4; k++) {
      conc_LowRight_4b_Boxes(k, l, m);
      conc_b_Boxes(list, l, m);
    }
  } break;

  default:
    cout << "No box in d_2x4_Boxes" << endl;
  }
}

void Pressure_Lattice::conc_4x4_Boxes(Particle **list, int position, int l,
                                      int m) {
  switch (position) {
  case 1: {
    for (k = 0; k < 4; k++) {
      conc_LowLeft_4a_Boxes(k, l, m);
      conc_a_Boxes(list, l, m);
    }
  } break;

  case 2: {
    for (k = 0; k < 4; k++) {
      conc_LowRight_4b_Boxes(k, l, m);
      conc_b_Boxes(list, l, m);
    }
  } break;

  case 3: {
    for (k = 0; k < 4; k++) {
      conc_UppLeft_4c_Boxes(k, l, m);
      conc_c_Boxes(list, l, m);
    }
  } break;

  case 4: {
    for (k = 0; k < 4; k++) {
      conc_UppRight_4d_Boxes(k, l, m);
      conc_d_Boxes(list, l, m);
    }
  } break;

  default:
    cout << "No box in d_4x4_Boxes" << endl;
  }
}

void Pressure_Lattice::conc_LowLeft_4a_Boxes(int pos_count, int q, int r) {
  int res;

  res = (dim - 3) * 2;

  switch (pos_count) {
  case 0:
    box = ((q - 1) * (4 * res)) - (4 * res) - 2 + ((r - 1) * 2);
    break;
  case 1:
    box = ((q - 1) * (4 * res)) - (4 * res) - 1 + ((r - 1) * 2);
    break;
  case 2:
    box = ((q - 1) * (4 * res)) - (2 * res) - 2 + ((r - 1) * 2);
    break;
  case 3:
    box = ((q - 1) * (4 * res)) - (2 * res) - 1 + ((r - 1) * 2);
    break;

  default:
    cout << "No box in d_LowLeft_4a_Boxes " << endl;
  }
}

void Pressure_Lattice::conc_LowRight_4b_Boxes(int pos_count, int q, int r) {
  int res;

  res = (dim - 3) * 2;
  switch (pos_count) {
  case 0:
    box = ((q - 1) * (4 * res)) - (4 * res) + ((r - 1) * 2);
    break;
  case 1:
    box = ((q - 1) * (4 * res)) - (4 * res) + 1 + ((r - 1) * 2);
    break;
  case 2:
    box = ((q - 1) * (4 * res)) - (2 * res) + ((r - 1) * 2);
    break;
  case 3:
    box = ((q - 1) * (4 * res)) - (2 * res) + 1 + ((r - 1) * 2);
    break;

  default:
    cout << "No box in d_LowRight_4b_Boxes " << endl;
  }
}

void Pressure_Lattice::conc_UppLeft_4c_Boxes(int pos_count, int q, int r) {
  int res;

  res = (dim - 3) * 2;

  switch (pos_count) {
  case 0:
    box = ((q - 1) * (4 * res)) - 2 + ((r - 1) * 2);
    break;
  case 1:
    box = ((q - 1) * (4 * res)) - 1 + ((r - 1) * 2);
    break;
  case 2:
    box = ((q - 1) * (4 * res)) + (2 * res) - 2 + ((r - 1) * 2);
    break;
  case 3:
    box = ((q - 1) * (4 * res)) + (2 * res) - 1 + ((r - 1) * 2);
    break;

  default:
    cout << "No box in d_UppLeft_4c_Boxes " << endl;
  }
}

void Pressure_Lattice::conc_UppRight_4d_Boxes(int pos_count, int q, int r) {
  int res;

  res = (dim - 3) * 2;

  switch (pos_count) {
  case 0:
    box = ((q - 1) * (4 * res)) + ((r - 1) * 2);
    break;
  case 1:
    box = ((q - 1) * (4 * res)) + 1 + ((r - 1) * 2);
    break;
  case 2:
    box = ((q - 1) * (4 * res)) + (2 * res) + ((r - 1) * 2);
    break;
  case 3:
    box = ((q - 1) * (4 * res)) + (2 * res) + 1 + ((r - 1) * 2);
    break;

  default:
    cout << "No box in d_UppRight_4d_Boxes " << endl;
  }
}

void Pressure_Lattice::conc_a_Boxes(Particle **list, int n, int o) {
  int res;

  res = (dim - 3) * 2;

  if (list[box]) {
    help = list[box];

    x_diff = abs(help->xpos - Pressure_x[n - 1][o - 1]);
    y_diff = abs(help->ypos - Pressure_y[n - 1][o - 1]);
    smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n - 1][o]);
    y_diff = abs(help->ypos - Pressure_y[n - 1][o]);
    smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n][o - 1]);
    y_diff = abs(help->ypos - Pressure_y[n][o - 1]);
    smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n][o]);
    y_diff = abs(help->ypos - Pressure_y[n][o]);
    smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    conc_1 = (concent[n - 1][o - 1] * smooth_1st);
    conc_2 = (concent[n - 1][o] * smooth_2nd);
    conc_3 = (concent[n][o - 1] * smooth_3rd);
    conc_4 = (concent[n][o] * smooth_4th);

    help->conc = conc_1 + conc_2 + conc_3 + conc_4;

    while (help->next_inBox) {
      x_diff = abs(help->next_inBox->xpos - Pressure_x[n - 1][o - 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n - 1][o - 1]);
      smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n - 1][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n - 1][o]);
      smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o - 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o - 1]);
      smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o]);
      smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      conc_1 = (concent[n - 1][o - 1] * smooth_1st);
      conc_2 = (concent[n - 1][o] * smooth_2nd);
      conc_3 = (concent[n][o - 1] * smooth_3rd);
      conc_4 = (concent[n][o] * smooth_4th);

      help->next_inBox->conc = conc_1 + conc_2 + conc_3 + conc_4;

      help = help->next_inBox;
    }
  }
}

void Pressure_Lattice::conc_b_Boxes(Particle **list, int n, int o) {

  int res;

  res = (dim - 3) * 2;

  if (list[box]) {
    help = list[box];

    x_diff = abs(help->xpos - Pressure_x[n - 1][o]);
    y_diff = abs(help->ypos - Pressure_y[n - 1][o]);
    smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n - 1][o + 1]);
    y_diff = abs(help->ypos - Pressure_y[n - 1][o + 1]);
    smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n][o]);
    y_diff = abs(help->ypos - Pressure_y[n][o]);
    smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n][o + 1]);
    y_diff = abs(help->ypos - Pressure_y[n][o + 1]);
    smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    conc_1 = (concent[n - 1][o] * smooth_1st);
    conc_2 = (concent[n - 1][o + 1] * smooth_2nd);
    conc_3 = (concent[n][o] * smooth_3rd);
    conc_4 = (concent[n][o + 1] * smooth_4th);

    help->conc = conc_1 + conc_2 + conc_3 + conc_4;

    while (help->next_inBox) {
      x_diff = abs(help->next_inBox->xpos - Pressure_x[n - 1][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n - 1][o]);
      smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n - 1][o + 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n - 1][o + 1]);
      smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o]);
      smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o + 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o + 1]);
      smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      conc_1 = (concent[n - 1][o] * smooth_1st);
      conc_2 = (concent[n - 1][o + 1] * smooth_2nd);
      conc_3 = (concent[n][o] * smooth_3rd);
      conc_4 = (concent[n][o + 1] * smooth_4th);

      help->next_inBox->conc = conc_1 + conc_2 + conc_3 + conc_4;

      help = help->next_inBox;
    }
  }
}

void Pressure_Lattice::conc_c_Boxes(Particle **list, int n, int o) {
  int res;

  res = (dim - 3) * 2;

  if (list[box]) {
    help = list[box];

    x_diff = abs(help->xpos - Pressure_x[n][o - 1]);
    y_diff = abs(help->ypos - Pressure_y[n][o - 1]);
    smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n][o]);
    y_diff = abs(help->ypos - Pressure_y[n][o]);
    smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n + 1][o - 1]);
    y_diff = abs(help->ypos - Pressure_y[n + 1][o - 1]);
    smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n + 1][o]);
    y_diff = abs(help->ypos - Pressure_y[n + 1][o]);
    smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    conc_1 = (concent[n][o - 1] * smooth_1st);
    conc_2 = (concent[n][o] * smooth_2nd);
    conc_3 = (concent[n + 1][o - 1] * smooth_3rd);
    conc_4 = (concent[n + 1][o] * smooth_4th);

    help->conc = conc_1 + conc_2 + conc_3 + conc_4;

    while (help->next_inBox) {
      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o - 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o - 1]);
      smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o]);
      smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n + 1][o - 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n + 1][o - 1]);
      smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n + 1][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n + 1][o]);
      smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      conc_1 = (concent[n][o - 1] * smooth_1st);
      conc_2 = (concent[n][o] * smooth_2nd);
      conc_3 = (concent[n + 1][o - 1] * smooth_3rd);
      conc_4 = (concent[n + 1][o] * smooth_4th);

      help->next_inBox->conc = conc_1 + conc_2 + conc_3 + conc_4;

      help = help->next_inBox;
    }
  }
}

void Pressure_Lattice::conc_d_Boxes(Particle **list, int n, int o) {
  //~ cout << "p " << p << " "<< " k " << k << endl;

  int res;

  res = (dim - 3) * 2;

  if (list[box]) {
    help = list[box];
    //~ cout << box << " " << "xpos " << help->xpos << " " ;

    x_diff = abs(help->xpos - Pressure_x[n][o]);
    y_diff = abs(help->ypos - Pressure_y[n][o]);
    smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n][o + 1]);
    y_diff = abs(help->ypos - Pressure_y[n][o + 1]);
    smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n + 1][o]);
    y_diff = abs(help->ypos - Pressure_y[n + 1][o]);
    smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    x_diff = abs(help->xpos - Pressure_x[n + 1][o + 1]);
    y_diff = abs(help->ypos - Pressure_y[n + 1][o + 1]);
    smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

    conc_1 = (concent[n][o] * smooth_1st);
    conc_2 = (concent[n][o + 1] * smooth_2nd);
    conc_3 = (concent[n + 1][o] * smooth_3rd);
    conc_4 = (concent[n + 1][o + 1] * smooth_4th);

    help->conc = conc_1 + conc_2 + conc_3 + conc_4;

    while (help->next_inBox) {
      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o]);
      smooth_1st = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n][o + 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n][o + 1]);
      smooth_2nd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n + 1][o]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n + 1][o]);
      smooth_3rd = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      x_diff = abs(help->next_inBox->xpos - Pressure_x[n + 1][o + 1]);
      y_diff = abs(help->next_inBox->ypos - Pressure_y[n + 1][o + 1]);
      smooth_4th = (1 - (x_diff / delta_x)) * (1 - (y_diff / delta_y));

      conc_1 = (concent[n][o] * smooth_1st);
      conc_2 = (concent[n][o + 1] * smooth_2nd);
      conc_3 = (concent[n + 1][o] * smooth_3rd);
      conc_4 = (concent[n + 1][o + 1] * smooth_4th);

      help->next_inBox->conc = conc_1 + conc_2 + conc_3 + conc_4;

      help = help->next_inBox;
    }
  }
}

void Pressure_Lattice::setDarcyVelocity(Particle **list) // used now
{
  int res;

  res = (dim - 3) * 2;

  for (i = 1; i < dim - 2; i++) {
    for (j = 1; j < dim - 2; j++) {

      particle_in_box = 0;

      vel_darcyx[i][j] = 0;
      vel_darcyy[i][j] = 0;

      for (k = 0; k < 4; k++) // note that fluid lattice is dim + 3 - empty rows
                              // around lattice
      {
        if (k == 0)
          box = ((i - 1) * 2 * res) +
                ((j - 1) * 2); // res = resolution of particle lattice
        else if (k == 1)
          box = ((i - 1) * 2 * res) + ((j - 1) * 2) + 1;
        else if (k == 2)
          box = ((i - 1) * 2 * res) + (2 * res) + ((j - 1) * 2);
        else if (k == 3)
          box = ((i - 1) * 2 * res) + (2 * res) + ((j - 1) * 2) + 1;

        if (list[box]) {
          help = list[box];

          vel_darcyx[i][j] = vel_darcyx[i][j] + help->v_darcyx;
          vel_darcyy[i][j] = vel_darcyy[i][j] + help->v_darcyy;

          particle_in_box++;

          while (help->next_inBox) {
            vel_darcyx[i][j] = vel_darcyx[i][j] + help->next_inBox->v_darcyx;
            vel_darcyy[i][j] = vel_darcyy[i][j] + help->next_inBox->v_darcyy;

            particle_in_box++;

            help = help->next_inBox;
          }
        }
      }

      vel_darcyx[i][j] = vel_darcyx[i][j] / particle_in_box;
      vel_darcyy[i][j] = vel_darcyy[i][j] / particle_in_box;
      ;
    }
  }
}

void Pressure_Lattice::get_concentration(Particle **list) {
  int res;

  res = (dim - 3) * 2;

  for (i = 0; i < dim - 3; i++) {
    for (j = 0; j < dim - 3; j++) {

      for (k = 0; k < 4; k++) // note that fluid lattice is dim + 3 - empty rows
                              // around lattice
      {
        if (k == 0)
          box = (j * 2 * res) + (i * 2); // res = resolution of particle lattice
        else if (k == 1)
          box = (j * 2 * res) + (i * 2) + 1;
        else if (k == 2)
          box = (j * 2 * res) + (2 * res) + (i * 2);
        else if (k == 3)
          box = (j * 2 * res) + (2 * res) + (i * 2) + 1;

        if (list[box]) {
          help = list[box];

          help->conc = concent[j][i];

          while (help->next_inBox) {
            help->next_inBox->conc = concent[j][i];

            help = help->next_inBox;
          }
        }
      }
    }
  }
}

void Pressure_Lattice::SolveadvDisExplicitTransient(double fluid_v,
                                                    double dif_const) {
  // remember the x is j and y is i

  double D, dt, dx, dy, dx2, dy2;

  // D=0.000001;

  D = dif_const;

  dx = delta_x; //
  dy = delta_y; //
  dx2 = dx * dx;
  dy2 = dy * dy;
  // dt=0.1*dx2;//0,5*(dx)^2 for diffusion at least but less also possible
  dt = 0.002;

  // velocity darcy
  // for(i = 0; i < dim; i++)
  //{
  // for(j =0; j < dim; j++)
  //{
  // vel_darcyy[i][j]=10000000000000*kappa[i][j];
  // vel_darcyx[i][j]=0*kappa[i][j];
  // }
  // }

  // boundary condition (background concentration is already given)

  for (v = 0; v < dim; v++) {
    // concent[v][0]=0.;
    concent[0][v] = 0.5; // setting the lower boundary
    // concent[v][dim]=0.;
    concent[dim][v] = 0;
  }
  // for (v = 5; v < dim-5; v++)
  //{
  // concent[v][0]=0.;
  //	concent[0][v]=0.25; // setting the lower boundary
  // concent[v][dim]=0.;
  // concent[dim][v]=0;
  //}
  // for (v = 0; v < dim; v++)

  {
    // concent[50][50]=0.5;
    // concent[0][v]=0.5; // setting the lower boundary
    // concent[v][dim]=0.;
    // concent[dim][v]=0;
  }
  for (i = 1; i < dim - 1; i++) // y loop
  {
    for (j = 3; j < dim - 3; j++) // x loop
    {

      ddxc[i][j] = (concent[i + 1][j] - 2 * concent[i][j] + concent[i - 1][j]) /
                   dy2; // x has to vary in y, so takes neighbours in i (=y)
      ddyc[i][j] = (concent[i][j + 1] - 2 * concent[i][j] + concent[i][j - 1]) /
                   dx2; // y has to vary in x, so takes neighbours in j (=x)

      // forward approximation takes the difference between node and node behind
      // (with respect to velocity) and adds to node

      if (vel_darcyx[i][j] >= 0.0) // x direction positive, to the right
      {
        dxc[i][j] = (concent[i][j] - concent[i][j - 1]) / dx;
      } else // x direction negative, to the left
      {
        dxc[i][j] = (concent[i][j] - concent[i][j + 1]) / dx;
      }
      if (vel_darcyy[i][j] >= 0.0) // y direction positive, upwards
      {
        dyc[i][j] = (concent[i][j] - concent[i - 1][j]) / dy;
      } else // y direction positive, downwards
      {
        dyc[i][j] = (concent[i][j] - concent[i + 1][j]) / dy;
        // cout << "negative" << endl;
        // dyc[i][j] = 0.0;
      }
    }
  }

  // you now have all the values to fill in the matrix and update the
  // concentration. You have to do that after the calculation, otherwise you get
  // a trend in the calculation depending on where you start because some have
  // been updated and others not. Need to update in on step.

  for (i = 1; i < dim - 1; i++) // y
  {
    for (j = 3; j < dim - 3; j++) // x
    {
      // diffusion
      concent[i][j] = concent[i][j] + dt * D * (ddxc[i][j] + ddyc[i][j]);
      // advection
      if (vel_darcyy[i][j] > 0.00000001) {
        concent[i][j] =
            concent[i][j] - (dt * vel_darcyx[i][j] * fluid_v * dxc[i][j]);
        concent[i][j] =
            concent[i][j] - (dt * vel_darcyy[i][j] * fluid_v * dyc[i][j]);
      } else if (vel_darcyy[i][j] < -0.00000001) {
        concent[i][j] =
            concent[i][j] - (dt * vel_darcyx[i][j] * fluid_v * dxc[i][j]);
        concent[i][j] =
            concent[i][j] - (dt * vel_darcyy[i][j] * fluid_v * dyc[i][j]);
      }

      // original

      // concent[i][j] = concent[i][j]-(dt*1000000000000*kappa[i][j]*dyc[i][j]);
    }
  }

  for (v = 1; v < dim - 1;
       v++) // boundary condition for sides, not really wrapping... sort of
  {
    concent[v][2] = concent[v][3];
    concent[v][1] = concent[v][2];
    concent[v][0] = concent[v][1];

    concent[v][dim - 3] = concent[v][dim - 4];
    concent[v][dim - 2] = concent[v][dim - 3];
    concent[v][dim - 1] = concent[v][dim - 2];
  }
}

void Pressure_Lattice::ChangeConcentration(int box_x, int box_y,
                                           double change) {
  concent[box_y][box_x] = concent[box_y][box_x] - change;
  if (concent[box_y][box_x] < 0.0)
    concent[box_y][box_x] = 0.0;
}
