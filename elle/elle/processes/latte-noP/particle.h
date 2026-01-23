/******************************************************
 * Spring Code Latte 2.0
 *
 * Class Particle in particle.h
 *
 * Daniel Koehn and Jochen Arnold Feb. 2002
 * Oslo, Mainz
 *
 * Daniel Dec. 2003
 *
 * We thank Anders Malthe-Soerenssen for his enormous help
 * and introduction to these codes
 *
 * new Daniel  2004/5
 * koehn_uni-mainz.de
 ******************************************************/
#ifndef _E_particle_h
#define _E_particle_h

#include "unodesP.h" // include Unode class
#include <cmath>

/******************************************************
 * the particle class for the particles
 * contains pointers to neighbours and springs
 * contains function that can relax single particles
 * and function that can connect particle with neighbours
 * called from the lattice constructor
 *
 * Daniel spring 2002
 *******************************************************/

class Particle {
public:
  //------------------------------------------------------------
  // variables
  //------------------------------------------------------------

  double angle_of_internal_friction;
  int nb;        // my number
  double radius; // my radius
  double average_radius;
  double rad[8];      // different initial length for different springs
  double newrad[8];   // second length for springs
  int neig_spring[8]; // find backwards pointers of neighbours
  int boundary;       // what is my boundary if any (not yet active)
  double mbreak;      // what is my max breaking probability
  int neigh;          // counter for which neighbour breaks easiest
  bool fix_x;         // am I fixed in x ?
  bool fix_y;         // am I fixed in y ?

  double neigpos[8][2]; // equilibrium position of the neighbour particles

  float rad_par_fluid;
  float frac_rad_par_fluid;
  double area_par_fluid;

  bool novisc;
  int fracture_time;
  double fluid_pressure_gradient;
  double average_porosity;
  double average_permeability;

  double disconst;

  double crack_seal_counter;

  int healing;

  int cluster_nb;

  double age;

  //----------------------------------
  // Seismic Moment Tensor stuff
  //----------------------------------

  double movement_vector_n[8][2];
  double movement_vector_s[8][2];
  double movement_normal[8];
  double movement_shear[8];

  double fault_vector[8][2];
  double fault_midpoint[8][2];

  //---------------------------------

  double ratio_a, add[9], ratio[8];

  int fluid_particle[8];     // fluid particles in lattice gas
  int fluid_particle_new[8]; // fluid particles in lattice gas
  int fluid_particles;       // counter for fluid particles in lattice gas

  bool nofluid; // bool for fracture walls

  double rep_rad[9]; // repulsion radius

  double viscosity; // viscosity of particle
  double previous_viscosity;
  bool is_moho_particle;

  // for conjugate radii (see viscosity)
  double ax, ay, bx, by;

  bool neig_of_two_spinel_grains;
  double normal_stress[6];
  int react_time;
  double xx, yy;
  double increase_rate;
  double rate_factor;
  double spring_var;
  double dis_time;
  double gravforce;

  double SheetT;
  double Tdif;

  int particle_age;

  double v_darcyx, v_darcyy;

  double previous_rad[6];

  bool spring_boundary[8];

  bool is_boundary;         // particle of grain boundary
  bool is_lattice_boundary; // particle of box boundary
  bool is_left_lattice_boundary;
  bool is_right_lattice_boundary;
  bool is_top_lattice_boundary;
  bool is_bottom_lattice_boundary;

  float right_lattice_wall_pos;
  float left_lattice_wall_pos;

  bool no_break[8];       // not allowed to break
  bool is_phase_boundary; // is a phase boundary
  bool isUpperBox;        // upper part of box
  bool isLowerBox;        // lower part of box

  int elle_Node[32]; // Elle nodes that are close by

  int box_pos; // pos in repulsion box
  int damage;
  float fracture_reaction;
  bool done; // already added in force list (direct neighbour)

  bool isHole;                 // flag for Holes
  bool isHoleBoundary;         // flag for hole boundaries
  bool isHoleInternalBoundary; // flag for growth at boundaries

  double area;        // area of particle for reactions
  double mV;          // Molecular Volume
  double surf_free_E; // interfacial free energy
  double surfE;       // surface energy
  double newsurfE;    // mean surface energy across surface

  double prob; // probability to react

  double sxx; // stress in x
  double syy; // stress in y
  double sxy; // shear stress in xy

  double exx; // strain in x
  double eyy; // strain in y
  double exy; // strain in xy
  double e1;  // principal strain 1
  double e2;  // principal strain 2
  double ev;  // volumatric strain

  // weight force of particle
  double fg;

  double average_pressure, local_max;

  double calcite; // calcite-dolomite ratio

  double Fe;
  double Mg;
  double Fe_fluid;
  double Mg_fluid;
  double Fluid;

  double aperture;

  double bound_stress; // stress on boundary
  int mineral;         // number for mineral
  int phase;
  // viscous sheet
  double sheet_x, sheet_y;
  float sheet_young, sheet_visc;
  int adjust_area_p;

  double springconst;
  double springf[8];   // normal force spring constants
  double springs[8];   // shear  force spring constant
  double springv[8];   // my viscosities
  double break_Str[8]; // spring breaking strength (in normal stress)
  double abreak;
  int grain;          // number for my grain
  bool break_Sp[8];   // I break
  bool draw_break;    // extra flag for drawing particles broken
  Particle *neigP[8]; // pointers to eight possible neighbours

  Particle
      *neigP2[8]; // pointer to neighbours, not influenced by breaking springs

  Particle *fluid_neigP[8]; // neighbour list along fluid solid interface
  int fluid_neig_count;     // counter for neighbours
  double young;             // youngs modulus

  double E;             // the 'real' Young's modulus
  double poisson_ratio; // the 'real' Poisson ratio

  int open_springs; // how many open springs for surface energies
  double fluid_P;   // fluid pressure
  double spin_young, spin_springf[6], olivine_young, olivine_springf[6];
  double movement;

  double real_density, real_young, real_radius;

  Particle *leftNeighbour; // left neighbour along interface (looking at fluid)
  Particle
      *rightNeighbour; // right neighbour along interface (looking at fluid)
  Particle *nextP;     // pointer to next in list
  Particle *prevP;     // pointer backwards in list
  Particle *neig;      // help pointer
  Particle *neig2;     // second help pointer
  Particle *neig3;     // third help pointer

  Particle *next_inBox; // in case there is more than one

  Unode *p_Unode; // connected Unode pointer

  double xpos; // my x position
  double ypos; // my y position

  double
      velx; // velocity-x (v = dx/t, t ~= 1) of particles after Pressure class
  double
      vely; // velocity-y (v = dy/t, t ~= 1) of particles after Pressure class

  double F_P_x; // pressure force on particle in x-direction by Pressure class
  double F_P_y; // pressure force on particle in y-direction by Pressure class

  double oldxpos; // old x position
  double oldypos; // old y position

  double oldx; // second old x position
  double oldy; // second old y position
  double initialy;

  short i; // external counter ?

  // neue Variablen

  double conc;      // concentration
  double new_conc;  // new concentration
  double left_conc; // left over  concentration
  double metal_conc;

  // new variables for phase change
  bool potential_nucleus;
  bool gbm_condition;
  bool nucleation_condition;
  bool already_done;
  bool phase_change_required;
  bool transformation_required;
  bool marker;

  double eEl_volume;
  double eEl_distortion;
  double eEl; // elastic Energy

  double time;
  double deltaG;
  double temperature;
  double deltaG_nucleation;
  double deltaG_gbm;
  double deltaG_d;
  double deltaG_t;
  double mean_stress;

  double xy[6][2];

  // new variables for sorting nuclei
  bool local_maximum;

  // new variables to remember previous settings
  int previous_mineral;

  Particle *fixed_neig;

  double previous_prob;
  double previous_mV;
  double previous_temperature;
  double previous_young;
  double previous_deltaG;
  double previous_real_radius;
  double previous_radius;
  double previous_area;
  double last_prob;
  double previous_eEl;
  double previous_surfE;
  double previous_mean_stress;
  double previous_sigman;
  double previous_eEl_volume;
  double previous_eEl_distortion;
  double previous_boundary_state;
  double previous_phase_boundary_state;
  double delta_eEl_distortion;
  double delta_pV;
  double delta_surfE;
  double delta_eEl_volume;
  double delta_Ed; // strain-energy-change without volume-change
  double Ed;
  double delta_Ev;
  double boundary_angle[6];
  double spinel_content;
  double previous_spinel_content;
  double rate;
  double phi; // angle between boundary and sigman
  double previous_eEl_dilatation;

  bool previously_transformed;
  bool nucleus;
  bool right_now_transformed;

  double box_nb[6]; // associates springs to box_rad radii, to speed up
                    // repulsion

  double dl[9]; // contains the difference between alen&dd in terms of
                // boxrad-numbering
  double rep[9]; // saves the according rep_constant from the repulsion-function
  double spring_angle[6];
  double area_factor;

  double smallest_rad;
  int surfcount;

  double maxaxisangle, minaxisangle, maxaxislength, minaxislength;

  // double   alen_list[88];
  int merk;
  // Particle *neigh_list[88];

  bool visc_flag;      // debugging
  bool no_spring_flag; // debugging

  double rot_angle; // the rotation angle for beams
  double a, b, fixangle, maxis1angle, maxis2angle, maxis1length, maxis2length;

  double dir_rate[6]; // rate with different directions for different
                      // neighbour-particles
  double l[4]; // sizes of the square

  double save_rate[4];

  bool no_reaction;

  double mole_per_particle;

  double previous_springf[8], eEl_dilatation, original_springf,
      previous_original_springf, delta_eEl_dilatation, eEl_dif, prev_viscosity,
      prev_xpos, prev_ypos, prev_sxx, prev_syy, prev_sxy;

  double newly_overgrown;
  double overgrown_area;

  int previous_grain;

  bool transform;

  bool transformed_in_last_timestep;
  int transformation_timestep;

  bool wrapping;

  bool no_nucleation;

  //------------------------------------------------
  // functions
  //------------------------------------------------

  Particle();    // constructor
  ~Particle(){}; // destructor

  /************************************************
   * main relaxation routines
   *****************************************************/

  // the fast relaxation routine

  float GetXDist(Particle *center, Particle *peripher);

  bool Relax(double relaxthresh, Particle **list, int size, bool wall,
             double rightwall, double leftwall, double lowerwall,
             double upperwall, double wallconstant, bool length, int debug_nb,
             int visc_rel, bool sheet, bool grav, bool grav_Press, int part_nr,
             int numb);

  bool
  Relax(Particle **list, int size, int visc_rel, bool wall, double rightwall,
        double leftwall, double lowerwall, double upperwall,
        double wallconstant); // end relaxroutine, calculates stress and break

  double BreakBond(); // breaks a bond

  void ChangeBox(Particle **list,
                 int size); // in case particle out of box after move

  /*************************************
   * build lattice routines
   **************************************/

  Particle *FindNeighbourF(int neigDist); // help routine for Connect

  Particle *FindNeighbourB(int neigDist); // help routine for Connect

  void Connect(int lType, int lParticlex, int lParticley,
               bool wrapping); // intial connections

  void SetSprings(double E, double poisson);

  //****************************
  // viscoelastic code
  //****************************

  double Alen(int, double, double); // change spring lengths of particles

  void debug();

  void Tilt(); // tilt eliptical particles

  double Add_Gravitation();
  double Press_gravitation(int curr_nr, int thresh_num);

  bool use_gravitation;
};
#endif
