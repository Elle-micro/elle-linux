#include "phreeqc.h"
using namespace std;
using std::vector;
PhreeqC::PhreeqC()

{
  this->iphreeqc;
}

template <typename T> string ConvertToString(T value) {
  std::ostringstream strs;
  strs << value;
  std::string str = strs.str();
  return (str);
}

void setVector(vector<const char *> &phreeqc_elements) {
  // this vector contains all possible species and parameters phreeqc can read
  // we will use this to verify the input in the txt files
  phreeqc_elements.push_back("pH");
  phreeqc_elements.push_back("pe");
  phreeqc_elements.push_back("temp");
  phreeqc_elements.push_back("density");
  phreeqc_elements.push_back("Alkalinity");

  phreeqc_elements.push_back("Al");
  phreeqc_elements.push_back("Ba");
  phreeqc_elements.push_back("B");
  phreeqc_elements.push_back("Br");
  phreeqc_elements.push_back("Cd");
  phreeqc_elements.push_back("Ca");
  phreeqc_elements.push_back("C");
  phreeqc_elements.push_back("C(4)");
  phreeqc_elements.push_back("C(-4)");
  phreeqc_elements.push_back("Cl");
  phreeqc_elements.push_back("Cu");
  phreeqc_elements.push_back("Cu(2)");
  phreeqc_elements.push_back("Cu(1)");
  phreeqc_elements.push_back("F");
  phreeqc_elements.push_back("H(0)");
  phreeqc_elements.push_back("Fe");
  phreeqc_elements.push_back("Fe(2)");
  phreeqc_elements.push_back("Fe(3)");
  phreeqc_elements.push_back("Pb");
  phreeqc_elements.push_back("Li");
  phreeqc_elements.push_back("Mg");
  phreeqc_elements.push_back("Mn");
  phreeqc_elements.push_back("Mn(2)");
  phreeqc_elements.push_back("Mn(3)");
  phreeqc_elements.push_back("N");
  phreeqc_elements.push_back("N(5)");
  phreeqc_elements.push_back("N(3)");
  phreeqc_elements.push_back("N(0)");
  phreeqc_elements.push_back("N(-3)");
  phreeqc_elements.push_back("O(0)");
  phreeqc_elements.push_back("P");
  phreeqc_elements.push_back("K");
  phreeqc_elements.push_back("Si");
  phreeqc_elements.push_back("Na");
  phreeqc_elements.push_back("Sr");
  phreeqc_elements.push_back("S");
  phreeqc_elements.push_back("S(6)");
  phreeqc_elements.push_back("S(-2)");
  phreeqc_elements.push_back("Zn");
}

void PhreeqC::UpdateVector(VAR v,
                           vector<std::pair<const char *, double>> &vec) {
  Selected(v);

  FLUID ::iterator it;
  for (it = vec.begin(); it != vec.end(); it++) {
    if (it->first != "Alkalinity")
      it->second = ReturnValue(v, it->first);
    else
      it->second = ReturnValue(v, "Alk(eq/kgw)");
  }
}

// for debugging. This prints everything in v
void PhreeqC::Selected(VAR v) {
  // Output of Selected output
  std::cout << "selected-output:" << std::endl;
  for (int i = 0; i < iphreeqc.GetSelectedOutputRowCount(); ++i) {
    for (int j = 0; j < iphreeqc.GetSelectedOutputColumnCount(); ++j) {
      if (iphreeqc.GetSelectedOutputValue(i, j, &v) == VR_OK) {
        switch (v.type) {
        case TT_LONG:
          std::cout << v.lVal << " ";
          break;
        case TT_DOUBLE:
          std::cout << v.dVal << " " << "(" << i << " " << j << ")" << " ";
          break;
        case TT_STRING:
          std::cout << v.sVal << " ";
          break;
        }
      }
    }
    std::cout << std::endl;
  }
}

double PhreeqC::ReturnValue(VAR vs, string name) {
  // note this returns only the requested value in the last row of the vector
  double value = 0;
  for (int i = 0; i < iphreeqc.GetSelectedOutputRowCount(); ++i) {
    for (int j = 0; j < iphreeqc.GetSelectedOutputColumnCount(); ++j) {
      if (iphreeqc.GetSelectedOutputValue(i, j, &vs) == VR_OK) {
        switch (vs.type) {
        case TT_STRING:
          string S = vs.sVal;
          if (S.find(name) != std::string::npos) {
            if (iphreeqc.GetSelectedOutputValue(
                    iphreeqc.GetSelectedOutputRowCount() - 1, j, &vs) == VR_OK)
              value = vs.dVal;
            else
              ERROR();
          }
          break;
        }
      }
    }
  }
  return (value);
}

// test the phreeqc database---------------------------------------------
int PhreeqC::loader(const char *d_base) {
  // Try provided path first
  if (iphreeqc.LoadDatabase(d_base) != 0) {
#ifdef LATTE_DATA_DIR
    // If relative name failed, try resolving under LATTE_DATA_DIR
    std::string alt =
        std::string(LATTE_DATA_DIR) + std::string("/") + std::string(d_base);
    if (iphreeqc.LoadDatabase(alt.c_str()) == 0) {
      return EXIT_SUCCESS;
    }
#endif
    iphreeqc.OutputErrorString(); // print error message from phreeqc
    exit(EXIT_FAILURE);           // exit
  }
  return EXIT_SUCCESS;
}

// just an error handler-------------------------------------------------
void PhreeqC::ERROR() {
  iphreeqc.OutputAccumulatedLines(); // print lines which caused the crash
  iphreeqc.OutputErrorString();      // print error message from phreeqc
  exit(EXIT_FAILURE);                // exit the caluculation
}

//=====================basic============================================
// Get version and max threads-------------------------------------------
void PhreeqC::Version() {
  // OpenMP is optional; only query threads when compiled with OpenMP
#ifdef _OPENMP
  nb_threads = omp_get_max_threads();
#else
  nb_threads = 1;
#endif
  static const char *Vers = iphreeqc.GetVersionString();
  cout << "   Iphreeqc-version-" << Vers << endl;
  cout << "	Number of threads: " << nb_threads << endl;
}

// set the database
void PhreeqC::SetDatabase(const char *database) {
  loader(database);
  p_database = database;
}

// set the file containing the pore fluid
void PhreeqC::SetPoreFluid(const char *f1) {
  fstream fileStream;
  fileStream.open(f1);
  if (fileStream.fail()) {
    // Try again under LATTE_DATA_DIR if available
#ifdef LATTE_DATA_DIR
    std::string alt =
        std::string(LATTE_DATA_DIR) + std::string("/") + std::string(f1);
    fileStream.clear();
    fileStream.open(alt.c_str());
    if (!fileStream.fail()) {
      pore_fluid = strdup(alt.c_str());
    } else
#endif
    {
      cerr << "Opening " << f1 << " failed!" << endl;
      exit(EXIT_FAILURE);
    }
  } else {
    pore_fluid = f1;
  }
}

// set the file containng infiltrating fluid
void PhreeqC::SetInfFluid(const char *f2) {
  fstream fileStream;
  fileStream.open(f2);
  if (fileStream.fail()) {
    // Try again under LATTE_DATA_DIR if available
#ifdef LATTE_DATA_DIR
    std::string alt =
        std::string(LATTE_DATA_DIR) + std::string("/") + std::string(f2);
    fileStream.clear();
    fileStream.open(alt.c_str());
    if (!fileStream.fail()) {
      inf_fluid = strdup(alt.c_str());
    } else
#endif
    {
      cerr << "Opening " << f2 << " failed!" << endl;
      exit(EXIT_FAILURE);
    }
  } else {
    inf_fluid = f2;
  }
}

//==========================Initialisation==============================
void PhreeqC::InitializeFluidVector() {
  vector<const char *> phreeqc_elements;
  vector<const char *> simulation_elements;

  setVector(phreeqc_elements);
  vector<const char *> FILES;
  FILES.push_back(pore_fluid);
  FILES.push_back(inf_fluid);

  for (int F = 0; F < 2; F++) {
    int n = 0;
    string str, y;
    size_t nPos = 0, nPos2 = 0;
    ifstream file(FILES[F]);

    if (file.is_open()) {
      while (getline(file, str)) {
        const char *c = str.c_str();

        // skip lines that start with whitespace or are empty
        if (!str.empty() && !isspace(static_cast<unsigned char>(c[0]))) {
          size_t nPos = str.find_first_of(' ');
          size_t nPos2 = str.find_first_of('	');
          size_t POS;

          if (nPos < nPos2) {
            y = (nPos != string::npos ? str.substr(0, nPos) : str);
            POS = nPos;
          }
          if (nPos > nPos2) {
            y = (nPos2 != string::npos ? str.substr(0, nPos2) : str);
            POS = nPos2;
          }

          if (n == 0) {
            if (y != "units") {
              cout << "ERROR: units not defined in input file " << FILES[F]
                   << endl;
              exit(EXIT_FAILURE);
            } else {
              if (units.length() == 0)
                units = str.substr(POS);
              else if (!boost::iequals(units, str.substr(POS))) {
                cout << "ERROR: Units do not match in input files " << units
                     << " " << str.substr(POS) << endl;
                exit(EXIT_FAILURE);
              }
              n++;
            }
          } else {
            for (int e = 0; e < phreeqc_elements.size(); e++) {
              if (boost::iequals(y, phreeqc_elements[e])) {
                if (std::find(simulation_elements.begin(),
                              simulation_elements.end(),
                              phreeqc_elements[e]) != simulation_elements.end())
                  continue;
                else
                  simulation_elements.push_back(phreeqc_elements[e]);
              }
            }
          }
        }
      }
    } else
      cout << "ERROR: cannot open file " << FILES[F] << endl;
  }
  // create a vector of pairs containing the parameters of the phreeqc input and
  boost::array<double, dimX * dimY> conc_1D = {0};
  Conc_matrices.reserve(simulation_elements.size());
  Conc_inf.reserve(simulation_elements.size());
  for (int i = 0; i < simulation_elements.size(); i++) {
    Conc_matrices.push_back(make_pair(simulation_elements.at(i), conc_1D));
    Conc_inf.push_back(make_pair(simulation_elements.at(i), 0));
    Conc_por.push_back(make_pair(simulation_elements.at(i), 0));
  }
  int i;
  string str;
  ifstream file(pore_fluid);
  if (file.is_open()) {
    getline(file, str);
    while (getline(file, str)) {
      istringstream is(str);
      double value;
      string name;
      is >> name >> value;
      for (FLUID ::iterator it = Conc_por.begin(); it != Conc_por.end(); it++) {
        if (boost::iequals(it->first, name))
          (*it).second = value;
      }
    }
  }
  file.close();
  ifstream file2(inf_fluid);
  if (file2.is_open()) {
    getline(file2, str);
    while (getline(file2, str)) {
      istringstream is(str);
      double value;
      string name;
      is >> name >> value;
      for (FLUID ::iterator it = Conc_inf.begin(); it != Conc_inf.end(); it++) {
        if (boost::iequals(it->first, name))
          (*it).second = value;
      }
    }
  }
  file2.close();
  AdjustUnits(Conc_por);
  AdjustUnits(Conc_inf);
  units = "mol/kgw";
}

void PhreeqC::SetFluidMatrices(vector<std::pair<const char *, double>> &fluid,
                               int boundary) {
  cout << "Setting matrices" << endl;
  for (FLUID ::iterator itt = fluid.begin(); itt != fluid.end(); itt++) {
    const char *name = itt->first;
    double value = itt->second;
    double value2 = FindConcentration(name, Conc_inf);

    cout << "setting value for " << name << " " << value << " and " << value2
         << endl;
    for (FLUID_array ::iterator it = Conc_matrices.begin();
         it != Conc_matrices.end(); it++) {
      if (boost::iequals(it->first, name)) {
        boost::array<double, dimX * dimY> conc_1D;
        for (int i = 0; i < dimX; i++) {
          for (int j = 0; j < dimY; j++) {
            switch (boundary) {
            case 0:
              if (j == 0)
                conc_1D[i * dimY + j] = value2;
              else
                conc_1D[i * dimY + j] = value;
              break;

            case 1:
              if (j < 2)
                conc_1D[i * dimY + j] = value2;
              else
                conc_1D[i * dimY + j] = value;
              break;

            case 2:
              if (i == dimX)
                conc_1D[i * dimY + j] = value2;
              else
                conc_1D[i * dimY + j] = value2;
              break;

            case 3:
              if (j == dimY)
                conc_1D[i * dimY + j] = value2;
              else
                conc_1D[i * dimY + j] = value2;
              break;
            default:
              cout << "Boundary not defined (0,1,2 or 3)";
              exit(EXIT_FAILURE);
            }
          }
        }
        it->second = conc_1D;
      }
    }
  }
}

void PhreeqC::AddPoreFluidLayer(vector<std::pair<const char *, double>> fluid,
                                int y_min, int y_max) {
  cout << "Adding pore fluid layer at " << y_min << " " << y_max << endl;
  for (FLUID ::iterator itt = fluid.begin(); itt != fluid.end(); itt++) {
    const char *name = itt->first;
    double value2 = FindConcentration(name, fluid);

    for (FLUID_array ::iterator it = Conc_matrices.begin();
         it != Conc_matrices.end(); it++) {
      if (boost::iequals(it->first, name)) {
        boost::array<double, dimX * dimY> conc_1D;
        for (int i = 0; i < dimX; i++) {
          for (int j = 0; j < dimY; j++) {
            if (j > y_min && j < y_max)
              conc_1D[i * dimY + j] = value2;
            else
              conc_1D[i * dimY + j] = it->second[i * dimY + j];
          }
        }
        it->second = conc_1D;
      }
    }
  }
}

double PhreeqC::TermalDiffusivity(double thermal_conductivity,
                                  double heat_capacity, double density) {
  double alpha = thermal_conductivity / (heat_capacity * density);
  return (alpha);
}

double PhreeqC::Diffusivity(double D, double Ea, double T) {
  double R = 8.314462;
  double d = D * exp(Ea / (T * R));
  return (d);
}

void PhreeqC::ChargeBalance(vector<std::pair<const char *, double>> &fluid) {
  // obtain charge balace for given fluid compositionx
  VAR vs;
  std::ostringstream oss;
  oss << "Solution 1" << endl;
  oss << "units " << units.c_str() << endl;

  for (FLUID ::iterator it = fluid.begin(); it != fluid.end(); it++) {
    if (it->first == "pH")
      oss << it->first << " " << ConvertToString(it->second).c_str()
          << " charge" << endl;
    else
      oss << it->first << " " << ConvertToString(it->second).c_str() << endl;
  }
  oss << "SELECTED_OUTPUT  						"
      << endl;
  oss << "-ph				 	   true			 "
         "	"
      << endl;
  oss << "-temperature		   true				" << endl;
  oss << "-alkalinity   		   true				"
      << endl;
  // define what is written into output vector
  oss << "  -totals ";

  for (FLUID ::iterator it = fluid.begin(); it != fluid.end(); it++)
    if (it->first != "ph" && it->first != "temp" && it->first != "Alkalinity")
      oss << it->first << " ";
  oss << endl;
  oss << "END								"
         "		"
      << endl;

  if (iphreeqc.RunString(oss.str().c_str()) != 0)
    ERROR();

  ::VarInit(&vs);
  UpdateVector(vs, fluid);
  ::VarClear(&vs);
}

void PhreeqC::AdjustUnits(vector<std::pair<const char *, double>> &fluid) {
  // this converts initial units into mol/kgw which si the standart unit for
  // phreeqc
  VAR vs;
  std::ostringstream oss;
  oss << "Solution 1" << endl;
  oss << "units " << units.c_str() << endl;

  for (FLUID ::iterator it = fluid.begin(); it != fluid.end(); it++) {
    if (it->first == "pH")
      oss << it->first << " " << ConvertToString(it->second).c_str()
          << endl; //" charge" << endl;
    else if (it->first == "Alkalinity")
      oss << it->first << " " << ConvertToString(it->second).c_str()
          << " as HCO3" << endl;
    else
      oss << it->first << " " << ConvertToString(it->second).c_str() << endl;
  }
  oss << "SELECTED_OUTPUT  						"
      << endl;
  oss << "-ph				 	   true			 "
         "	"
      << endl;
  oss << "-temperature		   true				" << endl;
  oss << "-alkalinity   		   true				"
      << endl;
  // define what is written into output vector
  oss << "  -totals ";

  for (FLUID ::iterator it = fluid.begin(); it != fluid.end(); it++)
    if (it->first != "pH" && it->first != "temp" && it->first != "Alkalinity")
      oss << it->first << " ";
  oss << endl;
  oss << "END								"
         "		"
      << endl;

  if (iphreeqc.RunString(oss.str().c_str()) != 0)
    ERROR();

  ::VarInit(&vs);
  UpdateVector(vs, fluid);
  ::VarClear(&vs);
}

// concentr
void PhreeqC::GetConcentration(const char *type, double (&Conc)[dimX][dimY]) {
  boost::array<double, (dimX * dimY)> Conc_matrix;
  for (FLUID_array ::iterator it = Conc_matrices.begin();
       it != Conc_matrices.end(); it++) {
    if (it->first == type)
      Conc_matrix = it->second;
  }
  if (Conc_matrix.size() > 0) {
    for (int x = 0; x < dimX; x++)
      for (int y = 0; y < dimY; y++) {
        Conc[x][y] = Conc_matrix[x * dimY + y];
      }
  } else
    cout << "WARNING: could not find " << type << " in fluid array" << endl;
}

void PhreeqC::SetConcentration(const char *type, double Conc[dimX][dimY]) {
  for (FLUID_array ::iterator it = Conc_matrices.begin();
       it != Conc_matrices.end(); it++) {
    if (it->first == type) {
      for (int x = 0; x < dimX; x++)
        for (int y = 0; y < dimY; y++) {
          it->second[x * dimY + y] = Conc[x][y];
        }
    }
  }
}

double
PhreeqC::FindConcentration(const char *type,
                           vector<std::pair<const char *, double>> fluid) {
  double value = 0;
  for (FLUID ::iterator it = fluid.begin(); it != fluid.end(); it++) {
    if (boost::iequals(it->first, type))
      value = it->second;
  }
  return (value);
}

PhreeqC::FLUID PhreeqC::GetFluid(int x, int y, PhreeqC::FLUID_array fluid_c) {
  FLUID fluid;
  for (FLUID_array ::iterator it = fluid_c.begin(); it != fluid_c.end(); it++)
    fluid.push_back(make_pair(it->first, it->second[x * dimY + y]));
  return (fluid);
}

//==================================Main================================
double PhreeqC::dissCalcite(std::vector<double> PFluid, double radius,
                            double porosity, int NB, double &AdCalcite) {
  // not that concentration is phreeqc are given in mol/kgw
  // the calculatins performed are therefore for a volume of 1l
  // appropriate scaling factor have to be applyed after the simulation in order
  // to convert back to the respective scale of teh particles

  /*
  *#######
  #Calcite
  ########
  #
  # Plummer, L.N., Wigley, T.M.L., and Parkhurst, D.L., 1978,
  # American Journal of Science, v. 278, p. 179-216.
  #
  # Example of KINETICS data block for calcite rate:
  #
  #       KINETICS 1
  #       Calcite
  #	       -tol    1e-8
  #	       -m0     3.e-3
  #	       -m      3.e-3
  #	       -parms  5.0      0.6
  Calcite
    -start
     1 REM	Modified from Plummer and others, 1978
     2 REM	M = current moles of calcite
     3 REM	M0 = initial moles of calcite
     4 REM	parm(1) = Area/Volume, cm^2/L   (or cm^2 per cell)
     5 REM	parm(2) = exponent for M/M0 for surface area correction
     10  REM rate = 0 if no calcite and undersaturated
     20    si_cc = SI("Calcite")
     30    if (M <= 0  and si_cc < 0) then goto 300
     40    k1 = 10^(0.198 - 444.0 / TK )
     50    k2 = 10^(2.84 - 2177.0 / TK )
     60    if TC <= 25 then k3 = 10^(-5.86 - 317.0 / TK )
     70    if TC > 25 then k3 = 10^(-1.1 - 1737.0 / TK )
     80  REM surface area calculation
     90    t = 1
     100   if M0 > 0 then t = M/M0
     110   if t = 0 then t = 1
     120   area = PARM(1) * (t)^PARM(2)
     130   rf = k1 * ACT("H+") + k2 * ACT("CO2") + k3 * ACT("H2O")
     140 REM 1e-3 converts mmol to mol
     150   rate = area * 1e-3 * rf * (1 - 10^(2/3*si_cc))
     160   moles = rate * TIME
     170 REM do not dissolve more calcite than present
     180   if (moles > M) then moles = M
     190   if (moles >= 0) then goto 300
     200 REM do not precipitate more Ca or C(4) than present
     210   temp = TOT("Ca")
     220   mc  = TOT("C(4)")
     230   if mc < temp then temp = mc
     240   if -moles > temp then moles = -temp
     300 SAVE moles
    -end
  */
  iphreeqc.ClearAccumulatedLines();
  iphreeqc.SetOutputFileOn(false);
  iphreeqc.SetOutputFileName("CaDiss10");

  int I = 0, i = 0;
  float Vs = 36.92; // molar volume of calcite (cm^3 / mol)
  double initCa, modCa, Vp, Solid, J;
  double V_solid, V_fluid, r_fluid, A_fluid, AV;
  const char *Y;
  const char *Y2;
  const char *X;
  const char *X2;

  stringstream float2str1;
  stringstream float2str3;
  string name;
  string x, x2, y, y2;
  vector<string> Name;
  vector<string>::iterator N;

  // built element names string from vector----------------------------
  for (N = Name.begin(); N < Name.end(); N++) {
    string J = *N;
    if (I > 2)
      name += J;
    I++;
  }
  // calculate A/V---------------------------------------------------------
  // assuming sperical volume of particle

  V_solid = ((1 - porosity));
  V_fluid = (porosity);
  r_fluid = cbrt(((3.0 / 4.0) * V_fluid) / 3.14);
  A_fluid = 4.0 * 3.14 * (r_fluid * r_fluid);
  // AV = A_fluid / V_solid;

  // A_fluid = 1 - (3.14 * pow(radius * 200 ,2)) ;
  // A_fluid = (A_fluid / 500);
  // V_solid = pow( 500 * ((1 - radius)),3);

  // A_fluid = porosity * 25;
  // A_fluid = 4 * sqrt(A_fluid);
  // V_solid = (1-porosity) * 25;

  AV = 0.01 * (A_fluid / V_solid);

  float2str3 << AV;
  y = "-parms " + float2str3.str() + " 0.6";
  Y = y.c_str();

  // calculate moles of calcite present in particle------------------------
  V_solid /= Vs; // molV(calcite) = 36.92 cm^3/mol;
  float2str1 << V_solid;
  y2 = "-m " + float2str1.str();
  Y2 = y2.c_str();

  // create phreeq c simulation input======================================
  iphreeqc.AccumulateLine(
      "Solution 1_mod					"
      "		"); // defining solution in particle
  iphreeqc.AccumulateLine(
      "units mol/kgw							");
  for (F = PFluid.begin(); F < PFluid.end(); F++) {
    J = *F;
    stringstream float2str0(stringstream::in | stringstream::out);
    float2str0 << J;
    x = Name[i] + float2str0.str();
    X = x.c_str();
    iphreeqc.AccumulateLine(X);
    i++;
  }
  iphreeqc.AccumulateLine(
      "KINETICS 1							"
      "	 "); // define the reaction with:
  iphreeqc.AccumulateLine(
      "Calcite	    						 "); // name of
                                                                     // mineral
  iphreeqc.AccumulateLine("-tol 1e-8	    				"
                          "		 "); // define tolerance for interation
  iphreeqc.AccumulateLine(
      Y); // kinetic parameters (A/V and correction exponent (0.6))
  iphreeqc.AccumulateLine(
      "	-steps 3600 in 1		    		 "); // define the time
                                                             // :-)
  iphreeqc.AccumulateLine(Y2); // current moles of calcite in particle
  iphreeqc.AccumulateLine(
      "SELECTED_OUTPUT 						 ");
  iphreeqc.AccumulateLine("-simulation           false			 ");
  iphreeqc.AccumulateLine("-state                false			 ");
  iphreeqc.AccumulateLine("-solution             false			 ");
  iphreeqc.AccumulateLine("-distance             false			 ");
  iphreeqc.AccumulateLine("-time                 false			 ");
  iphreeqc.AccumulateLine("-step                 false			 ");
  iphreeqc.AccumulateLine(
      "-pe				   false			 ");
  iphreeqc.AccumulateLine(
      "-ph				   false			 ");
  iphreeqc.AccumulateLine(
      "-temperature		   false			 ");
  iphreeqc.AccumulateLine(
      "-alkalinity   		   false			 ");
  iphreeqc.AccumulateLine("-kinetic_reactants    Calcite		"
                          "	 ");

  // define what is written into output vector
  x2 = "  -totals " + name;
  X2 = x2.c_str();
  iphreeqc.AccumulateLine(X2);
  iphreeqc.AccumulateLine("END						"
                          "			 ");
  if (iphreeqc.RunAccumulated() != 0)
    ERROR();

  // access values in output vector++++++++++++++++++++++++++++++++++++++++
  double dCa, V, V_prev, rad = 0;
  int JJ;

  // calculate dissolved volume--------------------------------------------
  VAR v;
  ::VarInit(&v);

  for (int i = 0; i < iphreeqc.GetSelectedOutputRowCount(); ++i) {
    for (int j = 0; j < iphreeqc.GetSelectedOutputColumnCount(); ++j) {
      if (iphreeqc.GetSelectedOutputValue(i, j, &v) == VR_OK) {
        switch (v.type) {
        case TT_STRING:
          string S = v.sVal;
          if (S == "dk_Calcite")
            JJ = j;
          break;
        }
      }
    }
  }

  if (iphreeqc.GetSelectedOutputValue(1, JJ, &v) == VR_OK)
    ;
  else
    ERROR();
  switch (v.type) {
  case TT_DOUBLE:
    dCa = v.dVal;
    break;
  }
  ::VarClear(&v);
  // if (dCa != 0 ) cout<<dCa<<endl;
  // Needs scaling=========================================================
  double scaling_factor = 1;
  // calculate the new radius of the particle-----------------------------
  if (dCa != 0) {
    dCa *= scaling_factor;
    // Add amount of dissolved calcite to varibale--------------------------
    AdCalcite -= (dCa);
    // calculate new radius--------------------------------------------------

    // V = dCa * Vs;
    // V = dCa ;
    // V_prev = 4.0/3.0 * 3.14 * pow(radius,3);
    //	if (NB == 0) cout <<V_solid<<" "<<V<<endl;
    //	V_prev =  V_prev +  V;
    //	rad = (V_prev / ( 4.0/3.0 * 3.14));
    //	rad = cbrt(rad);
    rad = radius + dCa;
    if (rad < 0 || rad == 0)
      rad = radius;
    // if (NB == 0) cout<<radius<<" "<<rad<<endl;
    return (rad);
  } else
    return (radius);
}

// calculate SI=========================================================
double PhreeqC::SI(int time_step,
                   vector<std::pair<const char *, double>> c_fluid,
                   string phase, bool write_file) {
  VAR vs;
  double SI = 0;
  string si = "si_" + phase;
  string filename = si + boost::lexical_cast<std::string>(time_step);

  iphreeqc.SetOutputFileOn(write_file);
  iphreeqc.SetOutputFileName(filename.c_str());

  std::ostringstream oss;
  oss << "Solution 1_mod " << endl;
  oss << "units " << units.c_str() << endl;

  for (FLUID ::iterator it = c_fluid.begin(); it != c_fluid.end(); it++) {
    if (it->first == "Alkalinity")
      oss << it->first << " " << ConvertToString(it->second).c_str()
          << " as HCO3" << endl;
    else
      oss << it->first << " " << ConvertToString(it->second).c_str() << endl;
  }

  oss << "SELECTED_OUTPUT  	" << endl;
  oss << "-user_punch true" << endl;
  oss << "-si " << phase.c_str() << endl;
  oss << "END" << endl;

  if (iphreeqc.RunString(oss.str().c_str()) != 0)
    ERROR();

  double rho;
  ::VarInit(&vs);
  SI = ReturnValue(vs, si);
  ::VarClear(&vs);

  return (SI);
}

void PhreeqC::OUTPUT(int time_step,
                     vector<std::pair<const char *, double>> c_fluid,
                     vector<std::pair<std::string, double>> &local_result,
                     bool write_file) {
  VAR vs;
  double SI = 0;

  string filename = boost::lexical_cast<std::string>(time_step);

  iphreeqc.SetOutputFileOn(write_file);
  iphreeqc.SetOutputFileName(filename.c_str());

  std::ostringstream oss;
  oss << "Solution 1_mod " << endl;
  oss << "units " << units.c_str() << endl;

  for (FLUID ::iterator it = c_fluid.begin(); it != c_fluid.end(); it++) {
    if (it->first == "Alkalinity")
      oss << it->first << " " << ConvertToString(it->second).c_str()
          << " as HCO3" << endl;
    else
      oss << it->first << " " << ConvertToString(it->second).c_str() << endl;
  }

  oss << "SELECTED_OUTPUT  	" << endl;
  oss << "-user_punch true" << endl;
  for (vector<pair<std::string, double>>::const_iterator it =
           local_result.begin();
       it != local_result.end(); it++) {
    if (it->first != "density")
      oss << "-si " << it->first << endl;
  }
  oss << "USER_PUNCH" << endl;
  oss << "-head density" << endl;
  oss << "10 PUNCH RHO" << endl;
  oss << "END" << endl;

  if (iphreeqc.RunString(oss.str().c_str()) != 0)
    ERROR();

  ::VarInit(&vs);
  for (int i = 0; i < local_result.size(); i++) {
    if (local_result.at(i).first != "density") {
      string s = "si_" + local_result.at(i).first;
      local_result.at(i).second = ReturnValue(vs, s);
    } else
      local_result.at(i).second = ReturnValue(vs, local_result.at(i).first);
  }
  ::VarClear(&vs);
}

void PhreeqC::Equilibrate(vector<std::pair<const char *, double>> &fluid,
                          string phase) {
  string name = "equilibrium_" + phase;
  VAR vs;
  iphreeqc.SetOutputFileOn(true);
  iphreeqc.SetOutputFileName(name.c_str());
  std::ostringstream oss;
  oss << "Solution 1_mod " << endl;
  oss << "units " << units.c_str() << endl;

  for (FLUID ::iterator it = fluid.begin(); it != fluid.end(); it++) {
    if (it->first == "pH")
      oss << it->first << " " << ConvertToString(it->second).c_str() << endl;
    else if (it->first == "Alkalinity")
      oss << it->first << " " << ConvertToString(it->second).c_str()
          << " as HCO3" << endl;
    else
      oss << it->first << " " << ConvertToString(it->second).c_str() << endl;
  }

  oss << "EQUILIBRIUM_PHASES 1					" << endl;
  oss << phase.c_str() << " 0						"
      << endl;
  oss << "SELECTED_OUTPUT  						"
      << endl;
  oss << "-ph				 	   true			 "
         "	"
      << endl;
  oss << "-temperature		   true				" << endl;
  oss << "-alkalinity   		   true				"
      << endl;
  // define what is written into output vector
  oss << "  -totals ";

  for (FLUID ::iterator it = fluid.begin(); it != fluid.end(); it++)
    if (it->first != "ph" && it->first != "temp" && it->first != "Alkalinity")
      oss << it->first << " ";
  oss << endl;
  oss << "END								"
         "		"
      << endl;

  if (iphreeqc.RunString(oss.str().c_str()) != 0)
    ERROR();

  ::VarInit(&vs);
  UpdateVector(vs, fluid);
  ::VarClear(&vs);
}

/*
 * This is our first dummy reaction. The temperture dependence in dumped in k.
 *
 * r = k * (1 - SI) * A0/V * (m/m0)^n
 *
 * k: specific rate at given temperature [mol/m^2/s] (note we only specifcy 'a'
 * in 10^a) n: is function of grain size (n = 2/3 for "monodisperse population
 * of uniformly dissolving or growing spheres or cubes") A0: area [m^2] V:
 * volume [dm^3] m0: initial moles present [mol] time_step: number of time steps
 * tot_time: total time to react [s]
 * mineral: the mineral phase (e.g.Quartz)
 * formula: the formula of the mineral (e.g. SiO2)
 */
void PhreeqC::React(PhreeqC::FLUID fluid, double k, double n, double A0,
                    double V, double m0, double time_step, double tot_time,
                    string mineral, string formula, bool write_file) {
  VAR vs;
  iphreeqc.SetOutputFileOn(write_file);
  iphreeqc.SetOutputFileName("react_qz");

  std::ostringstream oss;

  // this is how we define a rate expression for phreeqc
  oss << "Rates" << endl;
  oss << "	R1			  " << endl;
  oss << "	-start " << endl;
  oss << "	10 A0 = parm(1)" << endl;
  oss << "	20 V = parm(2)" << endl;
  oss << "	30 rate = 10^" << ConvertToString(k).c_str() << " * (1-(SR"
      << mineral.c_str() << ")) * A0/V * (m/m0)^" << ConvertToString(n).c_str()
      << endl;
  oss << "	40 moles = rate * time" << endl;
  oss << "	50 save moles" << endl;
  oss << "	-end" << endl;

  // here we
  oss << "KINETICS 1 " << endl;
  oss << "R1 " << endl;
  oss << "	-formula " << formula.c_str() << endl;
  oss << "	-m0 " << ConvertToString(m0).c_str() << endl;
  oss << "	-parms " << ConvertToString(A0).c_str() << " "
      << ConvertToString(V).c_str() << endl;
  oss << "	-step " << ConvertToString(tot_time).c_str() << " in "
      << ConvertToString(time_step).c_str() << endl;
  oss << "	-runge_kutta 3 " << endl;
  oss << " 	-tol     1.e-8 " << endl;
  oss << "INCREMENTAL_REACTIONS true" << endl;

  oss << "Solution 1_mod " << endl;
  oss << "units " << units.c_str() << endl;

  for (FLUID ::iterator it = fluid.begin(); it != fluid.end(); it++)
    oss << it->first << " " << ConvertToString(it->second).c_str() << endl;

  oss << "SELECTED_OUTPUT  	" << endl;
  oss << "	-ph	true" << endl;
  oss << "	-totals Si" << endl;
  oss << "	-percent_error	true" << endl;
  oss << "	-kinetic_reactants R1" << endl;
  oss << "END" << endl;

  if (iphreeqc.RunString(oss.str().c_str()) != 0)
    ERROR();

  ::VarInit(&vs);
  //	cout <<  "total si in f: " << ReturnValue(vs, "Si(mol/kgw)") << endl;
  //	cout <<  "delta moles of SiO2: " << ReturnValue(vs, "dk_R1") << endl;
  ::VarClear(&vs);
}
