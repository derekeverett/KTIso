#ifndef SRC_ITA_
#define SRC_ITA_

#include "Propagate.cpp"
#include "InitialConditions.cpp"
#include "LandauMatch.cpp"
#include "EquationOfState.cpp"
#include "HydroValidity.cpp"
#include "Memoryf.cpp"
#include "FileIO.cpp"
//#include "FreezeOut.cpp"
//#include "cornelius-c++-1.3/cornelius.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>
//#include <algorithm>
#include <array>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

class ITA {
private:

public:
  ITA();
  ~ITA();

  int run_ita();

  // IS THIS VARIABLE NECESSARY
  int gridSize; //the total number of grid points in x, y, and eta : used for vector memory allocation

  //support to initilialize the energy density from a vector - useful for JETSCAPE
  //note units of argument should be GeV / fm^3
  //then we convert to fm^(-4)
  void initialize_from_vector(std::vector<float>);
  std::vector<float> init_energy_density;

  //support to write final hydro variables to vectors - useful for JETSCAPE
  //note we need to convert back to GeV / fm^3 units here
  void output_to_vectors(std::vector<double>&, //e
    std::vector<double>&, //p
    std::vector<double>&, //ut
    std::vector<double>&, //ux
    std::vector<double>&, //uy
    std::vector<double>&, //un
    std::vector<double>&, //pitt
    std::vector<double>&, //pitx
    std::vector<double>&, //pity
    std::vector<double>&, //pitn
    std::vector<double>&, //pixx
    std::vector<double>&, //pixy
    std::vector<double>&, //pixn
    std::vector<double>&, //piyy
    std::vector<double>&, //piyn
    std::vector<double>&, //pinn
    std::vector<double>&); //Pi

    std::vector<double> final_energy_density;
    std::vector<double> final_pressure;
    std::vector<double> final_ut;
    std::vector<double> final_ux;
    std::vector<double> final_uy;
    std::vector<double> final_un;
    std::vector<double> final_pitt;
    std::vector<double> final_pitx;
    std::vector<double> final_pity;
    std::vector<double> final_pitn;
    std::vector<double> final_pixx;
    std::vector<double> final_pixy;
    std::vector<double> final_pixn;
    std::vector<double> final_piyy;
    std::vector<double> final_piyn;
    std::vector<double> final_pinn;
    std::vector<double> final_Pi;

  };

  ITA::ITA() {

  }

  ITA::~ITA() {
  }

  //use this function to initialize energy density within JETSCAPE
  void ITA::initialize_from_vector(std::vector<float> energy_density_in) {
    init_energy_density = energy_density_in;
  }

  //use this function to return final hydro variables as vectors within JETSCAPE
  void ITA::output_to_vectors(std::vector<double> &energy_density_out,
    std::vector<double> &pressure_out,
    std::vector<double> &ut_out,
    std::vector<double> &ux_out,
    std::vector<double> &uy_out,
    std::vector<double> &un_out,
    std::vector<double> &pitt_out,
    std::vector<double> &pitx_out,
    std::vector<double> &pity_out,
    std::vector<double> &pitn_out,
    std::vector<double> &pixx_out,
    std::vector<double> &pixy_out,
    std::vector<double> &pixn_out,
    std::vector<double> &piyy_out,
    std::vector<double> &piyn_out,
    std::vector<double> &pinn_out,
    std::vector<double> &Pi_out) {
      energy_density_out = final_energy_density;
      pressure_out = final_pressure;
      ut_out = final_ut;
      ux_out = final_ux;
      uy_out = final_uy;
      un_out = final_un;
      pitt_out = final_pitt;
      pitx_out = final_pitx;
      pity_out = final_pity;
      pitn_out = final_pitn;
      pixx_out = final_pixx;
      pixy_out = final_pixy;
      pixn_out = final_pixn;
      piyy_out = final_piyy;
      piyn_out = final_piyn;
      pinn_out = final_pinn;
      Pi_out = final_Pi;
    }

    //where the magic happens
    int ITA::run_ita() {
      printf("Welcome to ITA\n");
      //declare parameter struct
      struct parameters params;
      //set default parameters in case of missing ita_input file
      params.OUTPUTFORMAT = 2;
      params.IC_ENERGY = 5;
      params.DIM_X = 101;
      params.DIM_Y = 101;
      params.DIM_PHIP = 200;
      params.DIM_T = 100;
      params.DX = 0.1;
      params.DY = 0.1;
      params.DT = 0.05;
      params.T0 = 0.1;
      params.EOS_TYPE = 1;
      params.E_FREEZE = 1.7;
      params.ALPHA = 10.0;
      params.COLLISIONS = 1;
      //read in chosen parameters from freestream_input if such a file exists
      readInParameters(params);
      //define some useful combinations
      params.DIM = params.DIM_X * params.DIM_Y;
      int DIM = params.DIM;
      int DIM_T = params.DIM_T;
      int DIM_VZ = params.DIM_VZ;
      int COLLISIONS = params.COLLISIONS;

      float t0 = params.T0;
      float dt = params.DT;
      int nt = params.DIM_T;
      float tf = t0 + (float)nt * dt;

      printf("Parameters are ...\n");
      printf("(DIM_X, DIM_Y, DIM_PHIP, DIM_VZ) = (%d, %d, %d, %d)\n", params.DIM_X, params.DIM_Y, params.DIM_PHIP, params.DIM_VZ);
      printf("(DX, DY, DT) = (%.2f fm, %.2f fm, %.2f fm/c)\n", params.DX, params.DY, params.DT);
      printf("T0 = %.2f fm/c\n", params.T0);
      printf("DIM_T = %d \n", params.DIM_T);
      printf("E_FREEZE = %.3f GeV / fm^3 \n", params.E_FREEZE);

      if (COLLISIONS)
      {
        printf("Evolving with collision term : ");
        printf("alpha = %.2f \n", params.ALPHA);
      }
      else printf("Evolving by freestreaming (no collisions) \n");


      if (params.EOS_TYPE == 1) printf("Using EoS : Conformal \n");
      else if (params.EOS_TYPE == 2) printf("Using EoS : Wuppertal-Budhapest \n");
      else { printf("Not a valid EoS! \n"); exit(-1); }

      const double hbarc = 0.197326938;

      //allocate and initialize memory
      printf("Allocating memory\n");
      //the ten independent components of the stress tensor
      float **stressTensor = NULL;
      stressTensor = calloc2dArrayf(stressTensor, 10, params.DIM);
      //the energy density
      float *energyDensity = NULL;
      energyDensity = (float *)calloc(params.DIM, sizeof(float));
      //the flow velocity
      float **flowVelocity = NULL;
      flowVelocity = calloc2dArrayf(flowVelocity, 4, params.DIM);
      //a table containing 10 rows for 10 independent combinations of p_(mu)p_(nu)
      float ***hypertrigTable = NULL;
      hypertrigTable = calloc3dArrayf(hypertrigTable, 10, params.DIM_PHIP, params.DIM_VZ);

      //the table containing roots and weights for quadrature in v_z
      float **vz_quad = NULL;
      vz_quad = calloc2dArrayf(vz_quad, params.DIM_VZ, 2);

      //read in table from file
      /*
      ifstream vz_quad_file;
      vz_quad_file.open("ita_tables/GaussHermite20pts.csv");
      if (vz_quad_file)
      {
        float root = 0.0;
        float weight = 0.0;
        for (int row = 0; row < params.DIM_VZ; row++)
        {
          vz_quad_file >> root >> weight;
          vz_quad[row][0] = root;
          vz_quad[row][1] = weight;
        }
      }
      else printf("No quadrature table found!");
      vz_quad_file.close();
      */

      //construct table on the fly, using cubic spacing in v_z from 0 to 1
      //note that the weights stored in this table are for integration over v_z
      //and should not be used as the differential delta v_z when performing the MacCormack step
      for (int ivz = 0; ivz < DIM_VZ; ivz++)
      {
        float dvz_linear = 2.0 / (2*DIM_VZ - 1.0);
        float vz_linear = 0.0 + (dvz_linear / 2.0) + (ivz * dvz_linear);
        float vz_cub = pow(vz_linear, 3.0);
        vz_quad[ivz][0] = vz_cub;
      }

      for (int ivz = 0; ivz < DIM_VZ - 1; ivz++)
      {
        float dvz_cub = vz_quad[ivz+1][0] - vz_quad[ivz][0];
        vz_quad[ivz][1] = dvz_cub;
      }

      //the density F()
      float ***density = NULL;
      density = calloc3dArrayf(density, params.DIM, params.DIM_PHIP, params.DIM_VZ); // function of x,y,eta and rapidity

      //the table containing the v_z spacing, for adaptive integration routine
      //float ***dvz_table = NULL;
      //dvz_table = calloc3dArrayf(dvz_table, params.DIM, params.DIM_PHIP, params.DIM_VZ);

      //the previous value of the density F()
      float ***density_p = NULL;
      density_p = calloc3dArrayf(density_p, params.DIM, params.DIM_PHIP, params.DIM_VZ);

      //the integral of  1 / (udotv)^2 / (solid angle); this should be one if energy matching is satisfied numerically
      float *energyMatchIntegral = NULL;
      energyMatchIntegral = (float *)calloc(params.DIM, sizeof(float));

      // u^2 must be 1 numerically
      float *unorm = NULL;
      unorm = (float *)calloc(params.DIM, sizeof(float));

      //initialize energy density
      initializeEnergyDensity(energyDensity, init_energy_density, params);
      //initialize the flow velocity
      initializeFlow(flowVelocity, params);
      //write initial energy density  to file
      writeScalarToFile(energyDensity, (char *)"initial_e", params);
      writeScalarToFileProjection(energyDensity, (char *)"initial_e_projection", params);

      //convert the energy density profile into the density profile F(t, x, y ; phip, xi) to be propagated
      initializeDensity(energyDensity, density_p, vz_quad, params);
      //calculate entries in trig table - time independent in cartesian case
      calculateHypertrigTable(hypertrigTable, vz_quad, params);

      //initialize the table of dvz
      //initializeDVZTable(dvz_table, params);

      //calculate total energy to check convergence
      calculateStressTensor(stressTensor, density_p, hypertrigTable, vz_quad, t0, params);
      float totalEnergy = 0.0;
      for (int is = 0; is < params.DIM; is++) totalEnergy += stressTensor[0][is];
      totalEnergy *= (params.DX * params.DY * t0);
      printf("Total energy before evolution : %f \n", totalEnergy);

      //useful for plotting the momentum dependence of distribution function
      float **F_vz_phip = NULL;
      F_vz_phip = calloc2dArrayf(F_vz_phip, params.DIM_VZ, params.DIM_PHIP);

      //The main time step loop
      printf("Evolving F(x,y;phip,vz) via ITA Eqns of Motion \n");

      //FREEZEOUT
      //initialize cornelius for freezeout surface finding
      /*
      int dim;
      float *lattice_spacing;
      dim = 3;
      lattice_spacing = new float[dim];
      lattice_spacing[0] = dt;
      lattice_spacing[1] = dx;
      lattice_spacing[2] = dy;

      float ****energy_density_evoution;
      energy_density_evoution = calloc4dArrayf(energy_density_evoution, 2, DIM_X, DIM_Y, 1);
      //make an array to store all the hydrodynamic variables
      //to be written to file once the freezeout surface is determined by the critical energy density
      int n_hydro_vars = 10; //u1, u2, u3, e, pi11, pi12, pi13, pi22, pi23, Pi (the temperature and pressure are calclated with EoS)
      float *****hydrodynamic_evoution;
      hydrodynamic_evoution = calloc5dArrayf(hydrodynamic_evoution, n_hydro_vars, 2, DIM_X, DIM_Y, 1);
      //for 2+1D simulations
      float ***hyperCube3D;
      hyperCube3D = calloc3dArrayf(hyperCube3D, 2, 2, 2);
      //open the freezeout surface file
      ofstream freezeoutSurfaceFile;
      freezeoutSurfaceFile.open("output/surface.dat");
      */
      //FREEZEOUT

      //write to file every write_freq steps
      int write_freq = 1;

      //MAIN TIME STEP LOOP
      for (int it = 0; it < DIM_T + 1; it++)
      {
        float t = t0 + it * dt;

        //the index in grid center
        int icenter = DIM / 2;

        if (it % 1 == 0) printf("Step %d of %d : t = %.3f : e = %.3f GeV/fm^3 \n" , it, DIM_T, t, energyDensity[icenter] * hbarc);

        //get momentum dependence at center of grid
        std::ofstream myfile;
        char filename[255] = "";
        sprintf(filename, "output/F_vz_phip_%.3f.dat", t);
        myfile.open(filename);
        for (int ivz = 0; ivz < params.DIM_VZ; ivz++)
        {
          for (int iphip = 0; iphip < params.DIM_PHIP; iphip++)
          {
            myfile << density_p[icenter][iphip][ivz] << " ";
          }
          myfile << "\n";
        }
        myfile.close();

        //calculate the ten independent components of the stress tensor by integrating over phi_p and vz
        calculateStressTensor(stressTensor, density_p, hypertrigTable, vz_quad, t, params);
        //solve the eigenvalue problem for the energy density and flow velocity
        solveEigenSystem(stressTensor, energyDensity, flowVelocity, params);

        //solve for shear stress also
        //then use CORNELIUS to construct a freezeout surface

        if (it % write_freq == 0)
        {
          char e_file[255] = "";
          char T00_file[255] = "";
          char ux_file[255] = "";
	        char un_file[255] = "";
          sprintf(e_file, "e_projection_%.3f", t);
          sprintf(T00_file, "T00_projection_%.3f", t);
          sprintf(ux_file, "ux_projection_%.3f", t);
	        sprintf(un_file, "un_projection_%.3f", t);
          writeScalarToFileProjection(energyDensity, e_file, params);
          writeVectorToFileProjection(stressTensor, T00_file, 0, params);
          writeVectorToFileProjection(flowVelocity, ux_file, 1, params);
	        writeVectorToFileProjection(flowVelocity, un_file, 3, params);
        }

        propagate(density, density_p, energyDensity, flowVelocity, vz_quad, t, params);

        float totalEnergy = 0.0;
        for (int is = 0; is < params.DIM; is++) totalEnergy += stressTensor[0][is];
        totalEnergy *= (params.DX * params.DY * t);
        printf("Total energy : %f \n", totalEnergy);

      } // for (int it = 0; it < DIM_T; it++)

      //variables to store the hydrodynamic variables after the Landau matching is performed
      //the pressure
      float *pressure = NULL;
      pressure = (float *)calloc(params.DIM, sizeof(float));
      //the bulk pressure Pi
      float *bulkPressure = NULL;
      bulkPressure = (float *)calloc(params.DIM, sizeof(float));
      //the shear stress tensor
      float **shearTensor = NULL;
      shearTensor = calloc2dArrayf(shearTensor, 10, params.DIM); //calculate 10 components, can check tracelessness/orthogonality for accuracy
      //calculate hydro variables
      calculatePressure(energyDensity, pressure, params);
      calculateBulkPressure(stressTensor, energyDensity, pressure, bulkPressure, params);
      calculateShearViscTensor(stressTensor, energyDensity, flowVelocity, pressure, bulkPressure, shearTensor, params);

      float totalEnergyAfter = 0.0;
      for (int is = 0; is < params.DIM; is++) totalEnergyAfter += stressTensor[0][is];
      totalEnergyAfter *= (params.DX * params.DY * tf);
      printf("Total energy after evolution : %f \n", totalEnergyAfter);

      //check which fraction of total energy lies within freezeout surface, which lies in 'corona'
      float totalEnergyAfterLRF = 0.0;
      for (int is = 0; is < params.DIM; is++) totalEnergyAfterLRF += energyDensity[is];
      totalEnergyAfterLRF *= (params.DX * params.DY);

      float totalEnergyInsideHypersurf = 0.0;
      for (int is = 0; is < params.DIM; is++)
      {
        if ( (energyDensity[is] * hbarc) > params.E_FREEZE) totalEnergyInsideHypersurf += energyDensity[is];
      }
      totalEnergyInsideHypersurf *= (params.DX * params.DY);
      printf("Fraction of LRF energy contained in Freezeout Hypersurface : %f \n", totalEnergyInsideHypersurf / totalEnergyAfterLRF);

      //////////////////////////////////HYDRO VALIDITY//////////////////////////////////
      //bulk inv reynolds #
      float *R_Pi_Inv = NULL;
      R_Pi_Inv = (float *)calloc(params.DIM, sizeof(float));
      //shear inv reynolds #
      float *R_pimunu_Inv = NULL;
      R_pimunu_Inv = (float *)calloc(params.DIM, sizeof(float));
      calculateBulkInvReynolds(pressure, bulkPressure, R_Pi_Inv, params);
      calculateShearInvReynolds(energyDensity, pressure, shearTensor, R_pimunu_Inv, params);
      writeScalarToFileProjection(R_Pi_Inv, (char *)"R_Pi_Inv_projection", params);
      writeScalarToFileProjection(R_pimunu_Inv, (char *)"R_pimunu_Inv_projection", params);
      //int ctr = (DIM_Y * DIM_ETA * ((DIM_X - 1) / 2)) + (DIM_ETA * ((DIM_Y - 1) / 2)) + ((DIM_ETA - 1) / 2);
      //printf("R_Pi_Inv at center : %f \n", R_Pi_Inv[ctr]);
      //printf("R_pimunu_Inv at center : %f \n", R_pimunu_Inv[ctr]);
      free(R_Pi_Inv);
      free(R_pimunu_Inv);
      //////////////////////////////////HYDRO VALIDITY//////////////////////////////////


      printf("writing hydro variables\n");

      writeScalarToFile(energyDensity, (char *)"e", params);
      writeScalarToFile(pressure, (char *)"p", params);
      writeScalarToFile(bulkPressure, (char *)"bulk_PI", params);
      writeScalarToFileProjection(energyDensity, (char *)"e_projection", params);
      writeScalarToFileProjection(pressure, (char *)"p_projection", params);
      writeScalarToFileProjection(bulkPressure, (char *)"bulk_PI_projection", params);

      writeVectorToFile(flowVelocity, (char *)"u_tau", 0, params);
      writeVectorToFile(flowVelocity, (char *)"u_x", 1, params);
      writeVectorToFile(flowVelocity, (char *)"u_y", 2,params);
      writeVectorToFile(flowVelocity, (char *)"u_eta", 3,params);

      writeVectorToFileProjection(flowVelocity, (char *)"u_tau_projection", 0,params);
      writeVectorToFileProjection(flowVelocity, (char *)"u_x_projection", 1,params);
      writeVectorToFileProjection(flowVelocity, (char *)"u_y_projection", 2,params);
      writeVectorToFileProjection(flowVelocity, (char *)"u_eta_projection", 3,params);

      writeVectorToFile(shearTensor, (char *)"pi_tau_tau", 0,params);
      writeVectorToFile(shearTensor, (char *)"pi_tau_x", 1,params);
      writeVectorToFile(shearTensor, (char *)"pi_tau_y", 2,params);
      writeVectorToFile(shearTensor, (char *)"pi_tau_eta", 3,params);
      writeVectorToFile(shearTensor, (char *)"pi_x_x", 4,params);
      writeVectorToFile(shearTensor, (char *)"pi_x_y", 5,params);
      writeVectorToFile(shearTensor, (char *)"pi_x_eta", 6,params);
      writeVectorToFile(shearTensor, (char *)"pi_y_y", 7,params);
      writeVectorToFile(shearTensor, (char *)"pi_y_eta", 8,params);
      writeVectorToFile(shearTensor, (char *)"pi_eta_eta", 9,params);

      writeVectorToFileProjection(shearTensor, (char *)"pi_tau_tau_projection", 0,params);
      writeVectorToFileProjection(shearTensor, (char *)"pi_tau_x_projection", 1,params);
      writeVectorToFileProjection(shearTensor, (char *)"pi_tau_y_projection", 2,params);
      writeVectorToFileProjection(shearTensor, (char *)"pi_tau_eta_projection", 3,params);
      writeVectorToFileProjection(shearTensor, (char *)"pi_x_x_projection", 4,params);
      writeVectorToFileProjection(shearTensor, (char *)"pi_x_y_projection", 5,params);
      writeVectorToFileProjection(shearTensor, (char *)"pi_x_eta_projection", 6,params);
      writeVectorToFileProjection(shearTensor, (char *)"pi_y_y_projection", 7,params);
      writeVectorToFileProjection(shearTensor, (char *)"pi_y_eta_projection", 8,params);
      writeVectorToFileProjection(shearTensor, (char *)"pi_eta_eta_projection", 9,params);

      //support for JETSCAPE - write hydro variables to vectors
      final_energy_density.resize(params.DIM);
      final_pressure.resize(params.DIM);
      final_ut.resize(params.DIM);
      final_ux.resize(params.DIM);
      final_uy.resize(params.DIM);
      final_un.resize(params.DIM);
      final_pitt.resize(params.DIM);
      final_pitx.resize(params.DIM);
      final_pity.resize(params.DIM);
      final_pitn.resize(params.DIM);
      final_pixx.resize(params.DIM);
      final_pixy.resize(params.DIM);
      final_pixn.resize(params.DIM);
      final_piyy.resize(params.DIM);
      final_piyn.resize(params.DIM);
      final_pinn.resize(params.DIM);
      final_Pi.resize(params.DIM);

      if ( (params.OUTPUTFORMAT == 2) || (params.OUTPUTFORMAT == 3) )
      {
        for (int is = 0; is < params.DIM; is++)
        {
          //converting back to GeV / fm^3 for use in JETSCAPE
          final_energy_density[is] = (double)energyDensity[is] * hbarc;
          final_pressure[is] = (double)pressure[is] * hbarc;
          final_ut[is] = (double)flowVelocity[0][is];
          final_ux[is] = (double)flowVelocity[1][is];
          final_uy[is] = (double)flowVelocity[2][is];
          final_un[is] = (double)flowVelocity[3][is];
          final_pitt[is] = (double)shearTensor[0][is] * hbarc;
          final_pitx[is] = (double)shearTensor[1][is] * hbarc;
          final_pity[is] = (double)shearTensor[2][is] * hbarc;
          final_pitn[is] = (double)shearTensor[3][is] * hbarc;
          final_pixx[is] = (double)shearTensor[4][is] * hbarc;
          final_pixy[is] = (double)shearTensor[5][is] * hbarc;
          final_pixn[is] = (double)shearTensor[6][is] * hbarc;
          final_piyy[is] = (double)shearTensor[7][is] * hbarc;
          final_piyn[is] = (double)shearTensor[8][is] * hbarc;
          final_pinn[is] = (double)shearTensor[9][is] * hbarc;
          final_Pi[is] = (double)bulkPressure[is] * hbarc;
        }
      }

      //free the memory
      free2dArrayf(stressTensor, 10);
      free(energyDensity);
      free2dArrayf(flowVelocity, 4);
      free(pressure);
      free(bulkPressure);
      free2dArrayf(shearTensor, 10);

      printf("Done... Goodbye!\n");
      return 1;

    }

    #endif  // SRC_ITA_
