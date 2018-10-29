//This file contains a wrapper class for freestream-milne

#ifndef SRC_ITA_
#define SRC_ITA_

//#include "Parameter.h"
#include "Propagate.cpp"
#include "InitialConditions.cpp"
#include "LandauMatch.cpp"
#include "EquationOfState.cpp"
#include "HydroValidity.cpp"
#include "Memoryf.cpp"
#include "FileIO.cpp"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#define PI 3.141592654f
#define HBARC 0.197326938 //used to convert units of input / output hydro vectors

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
      //set default parameters in case of missing freestream_input file
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
      //read in chosen parameters from freestream_input if such a file exists
      readInParameters(params);
      //define some useful combinations
      params.DIM = params.DIM_X * params.DIM_Y;

      int DIM_X = params.DIM_X;
      int DIM_Y = params.DIM_Y;
      int DIM = params.DIM;
      int DIM_T = params.DIM_T;
      int DIM_PHIP = params.DIM_PHIP;
      float t0 = params.T0;
      float dt = params.DT;

      printf("Parameters are ...\n");
      printf("(DIM_X, DIM_Y, DIM_PHIP) = (%d, %d, %d)\n", params.DIM_X, params.DIM_Y, params.DIM_PHIP);
      printf("(DX, DY, DT) = (%.2f fm, %.2f fm, %.2f fm/c)\n", params.DX, params.DY, params.DT);
      printf("T0 = %.2f fm/c\n", params.T0);
      printf("DIM_T = %d \n", params.DIM_T);
      printf("E_FREEZE = %.3f GeV / fm^3 \n", params.E_FREEZE);
      printf("GAMMA = %.2f \n", params.GAMMA);
      if (params.EOS_TYPE == 1) printf("Using EoS : Conformal \n");
      else if (params.EOS_TYPE == 2) printf("Using EoS : Wuppertal-Budhapest \n");
      else if (params.EOS_TYPE == 3) printf("Using EoS : Lattice QCD + HRG matched.\n");

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
      float **hypertrigTable = NULL;
      hypertrigTable = calloc2dArrayf(hypertrigTable, 10, params.DIM_PHIP);

      //the density F()
      float **density = NULL;
      density = calloc2dArrayf(density, params.DIM, params.DIM_PHIP); // function of x,y,eta and rapidity

      //the previous value of the density F()
      float **density_p = NULL;
      density_p = calloc2dArrayf(density, params.DIM, params.DIM_PHIP);

      //initialize energy density
      //define a lower bound on energy density for all cells to regulate numerical noise in flow velocity in dilute regions
      float lower_tolerance = 1.0e-7;

      printf("setting initial conditions on energy density : ");
      if (params.IC_ENERGY == 1)
      {
        initializeEllipticalGauss(energyDensity, 5.0, 5.0, params);
        printf("Smooth Oblate Gaussian \n");
      }
      else if (params.IC_ENERGY == 2)
      {
        initializeEllipticalMCGauss(energyDensity, 5.0, 5.0,params);
        printf("Fluctuating Oblate Gaussian \n");
      }
      else if (params.IC_ENERGY == 3)
      {
        readDensityFile(energyDensity, "initial_profiles/e", params);
        printf("Reading from energy density file in initial_profiles/ \n");
      }
      else if (params.IC_ENERGY == 4)
      {
        readEnergyDensitySuperMCBlock(energyDensity, params);
        printf("Reading from superMC energy density file in initial_profiles/ \n");
      }
      else if (params.IC_ENERGY == 5)
      {
        //read in initial energy density using the initiliaze_from_vector() function
        //note that this is not safe - if one passes an empty vector it will not throw an error
        //converting units of energy density from GeV / fm^3 to fm^(-4)
        printf("Reading energy density from initial energy density vector\n");
        //do a value copy
        //try adding a small value everywhere to regulate problems with flow velocity in dilute regions
        for (int i = 0; i < params.DIM; i++) energyDensity[i] = init_energy_density[i] / (float)HBARC + lower_tolerance;
      }
      else
      {
        printf("Not a valid initial Condition... Goodbye\n");
        return 0;
      }

      //initialize the flow velocity
      initializeFlow(flowVelocity, params);

      //write initial energy density  to file
      writeScalarToFile(energyDensity, "initial_e", params);
      writeScalarToFileProjection(energyDensity, "initial_e_projection", params);

      //calculate total energy to check convergence
      float totalEnergy = 0.0;
      for (int is = 0; is < params.DIM; is++) totalEnergy += energyDensity[is];
      totalEnergy *= (params.DX * params.DY);
      printf("Total energy before evolution : %f \n", totalEnergy);

      //convert the energy density profile into the density profile F(t, x, y ; phip) to be propagated
      initializeDensity(energyDensity, density_p, params);

      //The main time step loop
      printf("Evolving T^munu via ITA Collision Dynamics \n");
      for (int it = 1; it < DIM_T + 1; it++)
      {
        float t = t0 + it * dt;
        printf("Step %d of %d : t = %.3f \n" , it, DIM_T, t);
        calculateHypertrigTable(hypertrigTable, params);
        //calculate the ten independent components of the stress tensor by integrating over rapidity and phi_p
        calculateStressTensor(stressTensor, density_p, hypertrigTable, params);
        //solve the eigenvalue problem for the energy density and flow velocity
        solveEigenSystem(stressTensor, energyDensity, flowVelocity, params);
        char e_file[255] = "";
        sprintf(e_file, "e_projection_%.3f", t);
        char ut_file[255] = "";
        sprintf(ut_file, "ut_projection_%.3f", t);
        writeScalarToFileProjection(energyDensity, e_file, params);
        writeVectorToFileProjection(flowVelocity, ut_file, 0, params);
        //propagate the density forward by one time step according to ITA EQN of Motion
        propagate(density, density_p, energyDensity, flowVelocity, params);

        //swap the density and previous value

        //this doesn't work?
        //std::swap(density, density_p);

        //value swap is clumsy

        for (int is = 0; is < DIM; is++)
        {
          for (int iphip = 0; iphip < DIM_PHIP; iphip++)
          {
            density_p[is][iphip] = density[is][iphip];
          }
        }

      }


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
      for (int is = 0; is < params.DIM; is++) totalEnergyAfter += energyDensity[is];
      totalEnergyAfter *= (params.DX * params.DY);
      printf("Total energy after streaming : %f \n", totalEnergyAfter);

      //check which fraction of total energy lies within freezeout surface, which lies in 'corona'
      float totalEnergyInsideHypersurf = 0.0;
      for (int is = 0; is < params.DIM; is++)
      {
        if ( (energyDensity[is] * HBARC) > params.E_FREEZE) totalEnergyInsideHypersurf += energyDensity[is];
      }
      totalEnergyInsideHypersurf *= (params.DX * params.DY);
      printf("Fraction of energy contained in Freezeout Hypersurface : %f \n", totalEnergyInsideHypersurf / totalEnergyAfter);

      //////////////////////////////////HYDRO VALIDITY//////////////////////////////////
      //bulk inv reynolds #
      float *R_Pi_Inv = NULL;
      R_Pi_Inv = (float *)calloc(params.DIM, sizeof(float));
      //shear inv reynolds #
      float *R_pimunu_Inv = NULL;
      R_pimunu_Inv = (float *)calloc(params.DIM, sizeof(float));
      calculateBulkInvReynolds(pressure, bulkPressure, R_Pi_Inv, params);
      calculateShearInvReynolds(energyDensity, pressure, shearTensor, R_pimunu_Inv, params);
      writeScalarToFileProjection(R_Pi_Inv, "R_Pi_Inv_projection", params);
      writeScalarToFileProjection(R_pimunu_Inv, "R_pimunu_Inv_projection", params);
      //int ctr = (DIM_Y * DIM_ETA * ((DIM_X - 1) / 2)) + (DIM_ETA * ((DIM_Y - 1) / 2)) + ((DIM_ETA - 1) / 2);
      //printf("R_Pi_Inv at center : %f \n", R_Pi_Inv[ctr]);
      //printf("R_pimunu_Inv at center : %f \n", R_pimunu_Inv[ctr]);
      free(R_Pi_Inv);
      free(R_pimunu_Inv);
      //////////////////////////////////HYDRO VALIDITY//////////////////////////////////


      printf("writing hydro variables\n");

      writeScalarToFile(energyDensity, "e", params);
      writeScalarToFile(pressure, "p", params);
      writeScalarToFile(bulkPressure, "bulk_PI", params);
      writeScalarToFileProjection(energyDensity, "e_projection", params);
      writeScalarToFileProjection(pressure, "p_projection", params);
      writeScalarToFileProjection(bulkPressure, "bulk_PI_projection", params);

      writeVectorToFile(flowVelocity, "u_tau", 0, params);
      writeVectorToFile(flowVelocity, "u_x", 1, params);
      writeVectorToFile(flowVelocity, "u_y", 2,params);
      writeVectorToFile(flowVelocity, "u_eta", 3,params);

      writeVectorToFileProjection(flowVelocity, "u_tau_projection", 0,params);
      writeVectorToFileProjection(flowVelocity, "u_x_projection", 1,params);
      writeVectorToFileProjection(flowVelocity, "u_y_projection", 2,params);
      writeVectorToFileProjection(flowVelocity, "u_eta_projection", 3,params);


      writeVectorToFile(shearTensor, "pi_tau_tau", 0,params);
      writeVectorToFile(shearTensor, "pi_tau_x", 1,params);
      writeVectorToFile(shearTensor, "pi_tau_y", 2,params);
      writeVectorToFile(shearTensor, "pi_tau_eta", 3,params);
      writeVectorToFile(shearTensor, "pi_x_x", 4,params);
      writeVectorToFile(shearTensor, "pi_x_y", 5,params);
      writeVectorToFile(shearTensor, "pi_x_eta", 6,params);
      writeVectorToFile(shearTensor, "pi_y_y", 7,params);
      writeVectorToFile(shearTensor, "pi_y_eta", 8,params);
      writeVectorToFile(shearTensor, "pi_eta_eta", 9,params);

      writeVectorToFileProjection(shearTensor, "pi_tau_tau_projection", 0,params);
      writeVectorToFileProjection(shearTensor, "pi_tau_x_projection", 1,params);
      writeVectorToFileProjection(shearTensor, "pi_tau_y_projection", 2,params);
      writeVectorToFileProjection(shearTensor, "pi_tau_eta_projection", 3,params);
      writeVectorToFileProjection(shearTensor, "pi_x_x_projection", 4,params);
      writeVectorToFileProjection(shearTensor, "pi_x_y_projection", 5,params);
      writeVectorToFileProjection(shearTensor, "pi_x_eta_projection", 6,params);
      writeVectorToFileProjection(shearTensor, "pi_y_y_projection", 7,params);
      writeVectorToFileProjection(shearTensor, "pi_y_eta_projection", 8,params);
      writeVectorToFileProjection(shearTensor, "pi_eta_eta_projection", 9,params);

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
          final_energy_density[is] = (double)energyDensity[is] * HBARC;
          final_pressure[is] = (double)pressure[is] * HBARC;
          final_ut[is] = (double)flowVelocity[0][is];
          final_ux[is] = (double)flowVelocity[1][is];
          final_uy[is] = (double)flowVelocity[2][is];
          final_un[is] = (double)flowVelocity[3][is];
          final_pitt[is] = (double)shearTensor[0][is] * HBARC;
          final_pitx[is] = (double)shearTensor[1][is] * HBARC;
          final_pity[is] = (double)shearTensor[2][is] * HBARC;
          final_pitn[is] = (double)shearTensor[3][is] * HBARC;
          final_pixx[is] = (double)shearTensor[4][is] * HBARC;
          final_pixy[is] = (double)shearTensor[5][is] * HBARC;
          final_pixn[is] = (double)shearTensor[6][is] * HBARC;
          final_piyy[is] = (double)shearTensor[7][is] * HBARC;
          final_piyn[is] = (double)shearTensor[8][is] * HBARC;
          final_pinn[is] = (double)shearTensor[9][is] * HBARC;
          final_Pi[is] = (double)bulkPressure[is] * HBARC;
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

    }

    #endif  // SRC_ITA_
