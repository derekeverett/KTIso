#pragma once
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include "Parameter.h"
#define THETA_FUNCTION(X) ((float)X < (float)0 ? (float)0 : (float)1)

//this creates the initial F function, a function of spatial coordinates and momentum phip and velocity vz
//this function assumes F is initially isotropic in phi_p
void initializeDensity(float *energyDensity, float ***density, float **vz_quad, parameters params)
{
  int ntot = params.ntot;
  int nphip = params.nphip;
  int nvz = params.nvz;
  float dvz_2 = 2.0 / (float)(nvz - 1);

  float b = 0.3; //width of distribution in v_z
  float norm = sqrt(2.0 * M_PI) * b * erf( 1.0 / ( sqrt(2.0) * b) ); //normalization of gaussian on (-1, 1)

  if (nvz > 1)
  {

    //check that v_z distn is normalized
    /*
    float check_norm = 0.0;
    for (int ivz = 0; ivz < nvz; ivz++)
    {
      float vz = -1.0 + (float)ivz * dvz_2;
      //float vz = vz_quad[ivz][0];
      //float dvz = vz_quad[ivz][1];
      check_norm += exp(-vz * vz / (2.0 * b * b)) / norm * dvz_2;
      //check_norm += exp(-vz * vz / (2.0 * b * b)) / norm * dvz;
    }
    */

    //even function, integrated on even domain
    //check_norm *= 2.0;
    //std::cout << "Checking norm of v_z distribution : norm = " << check_norm << std::endl;

    #pragma omp parallel for
    for (int is = 0; is < ntot; is++)
    {
      float e0 = energyDensity[is];
      for (int iphip = 0; iphip < nphip; iphip++)
      {
        //try intializing dependence on v_z by a flat distribution
        //alternatively try gaussian or delta function ?
        for (int ivz = 0; ivz < nvz; ivz++)
        {
          //density[is][iphip][ivz] = val;
          float vz = -1.0 + (float)ivz * dvz_2;
          //float vz = vz_quad[ivz][0];

          //density[is][iphip][ivz] = e0 * exp(-vz * vz / 0.1); //this distribution is not normalized

          //density[is][iphip][ivz] = e0 * exp(-vz * vz / (2.0 * b * b)) / norm; //gaussian in v_z (P_l/P_T != 1)

          density[is][iphip][ivz] = e0 / 2.; //flat in v_z => PL/P_T = 1 Initially Isotropic Distribution
        }
      } //for (int iphip = 0; iphip < nphip; iphip++)
    } //for (int is = 0; is < ntot; is++)
  } //if (nvz > 1)

  else
  {
    #pragma omp parallel for
    for (int is = 0; is < ntot; is++)
    {
      float e0 = energyDensity[is];
      for (int iphip = 0; iphip < nphip; iphip++)
      {
        //density[is][iphip][0] = e0;
        density[is][iphip][0] = 2.0 * e0;
      } //for (int iphip = 0; iphip < nphip; iphip++)
    } //for (int is = 0; is < ntot; is++)
  }
}


//this function assumes the initial energy density has already been read, and
//initializes the phi_p and v_2 distribution following the files provided in
// random_fields/initial_psi_profiles and random_fields/initial_v2_profiles
//this function is only written for boost invariant mode !
void initializeDensity_color_domains(float *energyDensity, float ***density, float **vz_quad, parameters params)
{
  int ntot = params.ntot;
  int nphip = params.nphip;
  float dphip = (2.0 * M_PI) / float(params.nphip);

  //get the momentum modulation phase velocity
  //float w_D = params.w_D;

  float dx = params.dx;
  float dy = params.dy;
  int nx = params.nx;
  int ny = params.ny;


  //Print the random field parameters used
  std::fstream newfile;
  newfile.open("random_fields/random_field_param.txt",std::ios::in);
  std::cout << "Using initial momentum anisotropic fields with the parameters :" << std::endl;
  if (newfile.is_open()){   //checking whether the file is open
     std::string tp;
     while(getline(newfile, tp)){ //read data from file object and put it into string.
        std::cout << tp << "  "; //print the data of the string
     }
  std::cout<<'\n';
   }
  newfile.close(); //close the file object.
  //now loop over all cells in each patch
  std::ifstream blockFile1;
  blockFile1.open("random_fields/initial_psi_profiles/psi.dat");
  std::ifstream blockFile2;
  blockFile2.open("random_fields/initial_v2_profiles/v2.dat");
  float psi_p = 0.;
  float v2 = 0.;
  float integral_temp = 0.0;
  if (blockFile1.is_open())
  {
  if (blockFile2.is_open())
  {
  for (int ix = 0; ix < nx; ix++)
    {
    for (int iy = 0; iy < ny; iy++)
        {
        int is = (ny) * ix + iy; //the column packed index spanning x, y
        float e0 = energyDensity[is]; //get the value of the initial energy density

        blockFile1 >> psi_p; //get the value of the color domain phase angle psi from file
        blockFile2 >> v2; //get the value of the color domain eliptic flow magnitude from file
        //std::cout << "v2 value loaded :" << v2 <<std::endl;

        for (int iphip = 0; iphip < nphip; iphip++)
            {
            //float psi_p = (M_PI/2); //This is for testing. We will remove this soon
            float phi_p = iphip * dphip;
            if (v2 >= 0.5)
            {
              v2 = 0.49999; // put an anchor at 0.5 to avoid unphysical negative energies.
            }
            float fourier_sum = (1. + 2. * v2 * cos(2. * (phi_p-psi_p))); //No need to include a 2pi factor in the denominator.
            //float fourier_sum = 1. + 2. * 0 * cos(2. * (phi_p-0.));
            density[is][iphip][0] = 2.0 * e0 * fourier_sum; // multiply the initial energy density by the phi_p angular fct., remember factor of two for boost-invar F(v_z) ~ delta(v_z)
            //integral_temp +=  density[is][iphip][0];
            //integral_temp += (dphip / (4.0 * M_PI)) * (density[is][iphip][0])*(params.dx * params.dy * params.t0);

            } //for (int iphip = 0; iphip < nphip; iphip++)
        }
    }
}
else
{
  printf("Could not find random v2 profile in random_fields!");
}
blockFile2.close();
}
else
{
  printf("Could not find random psi profile in random_fields!");
}
blockFile1.close();
//std::cout << "Total integrated energy here is  " << integral_temp << std::endl;
}

//this creates the initial F function, a function of spatial coordinates and momentum phip and velocity vz
// assumes F is initially an anisotropic function of phi_p
void initializeDensityAniso(float *energyDensity, float ***density, float **vz_quad, parameters params)
{
  int ntot = params.ntot;
  int nphip = params.nphip;
  int nvz = params.nvz;
  float dvz_2 = 2.0 / (float)(nvz - 1);

  float b = 0.3; //width of distribution in v_z
  float norm = sqrt(2.0 * M_PI) * b * erf( 1.0 / ( sqrt(2.0) * b) ); //normalization of gaussian on (-1, 1)

  if (nvz > 1)
  {
    //check that v_z distn is normalized
    float check_norm = 0.0;
    for (int ivz = 0; ivz < nvz; ivz++)
    {
      float vz = -1.0 + (float)ivz * dvz_2;
      //float vz = vz_quad[ivz][0];
      //float dvz = vz_quad[ivz][1];
      check_norm += exp(-vz * vz / (2.0 * b * b)) / norm * dvz_2;
      //check_norm += exp(-vz * vz / (2.0 * b * b)) / norm * dvz;

    }

    //even function, integrated on even domain
    //check_norm *= 2.0;
    std::cout << "Checking norm of v_z distribution : norm = " << check_norm << std::endl;

    #pragma omp parallel for
    for (int is = 0; is < ntot; is++)
    {
      float e0 = energyDensity[is];
      for (int iphip = 0; iphip < nphip; iphip++)
      {
        //try intializing dependence on v_z by a flat distribution
        //alternatively try gaussian or delta function ?
        for (int ivz = 0; ivz < nvz; ivz++)
        {
          //density[is][iphip][ivz] = val;
          float vz = -1.0 + (float)ivz * dvz_2;
          //float vz = vz_quad[ivz][0];
          //density[is][iphip][ivz] = e0 * exp(-vz * vz / 0.1); //this distribution is not normalized
          density[is][iphip][ivz] = e0 * exp(-vz * vz / (2.0 * b * b)) / norm;
        }
      } //for (int iphip = 0; iphip < nphip; iphip++)
    } //for (int is = 0; is < ntot; is++)
  } //if (nvz > 1)

  else
  {
    #pragma omp parallel for
    for (int is = 0; is < ntot; is++)
    {
      float e0 = energyDensity[is];
      for (int iphip = 0; iphip < nphip; iphip++)
      {
        float phip = float(iphip) * (2.0 * M_PI) / float(nphip);

        //density[is][iphip][0] = e0;
        density[is][iphip][0] = 4.0 * e0 * cos(phip) * cos(phip);
      } //for (int iphip = 0; iphip < nphip; iphip++)
    } //for (int is = 0; is < ntot; is++)
  }
}

void initializeZero(float *density, parameters params)
{
  int ntot = params.ntot;
  for (int is = 0; is < ntot; is++)
  {
    density[is] = 0.0;
  }
}

void initializeFlowZero(float **flowVelocity, parameters params)
{
  int ntot = params.ntot;
  for (int is = 0; is < ntot; is++)
  {
    flowVelocity[0][is] = 1.0;
    flowVelocity[1][is] = 0.0;
    flowVelocity[2][is] = 0.0;
    flowVelocity[3][is] = 0.0;
  }
}

void readFlowVelocityBlock(float **flowVelocity, parameters params)
{
  int nx = params.nx;
  int ny = params.ny;
  int ntot = params.ntot;

  //first read in the transverse flow profile
  float temp = 0.0;
  std::ifstream blockFile;
  blockFile.open("initial_profiles/ux.dat");
  if (blockFile.is_open())
  {
    for (int iy = 0; iy < ny; iy++)
    {
      for (int ix = 0; ix < nx; ix++)
      {
        blockFile >> temp;
        int is = (ny) * ix + iy; //the column packed index spanning x, y
        flowVelocity[1][is] = temp;
      }
    }
  }

  else
  {
    printf("Could not find initial profile in initial_profiles!");
  }
  blockFile.close();

  blockFile.open("initial_profiles/uy.dat");
  if (blockFile.is_open())
  {
    for (int iy = 0; iy < ny; iy++)
    {
      for (int ix = 0; ix < nx; ix++)
      {
        blockFile >> temp;
        int is = (ny) * ix + iy; //the column packed index spanning x, y
        flowVelocity[2][is] = temp;
      }
    }
  }

  else
  {
    printf("Could not find initial profile in initial_profiles!");
  }
  blockFile.close();

  //now set u^eta to zero for boost invariance and u^t = 1 - (u^x*u^x + u^y*u^y)
  for (int is = 0; is < ntot; is++)
  {
    float ux = flowVelocity[1][is];
    float uy = flowVelocity[2][is];
    flowVelocity[0][is] = 1.0 - (ux*ux + uy*uy);
    flowVelocity[3][is] = 0.0;
  }

}


void initializeGauss(float *density, float b, parameters params) // b is the variance ('spherically' symmetric)
{
  int ntot = params.ntot;
  int nx = params.nx;
  int ny = params.ny;
  float dx = params.dx;
  float dy = params.dy;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (int is = 0; is < ntot; is++)
  {
    int ix = is / (ny );
    int iy = (is - (ny * ix));

    //does it work for even number of points?
    float x = (float)ix * dx  - ((float)(nx-1)) / 2.0 * dx;
    float y = (float)iy * dy  - ((float)(ny-1)) / 2.0 * dy;

    density[is] = e0 * exp( -(1.0 / b) * ( (x * x) + (y * y) ) );
  }
}

void initializeEllipticalGauss(float *density, float bx, float by, parameters params) // bx is the x variance etc...
{
  float hbarc = params.hbarc;
  float regulate = 1.0e-20;
  int ntot = params.ntot;
  int nx = params.nx;
  int ny = params.ny;
  float dx = params.dx;
  float dy = params.dy;

  float T0_GeV = 0.8; //initial tempeature in GeV
  float T0 = T0_GeV / hbarc; //initial temperature in fm^-1
  float e0 = energyDensityFromTemperature(T0);

  for (int is = 0; is < ntot; is++)
  {
    int ix = is / (ny);
    int iy = (is - (ny * ix));

    //does it work for even number of points?
    float x = (float)ix * dx  - ((float)(nx-1)) / 2.0 * dx;
    float y = (float)iy * dy  - ((float)(ny-1)) / 2.0 * dy;

    density[is] = e0 * exp(-(1.0 / (2.0 * bx * bx)) * (x * x)) * exp(-(1.0 / (2.0 * by * by)) * (y * y)) + regulate;
  }
}

void initializeMCGauss(float * density, float b, parameters params)
{
  int ntot = params.ntot;
  int nx = params.nx;
  int ny = params.ny;
  float dx = params.dx;
  float dy = params.dy;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (int is = 0; is < ntot; is++)
  {
    int ix = is / (ny);
    int iy = (is - (ny * ix));

    //does it work for even number of points?
    float x = (float)ix * dx  - ((float)(nx-1)) / 2.0 * dx;
    float y = (float)iy * dy  - ((float)(ny-1)) / 2.0 * dy;

    density[is] = e0 * ( (float)rand() / RAND_MAX ) * exp( -(1.0 / b) * ((x * x) + (y * y) ) );
  }
}

void initializeEllipticalMCGauss(float *density, float bx, float by, parameters params) // bx is the x variance etc...
{
  int ntot = params.ntot;
  int nx = params.nx;
  int ny = params.ny;
  float dx = params.dx;
  float dy = params.dy;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (int is = 0; is < ntot; is++)
  {
    int ix = is / (ny);
    int iy = (is - (ny * ix));

    //does it work for even number of points?
    float x = (float)ix * dx  - ((float)(nx-1)) / 2.0 * dx;
    float y = (float)iy * dy  - ((float)(ny-1)) / 2.0 * dy;

    density[is] = e0 * ( (float)rand() / RAND_MAX) * exp(-(1.0 / bx) * (x * x)) * exp(-(1.0 / by) * (y * y)) ;
  }
}


void readEnergyDensityBlock(float *density, parameters params)
{
  float e_scale_factor = 1.;
  //float lower_tolerance = 1.0e-3;
  int nx = params.nx;
  int ny = params.ny;

  //first read in the transverse profile
  float temp = 0.0;
  std::ifstream blockFile;
  blockFile.open("initial_profiles/e.dat");
  if (blockFile.is_open())
  {
    //skip the eight line header
    std::string line;
    for (int l = 0; l < 8; l++) getline(blockFile, line);

    for (int iy = 0; iy < ny; iy++)
    {
      for (int ix = 0; ix < nx; ix++)
      {
        blockFile >> temp;
        int is = (ny) * ix + iy; //the column packed index spanning x, y
        density[is] = temp * e_scale_factor;
      }
    }
  }

  else
  {
    printf("Could not find initial profile in initial_profiles!");
  }
  blockFile.close();
}

void initialize2Gaussians(float *density, float bx, float by, parameters params) // bx is the x variance etc...
{

  int ntot = params.ntot;
  int nx = params.nx;
  int ny = params.ny;
  float dx = params.dx;
  float dy = params.dy;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (int is = 0; is < ntot; is++)
  {
    int ix = is / (ny);
    int iy = is - (ny * ix);

    //does it work for even number of points?
    float x = (float)ix * dx  - ((float)(nx-1)) / 2.0 * dx;
    float y = (float)iy * dy  - ((float)(ny-1)) / 2.0 * dy;

    float x1 = -1.0;
    float y1 = 0.0;

    float x2 = 1.0;
    float y2 = 0.0;
    density[is] = e0 * (exp(-(1.0 / bx) * ((x-x1) * (x-x1))) * exp(-(1.0 / by) * ((y-y1) * (y-y1)))
      + exp(-(1.0 / bx) * ((x-x2) * (x-x2))) * exp(-(1.0 / by) * ((y-y2) * (y-y2)))  );
  }

}



void initializeHomogeneous(float *density, parameters params) // bx is the x variance etc...
{
  int ntot = params.ntot;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (int is = 0; is < ntot; is++)
  {
    density[is] = e0;
  }
}

int initializeEnergyDensity(float *energyDensity, std::vector<float> init_energy_density, parameters params)
{
  float hbarc = params.hbarc;

  //initialize energy density
  //define a lower bound on energy density for all cells to regulate numerical noise in flow velocity in dilute regions
  float lower_tolerance = 1.0e-7;

  int option = params.ic_energy;

  printf("setting initial conditions on energy density : ");
  if (option == 1)
  {
    initializeEllipticalGauss(energyDensity, 1.0, 2.0, params);
    printf("Smooth Oblate Gaussian \n");
  }
  else if (option == 2)
  {
    initializeEllipticalMCGauss(energyDensity, 0.5, 1.0, params);
    printf("Fluctuating Oblate Gaussian \n");
  }
  else if (option == 3)
  {
    readEnergyDensityBlock(energyDensity, params);
    printf("Reading from energy density file in initial_profiles/ \n");
  }
  else if (option == 5)
  {
    //read in initial energy density using the initiliaze_from_vector() function
    //note that this is not safe - if one passes an empty vector it will not throw an error
    //converting units of energy density from GeV / fm^3 to fm^(-4)
    printf("Reading energy density from initial energy density vector\n");
    //do a value copy
    //try adding a small value everywhere to regulate problems with flow velocity in dilute regions
    for (int i = 0; i < params.ntot; i++) energyDensity[i] = init_energy_density[i] / (float)hbarc + lower_tolerance;
  }
  else if (option == 6)
  {
    initializeHomogeneous(energyDensity, params);
    printf("Initializing energy density uniform in transverse plane \n");
  }
  else if (option == 7)
  {
    initialize2Gaussians(energyDensity, 1.0, 1.0, params);
    printf("Initializing energy density as two Guassians \n");
  }
  else
  {
    printf("Not a valid initial Condition... Goodbye\n");
    return -1;
  }

  return 1;
}

int initializeFlow(float **flowVelocity, parameters params)
{
  int option = params.ic_flow;

  printf("setting initial conditions on flow velocity : ");
  if (option == 1)
  {
    initializeFlowZero(flowVelocity, params);
    printf("zero \n");
  }
  else if (option == 2)
  {
    readFlowVelocityBlock(flowVelocity, params);
    printf("reading from initial_profiles \n");
  }
  else
  {
    printf("Not a valid initial Condition... Goodbye\n");
    return -1;
  }

  return 1;
}
