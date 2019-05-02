#pragma once
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include "Parameter.h"
#define THETA_FUNCTION(X) ((float)X < (float)0 ? (float)0 : (float)1)

void initializeZero(float *density, parameters params)
{
  int DIM = params.DIM;
  for (int is = 0; is < DIM; is++)
  {
    density[is] = 0.0;
  }
}

void initializeFlow(float **flowVelocity, parameters params)
{
  int DIM = params.DIM;
  for (int is = 0; is < DIM; is++)
  {
    flowVelocity[0][is] = 1.0;
    flowVelocity[1][is] = 0.0;
    flowVelocity[2][is] = 0.0;
    flowVelocity[3][is] = 0.0;
  }
}

void initializeGauss(float *density, float b, parameters params) // b is the variance ('spherically' symmetric)
{
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  float DX = params.DX;
  float DY = params.DY;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y );
    int iy = (is - (DIM_Y * ix));

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;

    density[is] = e0 * exp( -(1.0 / b) * ( (x * x) + (y * y) ) );
  }
}

void initializeEllipticalGauss(float *density, float bx, float by, parameters params) // bx is the x variance etc...
{
  float regulate = 1.0e-10;
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  float DX = params.DX;
  float DY = params.DY;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y);
    int iy = (is - (DIM_Y * ix));

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;

    density[is] = e0 * exp(-(3.0 / (8.0 * bx * bx)) * (x * x)) * exp(-(3.0 / (8.0 * by * by)) * (y * y)) + regulate;
  }
}

void initializeMCGauss(float * density, float b, parameters params)
{
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  float DX = params.DX;
  float DY = params.DY;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y);
    int iy = (is - (DIM_Y * ix));

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;

    density[is] = e0 * ( (float)rand() / RAND_MAX ) * exp( -(1.0 / b) * ((x * x) + (y * y) ) );
  }
}

void initializeEllipticalMCGauss(float *density, float bx, float by, parameters params) // bx is the x variance etc...
{
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  float DX = params.DX;
  float DY = params.DY;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y);
    int iy = (is - (DIM_Y * ix));

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;

    density[is] = e0 * ( (float)rand() / RAND_MAX) * exp(-(1.0 / bx) * (x * x)) * exp(-(1.0 / by) * (y * y)) ;
  }
}

void readEnergyDensitySuperMCBlock(float *density, parameters params)
{
  /*
  float lower_tolerance = 1.0e-3;

  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  float ETA_WIDTH = params.ETA_WIDTH;
  float ETA_FLAT = params.ETA_FLAT;
  float DETA = params.DETA;

  //first read in the transverse profile from superMC block data format
  float temp = 0.0;
  std::ifstream superMCFile;
  superMCFile.open("initial_superMC_ed/2.dat");
  if (superMCFile.is_open())
  {
    for (int ix = 0; ix < DIM_X; ix++)
    {
      for (int iy = 0; iy < DIM_Y; iy++)
      {
        superMCFile >> temp;
        for (int ieta = 0; ieta < DIM_ETA; ieta++) //copy the same value for all eta, then we will multiply by eta dependent function
        {
          int is = (DIM_Y * DIM_ETA) * ix + (DIM_ETA) * iy + ieta; //the column packed index spanning x, y, z
          density[is] = temp;
        }
      }
    }
  }

  else
  {
    printf("Could not find initial profile in initial_superMC_ed!");
  }

  superMCFile.close();

  //now multiply by an eta-dependent profile; etaWidth is the width of the eta profile
  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y * DIM_ETA);
    int iy = (is - (DIM_Y * DIM_ETA * ix))/ DIM_ETA;
    int ieta = is - (DIM_Y * DIM_ETA * ix) - (DIM_ETA * iy);

    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;
    //here we use a the same profile as GPU-VH (see arXiv:1608.06577v1 p. 38)
    float arg = (-1.0) * (abs(eta) - ETA_FLAT) * (abs(eta) - ETA_FLAT) / (2.0 * ETA_WIDTH * ETA_WIDTH);
    arg = arg * THETA_FUNCTION(abs(eta) - ETA_FLAT);
    density[is] = density[is] * exp(arg) + lower_tolerance;
  }
  */
}

void readEnergyDensityTRENTOBlock(float *density, parameters params)
{
  float lower_tolerance = 1.0e-3;

  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;

  //first read in the transverse profile from superMC block data format
  float temp = 0.0;
  std::ifstream superMCFile;
  superMCFile.open("initial_profiles/e.dat");
  if (superMCFile.is_open())
  {
    //skip the eight line (l) header
    std::string line;
    for (int l = 0; l < 12; l++) getline(superMCFile, line);
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int ix = 0; ix < DIM_X; ix++)
      {
        superMCFile >> temp;
        int is = (DIM_Y) * ix + iy; //the column packed index spanning x, y, z
        density[is] = temp;
      }
    }
  }

  else
  {
    printf("Could not find initial profile in initial_profiles!");
  }
  superMCFile.close();
}

void initialize2Gaussians(float *density, float bx, float by, float beta, parameters params) // bx is the x variance etc...
{
  /*
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y * DIM_ETA);
    int iy = (is - (DIM_Y * DIM_ETA * ix))/ DIM_ETA;
    int ieta = is - (DIM_Y * DIM_ETA * ix) - (DIM_ETA * iy);

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;
    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;

    float x1 = 3.0;
    float y1 = 0.0;
    float eta1 = 0.0;

    float x2 = -3.0;
    float y2 = 0.0;
    float eta2 = 0.0;
    density[is] = e0 * (exp(-(1.0 / bx) * ((x-x1) * (x-x1))) * exp(-(1.0 / by) * ((y-y1) * (y-y1))) * exp(-(1.0 / beta) * ((eta-eta1) * (eta-eta1)))
      + exp(-(1.0 / bx) * ((x-x2) * (x-x2))) * exp(-(1.0 / by) * ((y-y2) * (y-y2))) * exp(-(1.0 / beta) * ((eta-eta2) * (eta-eta2))) );
  }
  */
}

void readEnergyDensityTRENTO3DBlock(float *density, parameters params)
{
  /*
  float lower_tolerance = 1.0e-3;

  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  float ETA_WIDTH = params.ETA_WIDTH;
  float ETA_FLAT = params.ETA_FLAT;
  float DETA = params.DETA;

  //first read in the transverse profile from superMC block data format
  float temp = 0.0;
  std::ifstream superMCFile;
  superMCFile.open("initial_profiles/e.dat");
  if (superMCFile.is_open())
  {
    //skip the eight line (l) header
    std::string line;
    for (int l = 0; l < 8; l++) getline(superMCFile, line);
    for (int ix = 0; ix < DIM_X; ix++)
    {
      for (int iy = 0; iy < DIM_Y; iy++)
      {
        for (int ieta = 0; ieta < DIM_ETA; ieta++)
        {
          int is = (DIM_Y * DIM_ETA) * ix + (DIM_ETA) * iy + ieta; //the column packed index spanning x, y, z
          superMCFile >> temp;
          density[is] = temp + lower_tolerance;
        }
      }
    }
  }

  else
  {
    printf("Could not find initial profile in initial_profiles!");
  }
  superMCFile.close();
  */
}

void readDensityFile(float *density, char name[255], parameters params)
{
  /*
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;
  float xmin = (-1.0) * ((float)(DIM_X-1) / 2.0) * DX;
  float ymin = (-1.0) * ((float)(DIM_Y-1) / 2.0) * DY;
  float etamin = (-1.0) * ((float)(DIM_ETA-1) / 2.0) * DETA;
  float x, y, eta, value;

  char filename[255] = "";
  sprintf(filename, "%s.dat", name);
  std::ifstream infile;
  infile.open(filename);
  if (!infile)
  {
    printf("Couldn't open initial profile!\n");
    exit(1);
  }
  while (infile >> x >> y >> eta >> value)
  {
    int ix = (int)round((x - xmin) / DX);
    int iy = (int)round((y - ymin) / DY);
    int ieta = (int)round((eta - etamin) / DETA);
    int is = (DIM_Y * DIM_ETA * ix) + (DIM_ETA * iy) + ieta;
    density[is] = value;
  }
  infile.close();
  */
}

void initializeHomogeneous(float *density, parameters params) // bx is the x variance etc...
{
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  float DX = params.DX;
  float DY = params.DY;

  float e0 = 500.0; //energy norm factor in fm^(-4) : roughly 500 MeV Temperature

  for (int is = 0; is < DIM; is++)
  {
    density[is] = e0;
  }
}

int initializeEnergyDensity(float *energyDensity, std::vector<float> init_energy_density, parameters params)
{
  float hbarc = 0.197326938;

  //initialize energy density
  //define a lower bound on energy density for all cells to regulate numerical noise in flow velocity in dilute regions
  float lower_tolerance = 1.0e-7;

  int option = params.IC_ENERGY;

  printf("setting initial conditions on energy density : ");
  if (option == 1)
  {
    initializeEllipticalGauss(energyDensity, 1.0, 1.0, params);
    printf("Smooth Oblate Gaussian \n");
  }
  else if (option == 2)
  {
    initializeEllipticalMCGauss(energyDensity, 3.0, 3.0,params);
    printf("Fluctuating Oblate Gaussian \n");
  }
  else if (option == 3)
  {
    readDensityFile(energyDensity, (char *)"initial_profiles/e", params);
    printf("Reading from energy density file in initial_profiles/ \n");
  }
  else if (option == 4)
  {
    readEnergyDensitySuperMCBlock(energyDensity, params);
    printf("Reading from superMC energy density file in initial_profiles/ \n");
  }
  else if (option == 5)
  {
    //read in initial energy density using the initiliaze_from_vector() function
    //note that this is not safe - if one passes an empty vector it will not throw an error
    //converting units of energy density from GeV / fm^3 to fm^(-4)
    printf("Reading energy density from initial energy density vector\n");
    //do a value copy
    //try adding a small value everywhere to regulate problems with flow velocity in dilute regions
    for (int i = 0; i < params.DIM; i++) energyDensity[i] = init_energy_density[i] / (float)hbarc + lower_tolerance;
  }
  else if (option == 6)
  {
    initializeHomogeneous(energyDensity, params);
    printf("Initializing energy density uniform in transverse plane \n");
  }
  else
  {
    printf("Not a valid initial Condition... Goodbye\n");
    return -1;
  }

}
