#pragma once
#include <unistd.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Parameter.h"
#include <math.h>

void writeScalarToFile(float *var, char name[255], parameters params)
{
  int nx = params.nx;
  int ny = params.ny;
  float dx = params.dx;
  float dy = params.dy;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int ix = 0; ix < nx; ix++)
  {
    for (int iy = 0; iy < ny; iy++)
    {
        float x = (float)ix * dx  - (((float)(nx-1)) / 2.0 * dx);
        x = dx * roundf(x / dx);
        float y = (float)iy * dy  - (((float)(ny-1)) / 2.0 * dy);
        y = dy * roundf(y / dy);

        int is = (ny) * ix +  iy; //the column packed index spanning x, y,

        myfile << x << " " << y << " " << var[is] << "\n";
    }
  }
  myfile.close();
}

void writeVectorToFile(float **var, char name[255], int idx, parameters params)
{
  int nx = params.nx;
  int ny = params.ny;
  float dx = params.dx;
  float dy = params.dy;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int ix = 0; ix < nx; ix++)
  {
    for (int iy = 0; iy < ny; iy++)
    {

        float x = (float)ix * dx  - (((float)(nx-1)) / 2.0 * dx);
        x = dx * roundf(x / dx); //rounding for regularly spaced values
        float y = (float)iy * dy  - (((float)(ny-1)) / 2.0 * dy);
        y = dy * roundf(y / dy);

        int is = (ny) * ix + iy; //the column packed index spanning x, y,
        myfile << x << " " << y << " "  << var[idx][is] << "\n";
    }
  }
  myfile.close();
}

//this function writes the transverse density of a variable at z = 0
// as regularly spaced values
void writeScalarToFileProjection(float *var, char name[255], parameters params)
{
  int nx = params.nx;
  int ny = params.ny;

  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int iy = 0; iy < ny; iy++)
  {
    for (int ix = 0; ix < nx; ix++)
    {
      int is = (ny) * ix + iy; //the column packed index spanning x, y, eta
      myfile << var[is] << " "; //different columns for x values
    }
    myfile << "\n"; // different rows correspond to different y values
  }
  myfile.close();
}

void writeVectorToFileProjection(float **var, char name[255], int idx, parameters params)
{
  int nx = params.nx;
  int ny = params.ny;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int iy = 0; iy < ny; iy++)
  {
    for (int ix = 0; ix < nx; ix++)
    {
      int is = (ny) * ix + iy; //the column packed index spanning x, y,
      myfile << var[idx][is] << " "; //different columns for x values
    }
    myfile << "\n"; // different rows correspond to different y values
  }
  myfile.close();
}


/*
void readDensityFile(float *density, char name[255], parameters params)
{

  int nx = params.nx;
  int ny = params.ny;
  int ntot_ETA = params.ntot_ETA;
  float dx = params.dx;
  float dy = params.dy;
  float DETA = params.DETA;
  float xmin = (-1.0) * ((float)(nx-1) / 2.0) * dx;
  float ymin = (-1.0) * ((float)(ny-1) / 2.0) * dy;
  float etamin = (-1.0) * ((float)(ntot_ETA-1) / 2.0) * DETA;
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
    int ix = (int)round((x - xmin) / dx);
    int iy = (int)round((y - ymin) / dy);
    int ieta = (int)round((eta - etamin) / DETA);
    int is = (ny * ntot_ETA * ix) + (ntot_ETA * iy) + ieta;
    density[is] = value;
  }
  infile.close();

}
*/

void readInParameters(struct parameters &params)
{
  char dummyChar[255];
  int dummyInt;
  //long dummyLong;
  float dummyFloat;

  FILE *fileIn;
  std::stringstream paramsStream;
  paramsStream << "ita_input";
  fileIn = fopen(paramsStream.str().c_str(),"r");

  if (fileIn == NULL)
  {
    printf("Couldn't open parameters.dat . Using default values!\n");
  }

  else
  {
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.output_format = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.ic_energy = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.ic_flow = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.nx = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.ny = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.nphip = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.nvz = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.nt = dummyInt;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.dx = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.dy = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.dt = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.t0 = dummyFloat;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.eos_type = dummyInt;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.e_sw = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.eta_over_s = dummyFloat;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.collisions = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.sources = dummyInt;

    fclose(fileIn);
  }
}
