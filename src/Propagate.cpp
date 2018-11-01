#pragma once
#include <math.h>
#include <stdio.h>
#include "Parameter.h"

#ifdef _OPENACC
#include <accelmath.h>
#endif

#define PI 3.141592654f

//this creates the initial F function, a function of spatial coordinates and momentum angles
void initializeDensity(float *energyDensity, float **density, parameters params)
{
  int DIM = params.DIM;
  int DIM_PHIP = params.DIM_PHIP;

  //#pragma omp parallel for simd
  for (int is = 0; is < DIM; is++)
  {
    float val = energyDensity[is] / (2.0 * PI);
    for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      density[is][iphip] = val;
    }
  } //for (int is = 0; is < DIM; is++)
}

void propagate(float **density, float **density_p, float *energyDensity, float **flowVelocity, parameters params)
{

  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_PHIP = params.DIM_PHIP;
  float gamma = params.GAMMA;
  float dt = params.DT;
  float dx = params.DX;
  float dy = params.DY;

  //using is = DIM_Y * ix + iy
  int x_stride = DIM_Y;
  int y_stride = 1;

  //update the density moment F based on ITA Eqns of Motion
  //#pragma omp parallel for simd
  for (int is = 0; is < DIM; is++)
  {

    float u0 = flowVelocity[0][is];
    float ux = flowVelocity[1][is];
    float uy = flowVelocity[2][is];

    float eps = energyDensity[is];
    float eps_rt4 = pow(eps, 1.0/4.0);

    for (int iphip; iphip < DIM_PHIP; iphip++)
    {
      float F = density_p[is][iphip];
      float F_px, F_mx, F_py, F_my;

      if (is + x_stride < DIM) F_px = density_p[is + x_stride][iphip];
      else F_px = F;

      if (is - x_stride > 0 ) F_mx = density_p[is - x_stride][iphip];
      else F_mx = F;

      if (is + y_stride < DIM) F_py = density_p[is + y_stride][iphip];
      else F_py = F;

      if (is - y_stride > 0) F_my = density_p[is - y_stride][iphip];
      else F_my = F;

      //using central differences for now
      float dF_dx = (F_px - F_mx) / (2.0 * dx);
      float dF_dy = (F_py - F_my) / (2.0 * dy);

      float phip = float(iphip) * (2.0 * PI) / float(DIM_PHIP);
      float vx = cos(phip);
      float vy = sin(phip);

      float udotv = u0 - ux*vx - uy*vy;
      float udotv4 = pow(udotv, 4.0);
      float F_iso = eps / udotv4;
      float diff = F - F_iso;
      float G = -( vx*dF_dx + vy*dF_dy ) - gamma * eps_rt4 * (-udotv) * diff;

      //TEMPORARY REMOVE AFTER TESTING
      G = 0.0;
      //TEMPORARY
      
      density[is][iphip] = F + dt * G;
    }
  }
}
