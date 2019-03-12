#pragma once
#include <math.h>
#include <stdio.h>
#include "Parameter.h"
#include "NTScheme.cpp"
#include <iostream>

#ifdef _OPENACC
#include <accelmath.h>
#endif

//this creates the initial F function, a function of spatial coordinates and momentum angles
void initializeDensity(float *energyDensity, float **density, parameters params)
{
  int DIM = params.DIM;
  int DIM_PHIP = params.DIM_PHIP;

  //#pragma omp parallel for simd
  for (int is = 0; is < DIM; is++)
  {
    float val = energyDensity[is] / (2.0 * M_PI);
    for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      density[is][iphip] = val;
    }
  } //for (int is = 0; is < DIM; is++)
}


void propagateX(float **density, float **density_p, float *energyDensity, float **flowVelocity, parameters params)
{
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_PHIP = params.DIM_PHIP;
  float tau_iso = params.TAU_ISO;
  float dt = params.DT;
  float dx = params.DX;
  float dy = params.DY;
  float theta = params.THETA;

  //using is = DIM_Y * ix + iy
  int x_stride = DIM_Y;
  int y_stride = 1;

  //update the density moment F based on ITA Eqns of Motion
  //#pragma omp parallel for simd
  for (int is = 0; is < DIM; is++)
  {

    int is_r, is_l, is_t, is_b;
    if (is + x_stride < DIM) is_r = is + x_stride;
    else is_r = is;

    if (is - x_stride > 0 ) is_l = is - x_stride;
    else is_l = is;

    for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      float F = density_p[is][iphip];
      float F_px, F_mx, F_py, F_my;

      F_px = density_p[is_r][iphip];
      F_mx = density_p[is_l][iphip];

      float phip = float(iphip) * (2.0 * M_PI) / float(DIM_PHIP);
      float vx = cos(phip);
      float vy = sin(phip);

      /////////////////////////////////////
      // MacCormack Method
      float F_updated = 0.0;

      //first evolve in x direction
      //predictor step (forward differences) for gradients
      float F_pred_x = F - dt * (vx * (F_px - F) / dx);
      float F_pred_mx = F_mx - dt * (vx * (F - F_mx) / dx);
      //the average
      float F_avg_x = (F + F_pred_x) / 2.0;
      //corrector step for gradients
      float F_corr_x = F_avg_x - dt * (vx * (F_pred_x - F_pred_mx) / 2.0 / dx);

      F_updated = F_corr_x;

      // MacCormack Method
      /////////////////////////////////////

      //update the value of F(x; phip)
      density[is][iphip] = F_updated;

    } //for (int iphip; iphip < DIM_PHIP; iphip++)
  } //for (int is = 0; is < DIM; is++)

}

void propagateY(float **density, float **density_p, float *energyDensity, float **flowVelocity, parameters params)
{
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_PHIP = params.DIM_PHIP;
  float tau_iso = params.TAU_ISO;
  float dt = params.DT;
  float dx = params.DX;
  float dy = params.DY;
  float theta = params.THETA;

  //using is = DIM_Y * ix + iy
  int x_stride = DIM_Y;
  int y_stride = 1;

  //update the density moment F based on ITA Eqns of Motion
  //#pragma omp parallel for simd
  for (int is = 0; is < DIM; is++)
  {
    //flow
    float u0 = flowVelocity[0][is];
    float ux = flowVelocity[1][is];
    float uy = flowVelocity[2][is];

    float eps = energyDensity[is];

    //the isotropic moment F_iso(x;p)
    float F_iso = eps / ( u0*u0 + 0.5 * (ux*ux + uy*uy) );
    F_iso /= (2.0 * M_PI);

    int is_r, is_l, is_t, is_b;
    if (is + y_stride < DIM) is_t = is + y_stride;
    else is_t = is;

    if (is - y_stride > 0) is_b = is - y_stride;
    else is_b = is;

    for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      float F = density_p[is][iphip];
      float F_px, F_mx, F_py, F_my;

      F_py = density_p[is_t][iphip];
      F_my = density_p[is_b][iphip];

      float phip = float(iphip) * (2.0 * M_PI) / float(DIM_PHIP);
      float vx = cos(phip);
      float vy = sin(phip);

      //collision term
      float udotv = u0 - ux*vx - uy*vy;
      float delta_F = F - F_iso;
      float coll = -delta_F * udotv / tau_iso; //the collision term C[F]
      //TEMPORARY FOR TESTING
      coll = 0.0;
      //TEMPORARY FOR TESTING

      /////////////////////////////////////
      // MacCormack Method
      float F_updated = 0.0;

      //evolve in y direction
      //predictor step (forward differences) for gradients
      float F_pred_y = F - dt * (vy * (F_py - F) / dy);
      float F_pred_my = F_my - dt * (vy * (F - F_my) / dy);
      //the average
      float F_avg_y = (F + F_pred_y) / 2.0;
      //corrector step for gradients
      float F_corr_y = F_avg_y - dt * (vy * (F_pred_y - F_pred_my) / 2.0 / dy);

      F_updated = F_corr_y;

      // MacCormack Method
      /////////////////////////////////////

      // add the collision term
      F_updated = F_updated + dt * coll;

      //update the value of F(x; phip)
      density[is][iphip] = F_updated;

    } //for (int iphip; iphip < DIM_PHIP; iphip++)
  } //for (int is = 0; is < DIM; is++)

}
