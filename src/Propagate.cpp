#pragma once
#include <math.h>
#include <stdio.h>
#include "Parameter.h"
#include "NTScheme.cpp"
#include <iostream>

#ifdef _OPENACC
#include <accelmath.h>
#endif

//this creates the initial F function, a function of spatial coordinates and momentum phip and velocity vz
void initializeDensity(float *energyDensity, float ***density, parameters params)
{
  int DIM = params.DIM;
  int DIM_PHIP = params.DIM_PHIP;
  int DIM_VZ = params.DIM_VZ;

  #pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    float val = energyDensity[is];
    for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      //try intializing dependence on v_z by a flat distribution
      //alternatively try gaussian or delta function ?
      for (int ivz = 0; ivz < DIM_VZ; ivz++)
      {
        density[is][iphip][ivz] = val;
      }

    }
  } //for (int is = 0; is < DIM; is++)
}


void propagateX(float ***density, float ***density_p, float *energyDensity, float **flowVelocity, parameters params)
{
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_PHIP = params.DIM_PHIP;
  float dt = params.DT;
  float dx = params.DX;
  float dy = params.DY;
  int DIM_VZ = params.DIM_VZ;

  //using is = DIM_Y * ix + iy
  int x_stride = DIM_Y;
  int y_stride = 1;

  //update the density moment F based on ITA Eqns of Motion
  #pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    int is_r, is_l, is_t, is_b;
    if (is + x_stride < DIM) is_r = is + x_stride;
    else is_r = is;

    if (is - x_stride > 0 ) is_l = is - x_stride;
    else is_l = is;

    for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      float phip = float(iphip) * (2.0 * M_PI) / float(DIM_PHIP);
      float vx = cos(phip);
      float vy = sin(phip);

      for (int ivz = 0; ivz < DIM_VZ; ivz++)
      {
        float F = density_p[is][iphip][ivz];
        float F_px, F_mx, F_py, F_my;

        F_px = density_p[is_r][iphip][ivz];
        F_mx = density_p[is_l][iphip][ivz];

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
        density[is][iphip][ivz] = F_updated;
      }
    } //for (int iphip; iphip < DIM_PHIP; iphip++)
  } //for (int is = 0; is < DIM; is++)

}

void propagateY(float ***density, float ***density_p, float *energyDensity, float **flowVelocity, parameters params)
{
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_PHIP = params.DIM_PHIP;
  float dt = params.DT;
  float dx = params.DX;
  float dy = params.DY;
  int DIM_VZ = params.DIM_VZ;

  //using is = DIM_Y * ix + iy
  int x_stride = DIM_Y;
  int y_stride = 1;

  //update the density moment F based on ITA Eqns of Motion
  #pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    int is_r, is_l, is_t, is_b;
    if (is + y_stride < DIM) is_t = is + y_stride;
    else is_t = is;

    if (is - y_stride > 0) is_b = is - y_stride;
    else is_b = is;

    for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      float phip = float(iphip) * (2.0 * M_PI) / float(DIM_PHIP);
      float vx = cos(phip);
      float vy = sin(phip);

      for (int ivz = 0; ivz < DIM_VZ; ivz++)
      {
        float F = density_p[is][iphip][ivz];
        float F_px, F_mx, F_py, F_my;

        F_py = density_p[is_t][iphip][ivz];
        F_my = density_p[is_b][iphip][ivz];

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

        //update the value of F(x; phip)
        density[is][iphip][ivz] = F_updated;
      }
    } //for (int iphip; iphip < DIM_PHIP; iphip++)
  } //for (int is = 0; is < DIM; is++)
}


void propagateVz(float ***density, float ***density_p, float *energyDensity, float **flowVelocity, float tau, parameters params)
{
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_PHIP = params.DIM_PHIP;
  float alpha = params.ALPHA;
  float dt = params.DT;
  float dx = params.DX;
  float dy = params.DY;
  int DIM_VZ = params.DIM_VZ;
  float dvz = 2.0 / (float)(DIM_VZ);
  float dvz_2 = 2.0 / (float)(DIM_VZ - 1);

  //update the density moment F based on ITA Eqns of Motion
  #pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      for (int ivz = 0; ivz < DIM_VZ; ivz++)
      {
        float vz = -1.0 + (float)ivz * dvz_2; 

        int ivz_l = ivz - 1;
        int ivz_r = ivz + 1;
        if (ivz_l < 0) ivz_l = ivz;
        if (ivz_r > DIM_VZ - 1) ivz_r = ivz;

        float F = density_p[is][iphip][ivz];
        float F_pvz = density_p[is][iphip][ivz_r];
        float F_mvz = density_p[is][iphip][ivz_l];

        /////////////////////////////////////
        // MacCormack Method
        float F_updated = 0.0;
        float avz = -vz * (1.0 - vz*vz) / tau; //coefficient of partial F / partial v_z

        //predictor step (forward differences) for gradients
        float F_pred_vz = F - dt * (avz * (F_pvz - F) / dvz);
        float F_pred_mvz = F_mvz - dt * (avz * (F - F_mvz) / dvz);
        //the average
        float F_avg_vz = (F + F_pred_vz) / 2.0;
        //corrector step for gradients
        float F_corr_vz = F_avg_vz - dt * (avz * (F_pred_vz - F_pred_mvz) / 2.0 / dvz);

        F_updated = F_corr_vz;

        // MacCormack Method
        /////////////////////////////////////

        //add the geometric source term
        float geom_src = (4.0 * vz * vz) * F / tau;
        F_updated = F_updated + dt * geom_src;

        //update the value of F(x; phip)
        density[is][iphip][ivz] = F_updated;
      }
    } //for (int iphip; iphip < DIM_PHIP; iphip++)
  } //for (int is = 0; is < DIM; is++)
}

void propagateColl(float ***density, float ***density_p, float *energyDensity, float **flowVelocity, parameters params)
{
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_PHIP = params.DIM_PHIP;
  float alpha = params.ALPHA;
  float dt = params.DT;
  float dx = params.DX;
  float dy = params.DY;
  int DIM_VZ = params.DIM_VZ;
  //float dvz = 2.0 / (DIM_VZ - 1);

  //update the density moment F based on ITA Eqns of Motion
  #pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    float u0 = flowVelocity[0][is];
    float ux = flowVelocity[1][is];
    float uy = flowVelocity[2][is];

    float eps = energyDensity[is];

    // EoS : \eps = a T^4
    // tau_iso = \alpha / T
    //float a = 15.6269; // Nc=3, Nf=3
    float a = 13.8997; // Nc=3, Nf=2.5
    float T = powf( (eps/a), 0.25);
    float tau_iso = alpha / T;

    if (tau_iso < 10.0 * dt) printf("Warning: tau_iso = %f , energy density = %f , take smaller dt! \n", tau_iso, eps);

    for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      float phip = float(iphip) * (2.0 * M_PI) / float(DIM_PHIP);
      float vx = cos(phip);
      float vy = sin(phip);

      for (int ivz = 0; ivz < DIM_VZ; ivz++)
      {
        float F = density_p[is][iphip][ivz];

        float F_updated = 0.0;
        //collision term
        float udotv = u0 - ux*vx - uy*vy;
        float F_iso = eps / powf(udotv, 4.0); //the isotropic moment F_iso(x;p),  check factors of 4pi everywhere!!!
        float delta_F = F - F_iso;
        float coll = -delta_F * udotv / tau_iso; //the collision term C[F]

        F_updated = F + dt * coll;

        //update the value of F(x; phip)
        density[is][iphip][ivz] = F_updated;
      }

    } //for (int iphip; iphip < DIM_PHIP; iphip++)
  } //for (int is = 0; is < DIM; is++)
}
