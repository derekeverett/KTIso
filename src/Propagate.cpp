#pragma once
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <iostream>

#include "Parameter.h"
#include "NTScheme.cpp"
#include "EquationOfState.cpp"

#ifdef _OPENACC
#include <accelmath.h>
#endif

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

void updateDensity(float ***density, float ***density_p, parameters params)
{
  //do a value swap
  for (int is = 0; is < params.DIM; is++)
  {
    for (int iphip = 0; iphip < params.DIM_PHIP; iphip++)
    {
      for (int ivz = 0; ivz < params.DIM_VZ; ivz++)
      {
        density_p[is][iphip][ivz] = density[is][iphip][ivz];
      } //for (int ivz = 0; ivz < params.DIM_VZ; ivz++)
    } //for (int iphip = 0; iphip < params.DIM_PHIP; iphip++)
  } //for (int is = 0; is < params.DIM; is++)
}
void propagateX(float ***density, float ***density_p, float dt, parameters params)
{
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_PHIP = params.DIM_PHIP;
  float dx = params.DX;
  int DIM_VZ = params.DIM_VZ;

  //using is = DIM_Y * ix + iy
  int x_stride = DIM_Y;

  //update the density moment F based on ITA Eqns of Motion
  #pragma omp parallel for
  for (int ix = 2; ix < DIM_X - 2; ix++)
  {
    for (int iy = 2; iy < DIM_Y - 2; iy++)
    {
      int is = (DIM_Y * ix) + iy;
      int is_r = is + x_stride;
      int is_l = is - x_stride;

      for (int iphip = 0; iphip < DIM_PHIP; iphip++)
      {
        float phip = float(iphip) * (2.0 * M_PI) / float(DIM_PHIP);
        float vx = cos(phip);

        for (int ivz = 0; ivz < DIM_VZ; ivz++)
        {
          float F = density_p[is][iphip][ivz];
          float F_px = density_p[is_r][iphip][ivz];
          float F_mx = density_p[is_l][iphip][ivz];

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
          density[is][iphip][ivz] = fmax(F_updated, 1.0e-10);
        }
      } //for (int iphip; iphip < DIM_PHIP; iphip++)
    } //for (int iy = 0; iy < DIM_Y; iy++)
  } //for (int ix = 0; ix < DIM_X; ix++)

}

void propagateY(float ***density, float ***density_p, float dt, parameters params)
{
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_PHIP = params.DIM_PHIP;
  float dy = params.DY;
  int DIM_VZ = params.DIM_VZ;

  //using is = DIM_Y * ix + iy
  int y_stride = 1;

  //update the density moment F based on ITA Eqns of Motion
  #pragma omp parallel for
  for (int ix = 2; ix < DIM_X - 2; ix++)
  {
    for (int iy = 2; iy < DIM_Y - 2; iy++)
    {
      int is = (DIM_Y * ix) + iy;
      int is_t = is + y_stride;
      int is_b = is - y_stride;

      for (int iphip = 0; iphip < DIM_PHIP; iphip++)
      {
        float phip = float(iphip) * (2.0 * M_PI) / float(DIM_PHIP);
        float vy = sin(phip);

        for (int ivz = 0; ivz < DIM_VZ; ivz++)
        {
          float F = density_p[is][iphip][ivz];
          float F_py = density_p[is_t][iphip][ivz];
          float F_my = density_p[is_b][iphip][ivz];

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
          density[is][iphip][ivz] = fmax(F_updated, 1.0e-10);
        }
      } //for (int iphip; iphip < DIM_PHIP; iphip++)
    } //for (int iy = 2; iy < DIM_Y - 2; iy++)
  } //for (int ix = 2; ix < DIM_X - 2; ix++)
}


void propagateVz(float ***density, float ***density_p, float **vz_quad, float tau, float dt, parameters params)
{
  int DIM = params.DIM;
  int DIM_PHIP = params.DIM_PHIP;
  //float dt = params.DT;
  int DIM_VZ = params.DIM_VZ;
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
        //float vz = vz_quad[ivz][0];
        //float dvz = vz_quad[ivz][1];

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
        //note this coeff diverges as tau -> 0 !
        //this term is zero for v_z = 0, and for v_z = +/- 1
        float avz = -vz * (1.0 - vz*vz) / tau; //coefficient of partial F / partial v_z

        //predictor step (forward differences) for gradients
        float F_pred_vz = F - dt * (avz * (F_pvz - F) / dvz_2);
        //float F_pred_vz = F - dt * (avz * (F_pvz - F) / dvz);
        float F_pred_mvz = F_mvz - dt * (avz * (F - F_mvz) / dvz_2);
        //float F_pred_mvz = F_mvz - dt * (avz * (F - F_mvz) / dvz);
        //the average
        float F_avg_vz = (F + F_pred_vz) / 2.0;
        //corrector step for gradients
        float F_corr_vz = F_avg_vz - dt * (avz * (F_pred_vz - F_pred_mvz) / 2.0 / dvz_2);
        //float F_corr_vz = F_avg_vz - dt * (avz * (F_pred_vz - F_pred_mvz) / 2.0 / dvz);

        F_updated = F_corr_vz;

        // MacCormack Method
        /////////////////////////////////////

        //update the value of F(x; phip)
        density[is][iphip][ivz] = fmax(F_updated, 1.0e-10);

      }
    } //for (int iphip; iphip < DIM_PHIP; iphip++)
  } //for (int is = 0; is < DIM; is++)
}

//propagate the geometric source term
void propagateVzGeom(float ***density, float ***density_p, float **vz_quad, float tau, float dt, parameters params)
{
  int DIM = params.DIM;
  int DIM_PHIP = params.DIM_PHIP;
  //float dt = params.DT;
  int DIM_VZ = params.DIM_VZ;
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
        //float vz = vz_quad[ivz][0];

        int ivz_l = ivz - 1;
        int ivz_r = ivz + 1;
        if (ivz_l < 0) ivz_l = ivz;
        if (ivz_r > DIM_VZ - 1) ivz_r = ivz;

        float F = density_p[is][iphip][ivz];
        //this term is zero for v_z = 0, and largest for v_z = +/-1
        float geom_src = -(4.0 * vz * vz) * F / tau;

        //add geometric source term with RK2 forward step
        float k1 = geom_src;
        //estimate value at t + dt/2
        float y1 = F + k1 * (dt / 2.0);
        //estimate slope at t+dt/2
        float k2 = (y1 - F) / (dt / 2.0);
        //estimate value at t+dt
        float y2 = F + k2 * dt;

        //update the value of F(x; phip)
        density[is][iphip][ivz] = y2;
      }
    } //for (int iphip; iphip < DIM_PHIP; iphip++)
  } //for (int is = 0; is < DIM; is++)
}

//sets to zero the density on all edges
void propagateBoundaries(float ***density, parameters params)
{
  int DIM = params.DIM;
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_PHIP = params.DIM_PHIP;
  int DIM_VZ = params.DIM_VZ;

  //along boundaries in y
  for (int ix = 0; ix < DIM_X; ix++)
  {
    for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      for (int ivz = 0; ivz < DIM_VZ; ivz++)
      {
        int iyll = 0;
        int iyl = 1;
        int iy_b = 2;

        int iy_t = DIM_Y - 3;
        int iyr = DIM_Y - 2;
        int iyrr = DIM_Y - 1;

        int isll = DIM_Y * ix + iyll;
        int isl = DIM_Y * ix + iyl;
        int isr = DIM_Y * ix + iyr;
        int isrr = DIM_Y * ix + iyrr;

        int is_b = DIM_Y * ix + iy_b;
        int is_t = DIM_Y * ix + iy_t;

        density[isll][iphip][ivz] = density[is_b][iphip][ivz];
        density[isl][iphip][ivz] = density[is_b][iphip][ivz];

        density[isr][iphip][ivz] = density[is_t][iphip][ivz];
        density[isrr][iphip][ivz] = density[is_t][iphip][ivz];
      }
    }
  } //for (int ix = 0; ix < DIM_X; ix++)

  //along boundaries in x
  for (int iy = 0; iy < DIM_Y; iy++)
  {
    for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      for (int ivz = 0; ivz < DIM_VZ; ivz++)
      {
        int ixll = 0;
        int ixl = 1;
        int ix_b = 2;

        int ix_t = DIM_X - 3;
        int ixr = DIM_X - 2;
        int ixrr = DIM_X - 1;

        int isll = DIM_Y * ixll + iy;
        int isl = DIM_Y * ixl + iy;
        int isr = DIM_Y * ixr + iy;
        int isrr = DIM_Y * ixrr + iy;

        int is_b = DIM_Y * ix_b + iy;
        int is_t = DIM_Y * ix_t + iy;

        density[isll][iphip][ivz] = density[is_b][iphip][ivz];
        density[isl][iphip][ivz] = density[is_b][iphip][ivz];

        density[isr][iphip][ivz] = density[is_t][iphip][ivz];
        density[isrr][iphip][ivz] = density[is_t][iphip][ivz];
      }
    }
  } //for (int iy = 0; iy < DIM_Y; iy++)

  //along boundaries in vz
  if (DIM_VZ > 1)
  {
    for (int is = 0; is < DIM; is++)
    {
      for (int iphip = 0; iphip < DIM_PHIP; iphip++)
      {
        density[is][iphip][0] = density[is][iphip][2];
        density[is][iphip][1] = density[is][iphip][2];
        density[is][iphip][DIM_VZ - 2] = density[is][iphip][DIM_VZ - 3];
        density[is][iphip][DIM_VZ - 1] = density[is][iphip][DIM_VZ - 3];
      }
    } //if (DIM_VZ > 1)
  } // if (DIM_VZ > 1)
}

void propagateCoordSpace(float ***density, float ***density_p, float *energyDensity, float **flowVelocity,
                float **vz_quad, float dt, parameters params)
{
  //propagate the density forward by one time step according to ITA EQN of Motion
  propagateX(density, density_p, dt, params); //propagate x direction
  propagateBoundaries(density, params);

  //do a value swap
  for (int is = 0; is < params.DIM; is++)
  {
    for (int iphip = 0; iphip < params.DIM_PHIP; iphip++)
    {
      for (int ivz = 0; ivz < params.DIM_VZ; ivz++)
      {
        density_p[is][iphip][ivz] = density[is][iphip][ivz];
      }
    }
  }

  //std::swap(density, density_p); //swap the density and previous value
  propagateY(density, density_p, dt, params); //propagate y direction
  propagateBoundaries(density, params);

  //do a value swap
  for (int is = 0; is < params.DIM; is++)
  {
    for (int iphip = 0; iphip < params.DIM_PHIP; iphip++)
    {
      for (int ivz = 0; ivz < params.DIM_VZ; ivz++)
      {
        density_p[is][iphip][ivz] = density[is][iphip][ivz];
      }
    }
  }

  //std::swap(density, density_p);
}

void propagateMomentumSpace(float ***density, float ***density_p, float *energyDensity, float **flowVelocity,
  float **vz_quad, float t, float dt, parameters params)
{
  propagateVz(density, density_p, vz_quad, t, dt, params); // propagate v_z gradient term
  propagateBoundaries(density, params);

  //do a value swap
  for (int is = 0; is < params.DIM; is++)
  {
    for (int iphip = 0; iphip < params.DIM_PHIP; iphip++)
    {
      for (int ivz = 0; ivz < params.DIM_VZ; ivz++)
      {
        density_p[is][iphip][ivz] = density[is][iphip][ivz];
      }
    }
  }

  //std::swap(density, density_p);
  propagateVzGeom(density, density_p, vz_quad, t, dt, params); // propagate v_z geometric source term
  propagateBoundaries(density, params);

  //do a value swap
  for (int is = 0; is < params.DIM; is++)
  {
    for (int iphip = 0; iphip < params.DIM_PHIP; iphip++)
    {
      for (int ivz = 0; ivz < params.DIM_VZ; ivz++)
      {
        density_p[is][iphip][ivz] = density[is][iphip][ivz];
      }
    }
  }

  //std::swap(density, density_p);
}


//this function propagates coordinate space and momentum space according to freestreaming terms
void propagate(float ***density, float ***density_p, float ***density_i,
                float *energyDensity, float **flowVelocity, float **vz_quad, float t, parameters params)
{
  int DIM_VZ = params.DIM_VZ;
  //int COLLISIONS = params.COLLISIONS;

  float dt = params.DT;
  float dt_2 = dt / 2.0;
  //to propagate forward in time, use Strang operator splitting approach to split advection terms in coord space and momentum space
  if (DIM_VZ > 1)
  {
    propagateCoordSpace(density, density_p, energyDensity, flowVelocity, vz_quad, dt_2, params);
    propagateMomentumSpace(density, density_p, energyDensity, flowVelocity, vz_quad, t, dt, params);
    propagateCoordSpace(density, density_p, energyDensity, flowVelocity, vz_quad, dt_2, params);
  }

  else propagateCoordSpace(density, density_p, energyDensity, flowVelocity, vz_quad, dt, params);

  //if (COLLISIONS) propagateCollisionTerms(density, density_p, density_i, energyDensity, flowVelocity, vz_quad, t, dt, params);
}
