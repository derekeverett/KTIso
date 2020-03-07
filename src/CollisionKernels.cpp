#pragma once
#include <math.h>
#include <stdio.h>
#include "Parameter.h"
#include "NTScheme.cpp"
#include "EquationOfState.cpp"
#include <iostream>

#ifdef _OPENACC
#include <accelmath.h>
#endif

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

//This uses a forward Euler step
void propagateITAColl(float ***density, float ***density_p, float *energyDensity, float **flowVelocity, float dt, parameters params)
{
  int ntot = params.ntot;
  int nphip = params.nphip;
  float eta_over_s = params.eta_over_s;
  int nvz = params.nvz;
  //float dvz = 2.0 / (nvz - 1);
  float dvz_2 = 2.0 / (float)(nvz - 1);

  int warn_flag = 1;

  //update the density moment F based on ITA Eqns of Motion
  #pragma omp parallel for
  for (int is = 0; is < ntot; is++)
  {
    float u0 = flowVelocity[0][is];
    float ux = flowVelocity[1][is];
    float uy = flowVelocity[2][is];

    float eps = energyDensity[is];

    float T = temperatureFromEnergyDensity(eps);
    float tau_iso = 5. * eta_over_s / T;

    if ( (tau_iso < 3.0 * dt) && (warn_flag) )
    {
      printf("Warning: tau_iso = %f < 3*dt, energy density = %f , take smaller dt! \n", tau_iso, eps);
      warn_flag = 0;
    }
    for (int iphip = 0; iphip < nphip; iphip++)
    {
      float phip = float(iphip) * (2.0 * M_PI) / float(nphip);
      //float vx = cos(phip);
      //float vy = sin(phip);

      for (int ivz = 0; ivz < nvz; ivz++)
      {
        //float vz = -1.0 + (float)ivz * dvz_2;
        float vz = (nvz > 1) ? -1.0 + (float)ivz * dvz_2 : 0.0;
        float thetap = acos(vz);
        float sin_thetap = (nvz > 1) ? sin(thetap) : 1.0;
        float vx = sin_thetap * cos(phip);
        float vy = sin_thetap * sin(phip);

        float F = density_p[is][iphip][ivz];

        //collision term
        float udotv = u0 - ux*vx - uy*vy;
        float F_iso = eps / powf(udotv, 4.0); //the isotropic moment F_iso(x;p),  check factors of 4pi everywhere!!!
        if (nvz == 1) F_iso = eps / powf(udotv, 4.0) * 2.0;
        float delta_F = F - F_iso;
        float coll = -1.0 * delta_F * udotv / tau_iso; //the collision term C[F]

        //update the value of F(x; phip)
        density[is][iphip][ivz] = F + coll * dt;

        //if (iphip == 0 && is == (ntot - 1) / 2) std::cout << "k1 * dt = " << k1 * dt << std::endl;
      } //for (int ivz = 0; ivz < nvz; ivz++)
    } //for (int iphip; iphip < nphip; iphip++)
  } //for (int is = 0; is < ntot; is++)
}


//this uses forward Euler step. When combined with kernel above, together can be used as RK2 forward step
void propagateITACollConvexComb(float ***density, float ***density_i, float ***density_p, float *energyDensity,
  float **flowVelocity, float dt, parameters params)
{
  int ntot = params.ntot;
  int nphip = params.nphip;
  float eta_over_s = params.eta_over_s;
  int nvz = params.nvz;
  //float dvz = 2.0 / (nvz - 1);
  float dvz_2 = 2.0 / (float)(nvz - 1);

  int warn_flag = 1;

  //update the density moment F based on ITA Eqns of Motion
  #pragma omp parallel for
  for (int is = 0; is < ntot; is++)
  {
    float u0 = flowVelocity[0][is];
    float ux = flowVelocity[1][is];
    float uy = flowVelocity[2][is];

    float eps = energyDensity[is];

    float T = temperatureFromEnergyDensity(eps);
    float tau_iso = 5. * eta_over_s / T;

    if ( (tau_iso < 3.0 * dt) && (warn_flag) )
    {
      printf("Warning: tau_iso = %f < 3*dt, energy density = %f , take smaller dt! \n", tau_iso, eps);
      warn_flag = 0;
    }
    for (int iphip = 0; iphip < nphip; iphip++)
    {
      float phip = float(iphip) * (2.0 * M_PI) / float(nphip);
      //float vx = cos(phip);
      //float vy = sin(phip);

      for (int ivz = 0; ivz < nvz; ivz++)
      {
        //float vz = -1.0 + (float)ivz * dvz_2;
        float vz = (nvz > 1) ? -1.0 + (float)ivz * dvz_2 : 0.0;
        float thetap = acos(vz);
        float sin_thetap = (nvz > 1) ? sin(thetap) : 1.0;
        float vx = sin_thetap * cos(phip);
        float vy = sin_thetap * sin(phip);

        float F = density_i[is][iphip][ivz];

        //collision term
        float udotv = u0 - ux*vx - uy*vy;
        float F_iso = eps / powf(udotv, 4.0); //the isotropic moment F_iso(x;p),  check factors of 4pi everywhere!!!
        if (nvz == 1) F_iso = eps / powf(udotv, 4.0) * 2.0;
        //float delta_F = F - F_iso;
        //float coll = -1.0 * delta_F * udotv / tau_iso; //the collision term C[F]

        //update the value of F(x; phip)
        //density[is][iphip][ivz] = density_p[is][iphip][ivz] + (coll * dt);

        //using exact solution 
        float nu = udotv / tau_iso;
        if (nu * dt > 0.2) printf("Warning: nu * dt = %f > 0.2, energy density = %f , take smaller dt! \n", nu * dt, eps);
        float exp_weight = exp(-1.0 * nu * dt);
        //update the value of F(x; phip)
        density[is][iphip][ivz] = exp_weight * density_p[is][iphip][ivz] + (1.0 - exp_weight) * F_iso;

        //if (iphip == 0 && is == (ntot - 1) / 2) std::cout << "k1 * dt = " << k1 * dt << std::endl;
      } //for (int ivz = 0; ivz < nvz; ivz++)
    } //for (int iphip; iphip < nphip; iphip++)
  } //for (int is = 0; is < ntot; is++)
}


/*
void propagateCollisionTerms(float ***density, float ***density_p, float ***density_i,
                float *energyDensity, float **flowVelocity,
                float **vz_quad, float t, float dt, parameters params)
{
  //RK2
  //first find the current energy density and flow after advection updates
  calculateStressTensor(stressTensor, density_p, hypertrigTable, vz_quad, t, params);
  solveEigenSystem(stressTensor, energyDensity, flowVelocity, params);

  //guess step to estimate the collision term at C[F(t + dt/2)]
  propagateITAColl(density_i, density_p, energyDensity, flowVelocity, dt / 2.0, params);
  propagateBoundaries(density_i, params);

  //now propagate using the estimated C[F(t + dt/2)]
  propagateITACollConvexComb(density, density_i, density_p, energyDensity, flowVelocity, dt, params);
  propagateBoundaries(density, params);

  //do a value swap
  for (int is = 0; is < params.ntot; is++)
  {
    for (int iphip = 0; iphip < params.nphip; iphip++)
    {
      for (int ivz = 0; ivz < params.nvz; ivz++)
      {
        density_p[is][iphip][ivz] = density[is][iphip][ivz];
      }
    }
  }

}
*/

void propagateITACollRK4(float ***density, float ***density_i4, float ***density_i3, float ***density_i2, float ***density_p,
                        float *energyDensity, float **flowVelocity, float dt, parameters params)
{
  int ntot = params.ntot;
  int nphip = params.nphip;
  float eta_over_s = params.eta_over_s;
  int nvz = params.nvz;
  //float dvz = 2.0 / (nvz - 1);
  float dvz_2 = 2.0 / (float)(nvz - 1);

  int warn_flag = 1;

  //update the density moment F based on ITA Eqns of Motion
  #pragma omp parallel for
  for (int is = 0; is < ntot; is++)
  {
    float u0 = flowVelocity[0][is];
    float ux = flowVelocity[1][is];
    float uy = flowVelocity[2][is];

    float eps = energyDensity[is];

    float T = temperatureFromEnergyDensity(eps);
    float tau_iso = 5. * eta_over_s / T;

    if ( (tau_iso < 3.0 * dt) && (warn_flag) )
    {
      printf("Warning: tau_iso = %f < 3*dt, energy density = %f , take smaller dt! \n", tau_iso, eps);
      warn_flag = 0;
    }
    for (int iphip = 0; iphip < nphip; iphip++)
    {
      float phip = float(iphip) * (2.0 * M_PI) / float(nphip);
      //float vx = cos(phip);
      //float vy = sin(phip);

      for (int ivz = 0; ivz < nvz; ivz++)
      {
        //float vz = -1.0 + (float)ivz * dvz_2;
        float vz = (nvz > 1) ? -1.0 + (float)ivz * dvz_2 : 0.0;
        float thetap = acos(vz);
        float sin_thetap = (nvz > 1) ? sin(thetap) : 1.0;
        float vx = sin_thetap * cos(phip);
        float vy = sin_thetap * sin(phip);

        //collision term
        float udotv = u0 - ux*vx - uy*vy;
        float F_iso = eps / powf(udotv, 4.0); //the isotropic moment F_iso(x;p),  check factors of 4pi everywhere!!!
        if (nvz == 1) F_iso = eps / powf(udotv, 4.0) * 2.0;

        float F_1 = density_p[is][iphip][ivz];
        float F_2 = density_i2[is][iphip][ivz];
        float F_3 = density_i3[is][iphip][ivz];
        float F_4 = density_i4[is][iphip][ivz];

        float delta_F_1 = F_1 - F_iso;
        float delta_F_2 = F_2 - F_iso;
        float delta_F_3 = F_3 - F_iso;
        float delta_F_4 = F_4 - F_iso;

        float k1 = -1.0 * delta_F_1 * udotv / tau_iso;
        float k2 = -1.0 * delta_F_2 * udotv / tau_iso;
        float k3 = -1.0 * delta_F_3 * udotv / tau_iso;
        float k4 = -1.0 * delta_F_4 * udotv / tau_iso;

        //update the value of F(x; phip)
        float weight_sum_slopes = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        density[is][iphip][ivz] = density_p[is][iphip][ivz] + weight_sum_slopes * dt;

        //if (iphip == 0 && is == (ntot - 1) / 2) std::cout << "k1 * dt = " << k1 * dt << std::endl;
      } //for (int ivz = 0; ivz < nvz; ivz++)
    } //for (int iphip; iphip < nphip; iphip++)
  } //for (int is = 0; is < ntot; is++)
}


// a toy model for the collision term that explicitly conserves energy
void propagateToyColl(float ***density, float ***density_p, parameters params)
{
  int ntot = params.ntot;
  int nphip = params.nphip;
  float eta_over_s = params.eta_over_s;
  int nvz = params.nvz;

  //update the density moment F based on ITA Eqns of Motion
  #pragma omp parallel for
  for (int is = 0; is < ntot; is++)
  {
    for (int iphip = 0; iphip < nphip; iphip++)
    {
      float phip = float(iphip) * (2.0 * M_PI) / float(nphip);

      for (int ivz = 0; ivz < nvz; ivz++)
      {
        float F_rot;
        if (iphip == 0) F_rot = density_p[is][nphip - 1][ivz];
        else F_rot = density_p[is][iphip - 1][ivz];

        //toy collision term
        //rotate momentum dependence by dphi
        density[is][iphip][ivz] = F_rot;

        //if (iphip == 0 && is == (ntot - 1) / 2) std::cout << "k1 * dt = " << k1 * dt << std::endl;
      } //for (int ivz = 0; ivz < nvz; ivz++)
    } //for (int iphip; iphip < nphip; iphip++)
  } //for (int is = 0; is < ntot; is++)
}

//NOTE this will tend to isotropize phi_p in the global/lab frame coordinates, which is not what we want...
//this isotropizes momentum space at each step using the relaxation method (local averaging)
//in this case the parameter alpha is treated as the coupling strength, not the inverse coupling strength
//now only works for case of nvz = 1
//This Kernel conserves energy in the LAB frame (GOOD) but isotropizes LAB frame momentum (BAD)
void propagateRelaxMethodColl(float ***density, float ***density_p, float *energyDensity, float dt, parameters params)
{
  int ntot = params.ntot;
  int nphip = params.nphip;
  float eta_over_s = params.eta_over_s;
  int nvz = params.nvz;

  //update the density moment F based on ITA Eqns of Motion
  #pragma omp parallel for
  for (int is = 0; is < ntot; is++)
  {
    float eps = energyDensity[is];
    float T = temperatureFromEnergyDensity(eps);
    float k = dt * T / (5. * eta_over_s);

    //if ( (1.0 / alpha / T) < 5.0 * dt) std::cout << "take smaller dt! " << std::endl;

    for (int iphip = 0; iphip < nphip; iphip++)
    {
      for (int ivz = 0; ivz < nvz; ivz++)
      {
        float F = density_p[is][iphip][ivz];
        float F_left, F_right;
        //phi_p is periodic
        if (iphip == 0) F_left = density_p[is][nphip - 1][ivz];
        else F_left = density_p[is][iphip - 1][ivz];

        if (iphip == nphip - 1) F_right = density_p[is][0][ivz];
        else F_right = density_p[is][iphip + 1][ivz];

        density[is][iphip][ivz] = F + k * ( F_right + F_left - (2.0 * F) );

      } //for (int ivz = 0; ivz < nvz; ivz++)
    } //for (int iphip; iphip < nphip; iphip++)
  } //for (int is = 0; is < ntot; is++)
}


//this will isotropize (u . v)^4 F in the GLOBAL frame
//doesnt seem to conserve energy
void propagateRelaxMethodColl2(float ***density, float ***density_p, float *energyDensity, float **flowVelocity, float dt, parameters params)
{
  int ntot = params.ntot;
  int nphip = params.nphip;
  float eta_over_s = params.eta_over_s;
  int nvz = params.nvz;
  float dvz_2 = 2.0 / (float)(nvz - 1);

  //update the density moment F based on ITA Eqns of Motion
  #pragma omp parallel for
  for (int is = 0; is < ntot; is++)
  {
    float u0 = flowVelocity[0][is];
    float ux = flowVelocity[1][is];
    float uy = flowVelocity[2][is];

    float eps = energyDensity[is];
    float T = temperatureFromEnergyDensity(eps);
    float k = dt * T / (5. * eta_over_s);

    for (int iphip = 0; iphip < nphip; iphip++)
    {
      float phip = float(iphip) * (2.0 * M_PI) / float(nphip);
      //float vx = cos(phip);
      //float vy = sin(phip);

      for (int ivz = 0; ivz < nvz; ivz++)
      {
        float vz = (nvz > 1) ? -1.0 + (float)ivz * dvz_2 : 0.0;
        float thetap = acos(vz);
        float sin_thetap = (nvz > 1) ? sin(thetap) : 1.0;
        float vx = sin_thetap * cos(phip);
        float vy = sin_thetap * sin(phip);

        float F = density_p[is][iphip][ivz];
        float udotv = u0 - ux*vx - uy*vy;
        float udotv4 = powf(udotv, 4.0);
        float comb = udotv4 * F;

        float F_left, F_right;
        float phip_left, phip_right;

        //phi_p is periodic
        if (iphip == 0)
        {
          F_left = density_p[is][nphip - 1][ivz];
          phip_left = float(nphip - 1) * (2.0 * M_PI) / float(nphip);
        }
        else
        {
          F_left = density_p[is][iphip - 1][ivz];
          phip_left = float(iphip - 1) * (2.0 * M_PI) / float(nphip);
        }

        float vx_left = cos(phip_left);
        float vy_left = sin(phip_left);
        float udotv_left = u0 - ux * vx_left - uy * vy_left;
        float udotv4_left = powf(udotv_left, 4.0);
        float comb_left = udotv4_left * F_left;

        if (iphip == nphip - 1)
        {
          F_right = density_p[is][0][ivz];
          phip_right = float(0) * (2.0 * M_PI) / float(nphip);
        }
        else
        {
          F_right = density_p[is][iphip + 1][ivz];
          phip_right = float(iphip + 1) * (2.0 * M_PI) / float(nphip);
        }

        float vx_right = cos(phip_right);
        float vy_right = sin(phip_right);
        float udotv_right = u0 - ux * vx_right - uy * vy_right;
        float udotv4_right = powf(udotv_right, 4.0);
        float comb_right = udotv4_right * F_right;

        float comb_update = comb + k * (comb_left + comb_right - (2.0 * comb) );
        float F_update = comb_update / udotv4;

        density[is][iphip][ivz] = F_update;

      } //for (int ivz = 0; ivz < nvz; ivz++)
    } //for (int iphip; iphip < nphip; iphip++)
  } //for (int is = 0; is < ntot; is++)
}


//this uses the exact solution of (dF/dt) = (dF/dt)|coll to propagate F forward by one time step
void propagateITACollExact(float ***density, float ***density_p, float *energyDensity, float **flowVelocity, float dt, parameters params)
{
  int ntot = params.ntot;
  int nphip = params.nphip;
  float eta_over_s = params.eta_over_s;
  int nvz = params.nvz;
  float dvz_2 = 2.0 / (float)(nvz - 1);

  int warn_flag = 1;

  #pragma omp parallel for
  for (int is = 0; is < ntot; is++)
  {
    float u0 = flowVelocity[0][is];
    float ux = flowVelocity[1][is];
    float uy = flowVelocity[2][is];

    float eps = energyDensity[is];

    float T = temperatureFromEnergyDensity(eps);
    float tau_iso = 5. * eta_over_s / T;

    if ( (tau_iso < 3.0 * dt) && (warn_flag) )
    {
      printf("Warning: tau_iso = %f < 3*dt, energy density = %f , take smaller dt! \n", tau_iso, eps);
      warn_flag = 0;
    }
    for (int iphip = 0; iphip < nphip; iphip++)
    {
      float phip = float(iphip) * (2.0 * M_PI) / float(nphip);

      for (int ivz = 0; ivz < nvz; ivz++)
      {
        float vz = (nvz > 1) ? -1.0 + (float)ivz * dvz_2 : 0.0;
        float thetap = acos(vz);
        float sin_thetap = (nvz > 1) ? sin(thetap) : 1.0;
        float vx = sin_thetap * cos(phip);
        float vy = sin_thetap * sin(phip);

        float F = density_p[is][iphip][ivz];

        //collision term
        float udotv = u0 - ux*vx - uy*vy;
        float F_iso = eps / powf(udotv, 4.0); //the isotropic moment F_iso(x;p),  check factors of 4pi everywhere!!!
        //for special case when we assume dist function ~ delta(v_z)
        if (nvz == 1) F_iso = eps / powf(udotv, 4.0) * 2.0;
        float nu = udotv / tau_iso;

        if (nu * dt > 0.2) printf("Warning: nu * dt = %f > 0.2, energy density = %f , take smaller dt! \n", nu * dt, eps);
        float exp_weight = exp(-1.0 * nu * dt);

        //std::cout << "exp_weight = " << exp_weight << std::endl;

        //update the value of F(x; phip)
        density[is][iphip][ivz] = exp_weight * F + (1.0 - exp_weight) * F_iso;

        //if (iphip == 0 && is == (ntot - 1) / 2) std::cout << "k1 * dt = " << k1 * dt << std::endl;
      } //for (int ivz = 0; ivz < nvz; ivz++)
    } //for (int iphip; iphip < nphip; iphip++)
  } //for (int is = 0; is < ntot; is++)
}

//this boosts F into the LRF
// NOTE THIS IS BROKEN; vx, vy need sin_thetap factors!!
void propagateRelaxMethodCollLRF(float ***density, float ***density_p, float *energyDensity, float **flowVelocity, float dt, parameters params)
{
  int ntot = params.ntot;
  int nphip = params.nphip;
  float eta_over_s = params.eta_over_s;
  int nvz = params.nvz;
  float delta_phip = (2.0 * M_PI) / float(nphip);

  //update the density moment F based on ITA Eqns of Motion
  #pragma omp parallel for
  for (int is = 0; is < ntot; is++)
  {
    float u0 = flowVelocity[0][is];
    float ux = flowVelocity[1][is];
    float uy = flowVelocity[2][is];

    //TEMPORARY
    //check that when boost is trivial (LRF = LAB) that energy is conserved
    //u0 = 1.0;
    //ux = 0.0;
    //uy = 0.0;

    //TEMPORARY

    float gamma = u0;
    float beta_x = ux / u0;
    float beta_y = uy / u0;
    float beta2 = beta_x * beta_x + beta_y * beta_y;
    beta2 = beta2 + 1.0e-15; //this ensures that quantities are numerically defined when beta -> 0

    float eps = energyDensity[is]; //in LRF
    float T = temperatureFromEnergyDensity(eps); //in LRF
    float k = dt * T / (5. * eta_over_s); //should we lorentz dilate dt ???
    k = k / gamma; //time dilation

    //if ( (1.0 / alpha / T) < 5.0 * dt) std::cout << "take smaller dt! " << std::endl;

    //store the interpolation of F as a function of phi_p
    double phip_arr[nphip + 1];
    double F_phip_arr[nphip + 1];

    for (int iphip = 0; iphip < nphip + 1; iphip++)
    {
      if (iphip == nphip) F_phip_arr[iphip] = (double)density_p[is][0][0];
      else F_phip_arr[iphip] = (double)density_p[is][iphip][0];
      phip_arr[iphip] = double(iphip) * delta_phip;
    }

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *F_spline_phip = gsl_spline_alloc(gsl_interp_cspline, nphip + 1);
    gsl_spline_init(F_spline_phip, phip_arr, F_phip_arr, nphip + 1); //an interpolation of F(phi_p) in the LAB

    for (int iphip = 0; iphip < nphip; iphip++)
    {
      //isotropic difference in LRF
      float phip_LRF = phip_arr[iphip];
      float phip_LRF_l = phip_LRF - delta_phip;
      float phip_LRF_r = phip_LRF + delta_phip;

      if (phip_LRF < 0) phip_LRF += 2.0 * M_PI;
      if (phip_LRF_l < 0) phip_LRF_l += 2.0 * M_PI;
      if (phip_LRF_r < 0) phip_LRF_r += 2.0 * M_PI;

      if (phip_LRF > 2.0 * M_PI) phip_LRF -= 2.0 * M_PI;
      if (phip_LRF_l > 2.0 * M_PI) phip_LRF_l -= 2.0 * M_PI;
      if (phip_LRF_r > 2.0 * M_PI) phip_LRF_r -= 2.0 * M_PI;

      //v^mu in LRF, useful for computing boost formula
      float vx_LRF = cos(phip_LRF);
      float vy_LRF = sin(phip_LRF);

      float vx_LRF_l = cos(phip_LRF_l);
      float vy_LRF_l = sin(phip_LRF_l);

      float vx_LRF_r = cos(phip_LRF_r);
      float vy_LRF_r = sin(phip_LRF_r);

      //need the inverse transformation
      //this is three vector product beta^i v^i
      float beta_dot_v_LRF   = beta_x * vx_LRF + beta_y * vy_LRF;
      float beta_dot_v_LRF_l = beta_x * vx_LRF_l + beta_y * vy_LRF_l;
      float beta_dot_v_LRF_r = beta_x * vx_LRF_r + beta_y * vy_LRF_r;

      float num_x_LRF =  vx_LRF + ( (gamma - 1.0) / beta2) * beta_dot_v_LRF * beta_x - gamma * beta_x;
      float num_y_LRF =  vy_LRF + ( (gamma - 1.0) / beta2) * beta_dot_v_LRF * beta_y - gamma * beta_y;
      float den_LRF = gamma * (1.0 + beta_dot_v_LRF_l);

      float num_x_LRF_l =  vx_LRF_l + ( (gamma - 1.0) / beta2) * beta_dot_v_LRF_l * beta_x - gamma * beta_x;
      float num_y_LRF_l =  vy_LRF_l + ( (gamma - 1.0) / beta2) * beta_dot_v_LRF_l * beta_y - gamma * beta_y;
      float den_LRF_l = gamma * (1.0 + beta_dot_v_LRF_l);

      float num_x_LRF_r =  vx_LRF_r + ( (gamma - 1.0) / beta2) * beta_dot_v_LRF_r * beta_x - gamma * beta_x;
      float num_y_LRF_r =  vy_LRF_r + ( (gamma - 1.0) / beta2) * beta_dot_v_LRF_r * beta_y - gamma * beta_y;
      float den_LRF_r = gamma * (1.0 + beta_dot_v_LRF_r);

      //these angles phip are in the LAB frame
      float vx = num_x_LRF / den_LRF;
      float vy = num_y_LRF / den_LRF;
      float phip = atan2(vy, vx);

      float vx_l = num_x_LRF_l / den_LRF_l;
      float vy_l = num_y_LRF_l / den_LRF_l;
      float phip_l = atan2(vy_l, vx_l);

      float vx_r = num_x_LRF_r / den_LRF_r;
      float vy_r = num_y_LRF_r / den_LRF_r;
      float phip_r = atan2(vy_r, vx_r);

      if (phip < 0) phip += 2.0 * M_PI;
      if (phip_l < 0) phip_l += 2.0 * M_PI;
      if (phip_r < 0) phip_r += 2.0 * M_PI;

      if (phip > 2.0 * M_PI) phip -= 2.0 * M_PI;
      if (phip_l > 2.0 * M_PI) phip_l -= 2.0 * M_PI;
      if (phip_r > 2.0 * M_PI) phip_r -= 2.0 * M_PI;

      //std::cout << "num_x_LRF = " << num_x_LRF << " , num_y_LRF = " << num_y_LRF << "\n";
      //std::cout << "den_LRF = " << den_LRF << "\n";
      //std::cout << "phip_LRF = " << phip_LRF << " , phip = " << phip << "\n";

      for (int ivz = 0; ivz < nvz; ivz++)
      {
        //now transform F to the LRF
        float boost_fac = powf( gamma * (1.0 + beta_dot_v_LRF) , 4.0);
        float boost_fac_l = powf( gamma * (1.0 + beta_dot_v_LRF_l) , 4.0);
        float boost_fac_r = powf( gamma * (1.0 + beta_dot_v_LRF_r) , 4.0);

        float F_LRF = boost_fac * gsl_spline_eval(F_spline_phip, phip, acc);
        float F_LRF_l = boost_fac_l * gsl_spline_eval(F_spline_phip, phip_l, acc);
        float F_LRF_r = boost_fac_r * gsl_spline_eval(F_spline_phip, phip_r, acc);

        float F_LRF_update = F_LRF + k * ( F_LRF_r + F_LRF_l - (2.0 * F_LRF) );

        //now boost F_LRF back to LAB frame
        float F_update = F_LRF_update / boost_fac;

        density[is][iphip][ivz] = F_update;

      } //for (int ivz = 0; ivz < nvz; ivz++)
    } //for (int iphip; iphip < nphip; iphip++)

    gsl_spline_free(F_spline_phip);
    gsl_interp_accel_free(acc);
  } //for (int is = 0; is < ntot; is++)
}
