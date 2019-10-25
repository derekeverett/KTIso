//trigTable is a table with 10 rows for ten combinations or p^(mu)p_(nu) normalized by the energy
#pragma once
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
//#include <math.h>
#include "Parameter.h"
#include <iostream>

//this constructs the tensor v^mu v^nu
void calculateHypertrigTable(float ***hypertrigTable, float **vz_quad, parameters params)
{
  int nphip = params.nphip;
  int nvz = params.nvz;
  //float dvz = 2.0 / (float)(nvz);
  float dvz_2 = 2.0 / (float)(nvz - 1);

  #pragma omp parallel for
  for (int iphip = 0; iphip < nphip; iphip++)
  {
    float phip = float(iphip) * (2.0 * M_PI) / float(nphip);
    for (int ivz = 0; ivz < nvz; ivz++)
    {
      float vz = -1.0 + (float)ivz * dvz_2;
      //float vz = vz_quad[ivz][0];
      if (nvz == 1) vz = 0.0;

      float thetap = acos(vz);
      float sin_thetap = sin(thetap);
      float vx = sin_thetap * cos(phip);
      float vy = sin_thetap * sin(phip);

      /*
      hypertrigTable[0][iphip][ivz] = 1.0; //p^t, p^t component
      hypertrigTable[1][iphip][ivz] = cos(phip); //p^t, p^x
      hypertrigTable[2][iphip][ivz] = sin(phip); //p^t, p^y
      hypertrigTable[3][iphip][ivz] = vz; //p^t, p^z
      hypertrigTable[4][iphip][ivz] = cos(phip) * cos(phip); //p^x, p^x
      hypertrigTable[5][iphip][ivz] = cos(phip) * sin(phip); //p^x, p^y
      hypertrigTable[6][iphip][ivz] = cos(phip) * vz; //p^x, p^z
      hypertrigTable[7][iphip][ivz] = sin(phip) * sin(phip); //p^y, p^y
      hypertrigTable[8][iphip][ivz] = sin(phip) * vz; //p^y, p^z
      hypertrigTable[9][iphip][ivz] = vz * vz; //p^z, p^z
      */

      hypertrigTable[0][iphip][ivz] = 1.0; //p^t, p^t component
      hypertrigTable[1][iphip][ivz] = vx; //p^t, p^x
      hypertrigTable[2][iphip][ivz] = vy; //p^t, p^y
      hypertrigTable[3][iphip][ivz] = vz; //p^t, p^z
      hypertrigTable[4][iphip][ivz] = vx * vx; //p^x, p^x
      hypertrigTable[5][iphip][ivz] = vx * vy; //p^x, p^y
      hypertrigTable[6][iphip][ivz] = vx * vz; //p^x, p^z
      hypertrigTable[7][iphip][ivz] = vy * vy; //p^y, p^y
      hypertrigTable[8][iphip][ivz] = vy * vz; //p^y, p^z
      hypertrigTable[9][iphip][ivz] = vz * vz; //p^z, p^z
    }

  }
}

void calculateStressTensor(float **stressTensor, float ***density, float ***hypertrigTable, float **vz_quad, float t, parameters params)
{
  int nphip = params.nphip;
  int ntot = params.ntot;
  float d_phip = (2.0 * M_PI) / float(nphip);
  int nvz = params.nvz;
  float dvz = 2.0 / (float)(nvz);

  //Right now just a Riemann Sum over phip and trapezoid rule for vz
  for (int ivar = 0; ivar < 10; ivar++)
  {
    #pragma omp parallel for
    for (int is = 0; is < ntot; is++) //the column packed index for x, y and z
    {

      //case of F ~ delta(v_z), boost invariant
      if (nvz == 1)
      {
        float integral = 0.0;
        for (int iphip = 0; iphip < nphip; iphip++)
        {
          integral += density[is][iphip][0] * hypertrigTable[ivar][iphip][0];
        } // for (int iphip = 0; iphip < nphip; iphip++)
        //stressTensor[ivar][is] = integral * (d_phip / (2.0 * M_PI) ) / t; //should we divide by tau???
        stressTensor[ivar][is] = integral * (d_phip / (4.0 * M_PI) );
      } // if (nvz == 1)

      else
      {
        float integral = 0.0;
        for (int iphip = 0; iphip < nphip; iphip++)
        {
          for (int ivz = 0; ivz < nvz; ivz++)
          {
            if (density[is][iphip][ivz] < 0.0) printf("Warning : F < 0 \n");

            float vz_weight = 1.0;
            //float vz_weight = vz_quad[ivz][1];
            if ( (ivz == 0) || (ivz == nvz - 1) ) vz_weight = 0.5;
            integral += density[is][iphip][ivz] * hypertrigTable[ivar][iphip][ivz] * vz_weight;
          }
        }
        stressTensor[ivar][is] = integral * dvz * (d_phip / (2.0 * M_PI) );
        //vz only runs from 0 to 1, integration from -1 to 1 is even, so multiply by 2
        //stressTensor[ivar][is] = integral * (d_phip / (2.0 * M_PI) ) * 2.0;

      } //else
    } //for (int is = 0; is < ntot; is++)
  } // for (int ivar = 0; ivar < 10; ivar++)
}

void solveEigenSystem(float **stressTensor, float *energyDensity, float **flowVelocity, parameters params)
{
  int ntot = params.ntot;

  //float tolerance = 1.0e-5;
  float gamma_max = 100.0;

  #pragma omp parallel for
  for (int is = 0; is < ntot; is++)
  {
    gsl_matrix *Tmunu; //T^(mu,nu) with two contravariant indices; we need to lower an index
    //using the metric to find the eigenvectors of T^(mu)_(nu) with one contravariant and one contravariant index
    Tmunu = gsl_matrix_alloc(4,4);
    gsl_matrix *gmunu;
    gmunu = gsl_matrix_alloc(4,4);
    gsl_matrix_complex *eigen_vectors;
    eigen_vectors = gsl_matrix_complex_alloc(4,4);
    gsl_vector_complex *eigen_values;
    eigen_values = gsl_vector_complex_alloc(4);

    //set the values of the energy momentum tensor
    //try adding a small value everywhere to T^tt make flow velocity look nicer
    gsl_matrix_set(Tmunu, 0, 0, stressTensor[0][is]); //t,t
    gsl_matrix_set(Tmunu, 0, 1, stressTensor[1][is]); //t,x
    gsl_matrix_set(Tmunu, 0, 2, stressTensor[2][is]); //t,y
    gsl_matrix_set(Tmunu, 0, 3, stressTensor[3][is]); //t,z
    gsl_matrix_set(Tmunu, 1, 1, stressTensor[4][is]); //x,x
    gsl_matrix_set(Tmunu, 1, 2, stressTensor[5][is]); //x,y
    gsl_matrix_set(Tmunu, 1, 3, stressTensor[6][is]); //x,z
    gsl_matrix_set(Tmunu, 2, 2, stressTensor[7][is]); //y,y
    gsl_matrix_set(Tmunu, 2, 3, stressTensor[8][is]); //y,z
    gsl_matrix_set(Tmunu, 3, 3, stressTensor[9][is]); //z,z
    gsl_matrix_set(Tmunu, 1, 0, stressTensor[1][is]); //x,t
    gsl_matrix_set(Tmunu, 2, 0, stressTensor[2][is]); //y,t
    gsl_matrix_set(Tmunu, 3, 0, stressTensor[3][is]); //z,t
    gsl_matrix_set(Tmunu, 2, 1, stressTensor[5][is]); //y,x
    gsl_matrix_set(Tmunu, 3, 1, stressTensor[6][is]); //z,x
    gsl_matrix_set(Tmunu, 3, 2, stressTensor[8][is]); //z,y

    //set the values of the "metric"; not really the metric, but the numerical constants
    //which are multiplied by the elements of T^(mu,nu) to get the values of T^(mu)_(nu)
    //g_(mu.nu) = diag(1,-1,-1,-1)
    gsl_matrix_set(gmunu, 0, 0, 1.0); //t,t
    gsl_matrix_set(gmunu, 0, 1, -1.0); //t,x
    gsl_matrix_set(gmunu, 0, 2, -1.0); //t,y
    gsl_matrix_set(gmunu, 0, 3, -1.0); //t,z
    gsl_matrix_set(gmunu, 1, 0, 1.0); //x,t
    gsl_matrix_set(gmunu, 1, 1, -1.0); //x,x
    gsl_matrix_set(gmunu, 1, 2, -1.0); //x,y
    gsl_matrix_set(gmunu, 1, 3, -1.0); //x,z
    gsl_matrix_set(gmunu, 2, 0, 1.0); //y,t
    gsl_matrix_set(gmunu, 2, 1, -1.0); //y,x
    gsl_matrix_set(gmunu, 2, 2, -1.0); //y,y
    gsl_matrix_set(gmunu, 2, 3, -1.0); //y,z
    gsl_matrix_set(gmunu, 3, 0, 1.0); //z,t
    gsl_matrix_set(gmunu, 3, 1, -1.0); //z,x
    gsl_matrix_set(gmunu, 3, 2, -1.0); //z,y
    gsl_matrix_set(gmunu, 3, 3, -1.0); //z,z
    //lower one index of the stress tensor; save it to the same matrix to save memory
    gsl_matrix_mul_elements(Tmunu, gmunu); //result stored in Tmunu !this multiplies element-wise, not ordinary matrix multiplication!
    gsl_eigen_nonsymmv_workspace *eigen_workspace;
    eigen_workspace = gsl_eigen_nonsymmv_alloc(4);
    gsl_eigen_nonsymmv(Tmunu, eigen_values, eigen_vectors, eigen_workspace);
    gsl_eigen_nonsymmv_free(eigen_workspace);

    int eigenvalue_exists = 0;
    for (int i = 0; i < 4; i++)
    {
      gsl_complex eigenvalue = gsl_vector_complex_get(eigen_values, i);
      //eigenvalue condition taken from JF's suggestion, test for robustness in dilute region
      if ( GSL_REAL(eigenvalue) > 0.0 && fabs( GSL_IMAG(eigenvalue) ) < ( fabs(GSL_REAL(eigenvalue)) * 1.0e-30) ) //choose eigenvalue
      {
        gsl_complex v0 = gsl_matrix_complex_get(eigen_vectors, 0 , i);
        gsl_complex v1 = gsl_matrix_complex_get(eigen_vectors, 1 , i);
        gsl_complex v2 = gsl_matrix_complex_get(eigen_vectors, 2 , i);
        gsl_complex v3 = gsl_matrix_complex_get(eigen_vectors, 3 , i);

        if ( GSL_IMAG(v0) == 0 && ( 2.0 * GSL_REAL(v0) * GSL_REAL(v0) - 1.0  > 0.0 ) ) //choose timelike eigenvector
        //if ( 2.0 * GSL_REAL(v0) * GSL_REAL(v0) - 1.0  > 0.0 ) //choose timelike eigenvector
        {
          double minkowskiLength = GSL_REAL(v0)*GSL_REAL(v0) - (GSL_REAL(v1)*GSL_REAL(v1) + GSL_REAL(v2)*GSL_REAL(v2) + GSL_REAL(v3)*GSL_REAL(v3));
          double factor = 1.0 / sqrt(minkowskiLength);
          if (GSL_REAL(v0) < 0) factor=-factor;

          //ignore eigenvectors with gamma too large
          if ( (GSL_REAL(v0) * factor) < gamma_max)
          {
            eigenvalue_exists = 1;
            energyDensity[is] = GSL_REAL(eigenvalue);
            flowVelocity[0][is] = GSL_REAL(v0) * factor;
            flowVelocity[1][is] = GSL_REAL(v1) * factor;
            flowVelocity[2][is] = GSL_REAL(v2) * factor;
            flowVelocity[3][is] = GSL_REAL(v3) * factor;
          }
        } // if (GSL_IMAG(v0) == 0 && (2.0 * GSL_REAL(v0) * GSL_REAL(v0) - 1.0 - (GSL_REAL(v3) * GSL_REAL(v3) * (TAU * TAU - 1.0) )) > 0) //choose timelike eigenvector
      } // if (GSL_REAL(eigenvalue) > 0.0 && GSL_IMAG(eigenvalue) == 0) //choose eigenvalue
    } //for (int i = 0; i < 4; ...)

    if (eigenvalue_exists == 0)
    {
      //in dilute regions where we can't find a timelike eigenvector, set e = 0, u^t = 1, u^x=u^y=u^n=0
      energyDensity[is] = 1.0e-7;
      flowVelocity[0][is] = 1.0;
      flowVelocity[1][is] = 0.0;
      flowVelocity[2][is] = 0.0;
      flowVelocity[3][is] = 0.0;
    }
  } // for (int is; is < ntot; ...)
} //solveEigenSystem()

void calculateBulkPressure(float **stressTensor, float *energyDensity, float *pressure, float *bulkPressure, parameters params)
{
  int ntot = params.ntot;
  #pragma omp parallel for
  for (int is = 0; is < ntot; is++)
  {
    // PI = -1/3 * (T^(mu)_(mu) - epsilon) - p
    // T^(mu)_(mu) = T^(0,0) - T^(1,1) - T^(2,2) - (TAU^2)T^(3,3)
    float a =  stressTensor[0][is] - stressTensor[4][is] - stressTensor[7][is] - stressTensor[9][is];
    bulkPressure[is] = (-1.0/3.0) * (a - energyDensity[is]) - pressure[is];
  }
}
void calculateShearViscTensor(float **stressTensor, float *energyDensity, float **flowVelocity, float *pressure, float *bulkPressure, float **shearTensor, parameters params)
{
  int ntot = params.ntot;
  #pragma omp parallel for
  for (int is = 0; is < ntot; is++)
  {
    // pi^(mu,nu) = T^(mu,nu) - epsilon * u^(mu)u^(nu) + (P + PI) * (g^(mu,nu) - u^(mu)u^(nu))
    //calculate ten components : upper triangular part
    float b = energyDensity[is] + pressure[is] + bulkPressure[is];
    float c = pressure[is] + bulkPressure[is];
    shearTensor[0][is] = stressTensor[0][is] - flowVelocity[0][is] * flowVelocity[0][is] * b + c; //pi^(t,t)
    shearTensor[1][is] = stressTensor[1][is] - flowVelocity[0][is] * flowVelocity[1][is] * b; //pi^(t,x)
    shearTensor[2][is] = stressTensor[2][is] - flowVelocity[0][is] * flowVelocity[2][is] * b; //pi^(t,y)
    shearTensor[3][is] = stressTensor[3][is] - flowVelocity[0][is] * flowVelocity[3][is] * b; //pi^(t,z)
    shearTensor[4][is] = stressTensor[4][is] - flowVelocity[1][is] * flowVelocity[1][is] * b - c; //pi^(x,x)
    shearTensor[5][is] = stressTensor[5][is] - flowVelocity[1][is] * flowVelocity[2][is] * b; //pi^(x,y)
    shearTensor[6][is] = stressTensor[6][is] - flowVelocity[1][is] * flowVelocity[3][is] * b; //pi^(x,z)
    shearTensor[7][is] = stressTensor[7][is] - flowVelocity[2][is] * flowVelocity[2][is] * b - c; //pi^(y,y)
    shearTensor[8][is] = stressTensor[8][is] - flowVelocity[2][is] * flowVelocity[3][is] * b; //pi^(y,z)
    shearTensor[9][is] = stressTensor[9][is] - flowVelocity[3][is] * flowVelocity[3][is] * b - c ; //pi^(z,z)
  }
}

/*
calculateF_iso(float ***F_iso, float ***density, float *energyDensity, float **flowVelocity, parameters params)
{
  int ntot = params.ntot;
  int nphip = params.nphip;
  float alpha = params.ALPHA;
  int nvz = params.nvz;
  //float dvz = 2.0 / (nvz - 1);

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
    float tau_iso = alpha / T;

    if ( (tau_iso < 3.0 * dt) && (warn_flag) )
    {
      printf("Warning: tau_iso = %f < 3*dt, energy density = %f , take smaller dt! \n", tau_iso, eps);
      warn_flag = 0;
    }
    for (int iphip = 0; iphip < nphip; iphip++)
    {
      float phip = float(iphip) * (2.0 * M_PI) / float(nphip);
      float vx = cos(phip);
      float vy = sin(phip);

      for (int ivz = 0; ivz < nvz; ivz++)
      {
        float F = density[is][iphip][ivz];

        //collision term
        float udotv = u0 - ux*vx - uy*vy;
        float F_iso = eps / powf(udotv, 4.0); //the isotropic moment F_iso(x;p),  check factors of 4pi everywhere!!!

        //if (iphip == 0 && is == (ntot - 1) / 2) std::cout << "k1 * dt = " << k1 * dt << std::endl;
      } //for (int ivz = 0; ivz < nvz; ivz++)
    } //for (int iphip; iphip < nphip; iphip++)
  } //for (int is = 0; is < ntot; is++)

}
*/

float calculateLongitudinalWork(float **stressTensor, float t, float dt, parameters params)
{
  int ntot = params.ntot;
  float dx = params.dx;
  float dy = params.dy;

  float work = 0.0;
  for (int is = 0; is < ntot; is++)
  {
    float T_zz = stressTensor[9][is];
    work += T_zz;
  }

  //what is the right factor? we are propagating stress tensor cartesian components right?
  //float factor = t*t*dt*dx*dy;
  float factor = dt*dx*dy;
  work *= factor;

  return work;

}
