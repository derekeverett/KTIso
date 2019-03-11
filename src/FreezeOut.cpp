#include "cornelius-c++-1.3/cornelius.cpp"
#include <stdio.h>
#include <iostream>
#include <fstream>

//return a 4 dimensional linear interpolation inside the hypercube, given the values
//on the corners (a0000 through a1111) and edge lengths x0 through x3
float linearInterp4D(float x0, float x1, float x2, float x3,
                      float a0000, float a1000, float a0100, float a0010, float a0001,
                      float a1100, float a1010, float a1001,
                      float a0110, float a0101, float a0011,
                      float a1110, float a1101, float a0111, float a1011, float a1111)
{
  float result = 0;
  result = ((1-x0) * (1-x1) * (1-x2) * (1-x3) * a0000)
            + ((x0) * (1-x1) * (1-x2) * (1-x3) * a1000)
            + ((1-x0) * (x1) * (1-x2) * (1-x3) * a0100)
            + ((1-x0) * (1-x1) * (x2) * (1-x3) * a0010)
            + ((1-x0) * (1-x1) * (1-x2) * (x3) * a0001)
            + ((x0) * (x1) * (1-x2) * (1-x3) * a1100)
            + ((x0) * (1-x1) * (x2) * (1-x3) * a1010)
            + ((x0) * (1-x1) * (1-x2) * (x3) * a1001)
            + ((1-x0) * (x1) * (x2) * (1-x3) * a0110)
            + ((1-x0) * (x1) * (1-x2) * (x3) * a0101)
            + ((1-x0) * (1-x1) * (x2) * (x3) * a0011)
            + ((x0) * (x1) * (x2) * (1-x3) * a1110)
            + ((x0) * (x1) * (1-x2) * (x3) * a1101)
            + ((x0) * (1-x1) * (x2) * (x3) * a1011)
            + ((1-x0) * (x1) * (x2) * (x3) * a0111)
            + ((x0) * (x1) * (x2) * (x3) * a1111);

  return result;
}

float linearInterp3D(float x0, float x1, float x2,
                      float a000, float a100, float a010, float a001,
                      float a110, float a101, float a011, float a111)
{
  float result = 0;
  result = ((1-x0) * (1-x1) * (1-x2) * a000)
            + ((x0) * (1-x1) * (1-x2) * a100)
            + ((1-x0) * (x1) * (1-x2) * a010)
            + ((1-x0) * (1-x1) * (x2) * a001)
            + ((x0) * (x1) * (1-x2) * a110)
            + ((x0) * (1-x1) * (x2) * a101)
            + ((1-x0) * (x1) * (x2) * a011)
            + ((x0) * (x1) * (x2)  * a111);

  return result;
}

void swapAndSetHydroVariables(float ****energy_density_evoution, float *****hydrodynamic_evoution,
                              CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e,
                              FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz)
{
  #pragma omp parallel for collapse(3)
  for (int ix = 0; ix < nx; ix++)
  {
    for (int iy = 0; iy < ny; iy++)
    {
      for (int iz = 0; iz < nz; iz++)
      {
        int s = columnMajorLinearIndex(ix, iy, iz, nx, ny, nz);
        //previous hydro variable values written to zeroth index
        energy_density_evoution[0][ix-2][iy-2][iz-2] = energy_density_evoution[1][ix-2][iy-2][iz-2];

        for (int ivar = 0; ivar < 10; ivar++)
        {
          hydrodynamic_evoution[ivar][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[ivar][1][ix-2][iy-2][iz-2];
        }

        //current hydro variable values written to first index
        energy_density_evoution[1][ix-2][iy-2][iz-2] = e[s];

        hydrodynamic_evoution[0][1][ix-2][iy-2][iz-2] = ux[s];
        hydrodynamic_evoution[1][1][ix-2][iy-2][iz-2] = uy[s];
        hydrodynamic_evoution[3][1][ix-2][iy-2][iz-2] = e[s];
        hydrodynamic_evoution[4][1][ix-2][iy-2][iz-2] = pixx[s];
        hydrodynamic_evoution[5][1][ix-2][iy-2][iz-2] = pixy[s];
        hydrodynamic_evoution[6][1][ix-2][iy-2][iz-2] = pixn[s];
        hydrodynamic_evoution[7][1][ix-2][iy-2][iz-2] = piyy[s];
        hydrodynamic_evoution[8][1][ix-2][iy-2][iz-2] = piyn[s];
      }
    } //for (int iy = 2; iy < ny+2; iy++)
  } //for (int ix = 2; ix < nx+2; ix++)
}

//freezeout functions for vorticity and polarization studies

void swapAndSetHydroVariables_Vorticity(float ****energy_density_evoution, float *****hydrodynamic_evoution,
                              CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e,
                              FLUID_VELOCITY * const __restrict__ u, VORTICITY * const __restrict__ wmunu,
                              int nx, int ny, int nz)
{
  #pragma omp parallel for collapse(3)
  for (int ix = 2; ix < nx+2; ix++)
  {
    for (int iy = 2; iy < ny+2; iy++)
    {
      for (int iz = 2; iz < nz+2; iz++)
      {
        int s = columnMajorLinearIndex(ix, iy, iz, nx+4, ny+4, nz+4);
        //previous hydro variable values written to zeroth index
        energy_density_evoution[0][ix-2][iy-2][iz-2] = energy_density_evoution[1][ix-2][iy-2][iz-2];
        for (int ivar = 0; ivar < 16; ivar++)
        {
          hydrodynamic_evoution[ivar][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[ivar][1][ix-2][iy-2][iz-2];
        }

        //current hydro variable values written to first index
        energy_density_evoution[1][ix-2][iy-2][iz-2] = (float)e[s];

        hydrodynamic_evoution[0][1][ix-2][iy-2][iz-2] = (float)(u->ux[s]);
        hydrodynamic_evoution[1][1][ix-2][iy-2][iz-2] = (float)(u->uy[s]);
        hydrodynamic_evoution[2][1][ix-2][iy-2][iz-2] = (float)(u->un[s]);
        hydrodynamic_evoution[3][1][ix-2][iy-2][iz-2] = (float)(e[s]);
	#ifdef PIMUNU
        hydrodynamic_evoution[4][1][ix-2][iy-2][iz-2] = (float)(q->pixx[s]);
        hydrodynamic_evoution[5][1][ix-2][iy-2][iz-2] = (float)(q->pixy[s]);
        hydrodynamic_evoution[6][1][ix-2][iy-2][iz-2] = (float)(q->pixn[s]);
        hydrodynamic_evoution[7][1][ix-2][iy-2][iz-2] = (float)(q->piyy[s]);
        hydrodynamic_evoution[8][1][ix-2][iy-2][iz-2] = (float)(q->piyn[s]);
	#endif
	#ifdef PI
        hydrodynamic_evoution[9][1][ix-2][iy-2][iz-2] = (float)(q->Pi[s]);
	#endif
  #ifdef THERMAL_VORTICITY
        hydrodynamic_evoution[10][1][ix-2][iy-2][iz-2] = (float)(wmunu->wtx[s]);
        hydrodynamic_evoution[11][1][ix-2][iy-2][iz-2] = (float)(wmunu->wty[s]);
        hydrodynamic_evoution[12][1][ix-2][iy-2][iz-2] = (float)(wmunu->wtn[s]);
        hydrodynamic_evoution[13][1][ix-2][iy-2][iz-2] = (float)(wmunu->wxy[s]);
        hydrodynamic_evoution[14][1][ix-2][iy-2][iz-2] = (float)(wmunu->wxn[s]);
        hydrodynamic_evoution[15][1][ix-2][iy-2][iz-2] = (float)(wmunu->wyn[s]);
  #endif
      }
    } //for (int iy = 2; iy < ny+2; iy++)
  } //for (int ix = 2; ix < nx+2; ix++)
}


void writeEnergyDensityToHypercube4D(float ****hyperCube, float ****energy_density_evoution, int it, int ix, int iy, int iz)
{
  hyperCube[0][0][0][0] = energy_density_evoution[it][ix][iy][iz];
  hyperCube[1][0][0][0] = energy_density_evoution[it+1][ix][iy][iz];
  hyperCube[0][1][0][0] = energy_density_evoution[it][ix+1][iy][iz];
  hyperCube[0][0][1][0] = energy_density_evoution[it][ix][iy+1][iz];
  hyperCube[0][0][0][1] = energy_density_evoution[it][ix][iy][iz+1];
  hyperCube[1][1][0][0] = energy_density_evoution[it+1][ix+1][iy][iz];
  hyperCube[1][0][1][0] = energy_density_evoution[it+1][ix][iy+1][iz];
  hyperCube[1][0][0][1] = energy_density_evoution[it+1][ix][iy][iz+1];
  hyperCube[0][1][1][0] = energy_density_evoution[it][ix+1][iy+1][iz];
  hyperCube[0][1][0][1] = energy_density_evoution[it][ix+1][iy][iz+1];
  hyperCube[0][0][1][1] = energy_density_evoution[it][ix][iy+1][iz+1];
  hyperCube[1][1][1][0] = energy_density_evoution[it+1][ix+1][iy+1][iz];
  hyperCube[1][1][0][1] = energy_density_evoution[it+1][ix+1][iy][iz+1];
  hyperCube[1][0][1][1] = energy_density_evoution[it+1][ix][iy+1][iz+1];
  hyperCube[0][1][1][1] = energy_density_evoution[it][ix+1][iy+1][iz+1];
  hyperCube[1][1][1][1] = energy_density_evoution[it+1][ix+1][iy+1][iz+1];
}
void writeEnergyDensityToHypercube3D(float ***hyperCube, float ****energy_density_evoution, int it, int ix, int iy)
{
  hyperCube[0][0][0] = energy_density_evoution[it][ix][iy][0];
  hyperCube[1][0][0] = energy_density_evoution[it+1][ix][iy][0];
  hyperCube[0][1][0] = energy_density_evoution[it][ix+1][iy][0];
  hyperCube[0][0][1] = energy_density_evoution[it][ix][iy+1][0];
  hyperCube[1][1][0] = energy_density_evoution[it+1][ix+1][iy][0];
  hyperCube[1][0][1] = energy_density_evoution[it+1][ix][iy+1][0];
  hyperCube[0][1][1] = energy_density_evoution[it][ix+1][iy+1][0];
  hyperCube[1][1][1] = energy_density_evoution[it+1][ix+1][iy+1][0];
}
float interpolateVariable4D(float *****hydrodynamic_evoution, int ivar, int it, int ix, int iy, int iz, float tau_frac, float x_frac, float y_frac, float z_frac)
{
  float result = linearInterp4D(tau_frac, x_frac, y_frac, z_frac,
    hydrodynamic_evoution[ivar][it][ix][iy][iz], hydrodynamic_evoution[ivar][it+1][ix][iy][iz], hydrodynamic_evoution[ivar][it][ix+1][iy][iz], hydrodynamic_evoution[ivar][it][ix][iy+1][iz], hydrodynamic_evoution[ivar][it][ix][iy][iz+1],
    hydrodynamic_evoution[ivar][it+1][ix+1][iy][iz], hydrodynamic_evoution[ivar][it+1][ix][iy+1][iz], hydrodynamic_evoution[ivar][it+1][ix][iy][iz+1],
    hydrodynamic_evoution[ivar][it][ix+1][iy+1][iz], hydrodynamic_evoution[ivar][it][ix+1][iy][iz+1], hydrodynamic_evoution[ivar][it][ix][iy+1][iz+1],
    hydrodynamic_evoution[ivar][it+1][ix+1][iy+1][iz], hydrodynamic_evoution[ivar][it+1][ix+1][iy][iz+1], hydrodynamic_evoution[ivar][it][ix+1][iy+1][iz+1], hydrodynamic_evoution[ivar][it+1][ix][iy+1][iz+1], hydrodynamic_evoution[ivar][it+1][ix+1][iy+1][iz+1]);
    return result;
}

float interpolateVariable3D(float *****hydrodynamic_evoution, int ivar, int it, int ix, int iy, float tau_frac, float x_frac, float y_frac)
{
  float result = linearInterp3D(tau_frac, x_frac, y_frac,
    hydrodynamic_evoution[ivar][it][ix][iy][0], hydrodynamic_evoution[ivar][it+1][ix][iy][0], hydrodynamic_evoution[ivar][it][ix+1][iy][0], hydrodynamic_evoution[ivar][it][ix][iy+1][0],
    hydrodynamic_evoution[ivar][it+1][ix+1][iy][0], hydrodynamic_evoution[ivar][it+1][ix][iy+1][0], hydrodynamic_evoution[ivar][it][ix+1][iy+1][0], hydrodynamic_evoution[ivar][it+1][ix+1][iy+1][0]);
    return result;
}

void callFOFinder3p1D(int dim, int nx, int ny, int nz, int n, float t0, float dt, float t, float dx, float dy, float dz, float *lattice_spacing, float freezeoutEnergyDensity,
  float ****hyperCube4D, float ***hyperCube3D, float ****energy_density_evoution, float *****hydrodynamic_evoution,
  std::ofstream& freezeoutSurfaceFile, std::vector<FO_Element>& fo_surf)
{

  //besides writing centroid and normal to file, write all the hydro variables
  //#pragma omp parallel for collapse(3)
  for (int ix = 0; ix < nx-1; ix++)
  {
    for (int iy = 0; iy < ny-1; iy++)
    {
      for (int iz = 0; iz < nz-1; iz++)
      {
        Cornelius cor;
        cor.init(dim, freezeoutEnergyDensity, lattice_spacing);

        //write the values of energy density to all corners of the hyperCube
        writeEnergyDensityToHypercube4D(hyperCube4D, energy_density_evoution, 0, ix, iy, iz);
        //use cornelius to find the centroid and normal vector of each hyperCube
        cor.find_surface_4d(hyperCube4D);
        //write centroid and normal of each surface element to file
        for (int i = 0; i < cor.get_Nelements(); i++)
        {
          //declare a new fo cell to hold info, later push back to vector
          FO_Element fo_cell;

          float temp = 0.0; //temporary variable
          //first write the position of the centroid of surface element
          float cell_tau = t0 + ((float)n) * dt; //check if this is the correct time!
          float cell_x = (float)ix * dx  - (((float)(nx-1)) / 2.0 * dx);
          float cell_y = (float)iy * dy  - (((float)(ny-1)) / 2.0 * dy);
          float cell_z = (float)iz * dz  - (((float)(nz-1)) / 2.0 * dz);

          float tau_frac = cor.get_centroid_elem(i,0) / lattice_spacing[0];
          float x_frac = cor.get_centroid_elem(i,1) / lattice_spacing[1];
          float y_frac = cor.get_centroid_elem(i,2) / lattice_spacing[2];
          float z_frac = cor.get_centroid_elem(i,3) / lattice_spacing[3];

          float tau = cor.get_centroid_elem(i,0) + cell_tau;
          float x = cor.get_centroid_elem(i,1) + cell_x;
          float y = cor.get_centroid_elem(i,2) + cell_y;
          float eta = cor.get_centroid_elem(i,3) + cell_z;

          float ds0 = t * cor.get_normal_elem(i,0);
          float ds1 = t * cor.get_normal_elem(i,1);
          float ds2 = t * cor.get_normal_elem(i,2);
          float ds3 = t * cor.get_normal_elem(i,3);

          float ux = interpolateVariable4D(hydrodynamic_evoution, 0, 0, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
          float uy = interpolateVariable4D(hydrodynamic_evoution, 1, 0, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
          float un = interpolateVariable4D(hydrodynamic_evoution, 2, 0, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);

          float eps = interpolateVariable4D(hydrodynamic_evoution, 3, 0, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
          float T = effectiveTemperature(eps);
          float P = equilibriumPressure(eps);

          float pixx = interpolateVariable4D(hydrodynamic_evoution, 4, 0, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
          float pixy = interpolateVariable4D(hydrodynamic_evoution, 5, 0, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
          float pixn = interpolateVariable4D(hydrodynamic_evoution, 6, 0, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
          float piyy = interpolateVariable4D(hydrodynamic_evoution, 7, 0, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
          float piyn = interpolateVariable4D(hydrodynamic_evoution, 8, 0, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);

          float Pi = interpolateVariable4D(hydrodynamic_evoution, 9, 0, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);

          #pragma omp critical
          freezeoutSurfaceFile << tau  << " " <<  x   << " " <<  y   << " " << eta  << " "
                               << ds0  << " " << ds1  << " " << ds2  << " " << ds3  << " "
                               <<  ux  << " " <<  uy  << " " <<  un  << " "
                               << eps  << " " <<  T   << " " <<  P   << " "
                               << pixx << " " << pixy << " " << pixn << " " << piyy << " " << piyn << " "
                               <<  Pi  << endl;

          //add the fo cell to fo surface
          //#pragma omp critical
          //if (SAVE_FO_SURF_VECTOR) fo_surf.push_back(fo_cell);

        } //for (int i = 0; i < cor.get_Nelements(); i++)
      } // for (int iz = 0; iz < dimZ; iz++)
    } // for (int iy = 0; iy < ny-1; iy++)
  } // for (int ix = 0; ix < nx-1; ix++)
}

void callFOFinder2p1D(int dim, int nx, int ny, int nz, int n, float t0, float dt, float t, float dx, float dy, float dz, float *lattice_spacing, float freezeoutEnergyDensity,
  float ****hyperCube4D, float ***hyperCube3D, float ****energy_density_evoution, float *****hydrodynamic_evoution,
  std::ofstream& freezeoutSurfaceFile, std::vector<FO_Element>& fo_surf)
  {
    //besides writing centroid and normal to file, write all the hydro variables
    //#pragma omp parallel for collapse(2)
    for (int ix = 0; ix < nx-1; ix++)
    {
      for (int iy = 0; iy < ny-1; iy++)
      {
        Cornelius cor;
        cor.init(dim, freezeoutEnergyDensity, lattice_spacing);
        //write the values of energy density to all corners of the hyperCube
        writeEnergyDensityToHypercube3D(hyperCube3D, energy_density_evoution, 0, ix, iy);
        //use cornelius to find the centroid and normal vector of each hyperCube
        cor.find_surface_3d(hyperCube3D);
        //write centroid and normal of each surface element to file
        for (int i = 0; i < cor.get_Nelements(); i++)
        {
          //declare a new fo cell to hold info, later push back to vector
          FO_Element fo_cell;

          //first write the position of the centroid of surface element
          float cell_tau = t0 + ((float)n) * dt; //check if this is the correct time!
          float cell_x = (float)ix * dx  - (((float)(nx-1)) / 2.0 * dx);
          float cell_y = (float)iy * dy  - (((float)(ny-1)) / 2.0 * dy);
          float cell_z = 0.0;

          float tau_frac = cor.get_centroid_elem(i,0) / lattice_spacing[0];
          float x_frac = cor.get_centroid_elem(i,1) / lattice_spacing[1];
          float y_frac = cor.get_centroid_elem(i,2) / lattice_spacing[2];

          //cell position
          float tau = cor.get_centroid_elem(i,0) + cell_tau;
          float x = cor.get_centroid_elem(i,1) + cell_x;
          float y = cor.get_centroid_elem(i,2) + cell_y;
          float eta = 0.0;

          //covariant surface normal vector
          float ds0 = t * cor.get_normal_elem(i,0);
          float ds1 = t * cor.get_normal_elem(i,1);
          float ds2 = t * cor.get_normal_elem(i,2);
          float ds3 = 0.0;

          //contravariant flow velocity
          float ux = interpolateVariable3D(hydrodynamic_evoution, 0, 0, ix, iy, tau_frac, x_frac, y_frac);
          float uy = interpolateVariable3D(hydrodynamic_evoution, 1, 0, ix, iy, tau_frac, x_frac, y_frac);
          float un = interpolateVariable3D(hydrodynamic_evoution, 2, 0, ix, iy, tau_frac, x_frac, y_frac);

          //energy density, Temperature, Pressure
          float eps = interpolateVariable3D(hydrodynamic_evoution, 3, 0, ix, iy, tau_frac, x_frac, y_frac);
          float T = effectiveTemperature(eps);
          float P = equilibriumPressure(eps);

          //contravariant components of shear stress
          float pixx = interpolateVariable3D(hydrodynamic_evoution, 4, 0, ix, iy, tau_frac, x_frac, y_frac);
          float pixy = interpolateVariable3D(hydrodynamic_evoution, 5, 0, ix, iy, tau_frac, x_frac, y_frac);
          float pixn = interpolateVariable3D(hydrodynamic_evoution, 6, 0, ix, iy, tau_frac, x_frac, y_frac);
          float piyy = interpolateVariable3D(hydrodynamic_evoution, 7, 0, ix, iy, tau_frac, x_frac, y_frac);
          float piyn = interpolateVariable3D(hydrodynamic_evoution, 8, 0, ix, iy, tau_frac, x_frac, y_frac);
          //bulk pressure
          float Pi = interpolateVariable3D(hydrodynamic_evoution, 9, 0, ix, iy, tau_frac, x_frac, y_frac);

          #pragma omp critical
          {
            freezeoutSurfaceFile << tau  << " " <<  x   << " " <<  y   << " " << eta  << " "
                                 << ds0  << " " << ds1  << " " << ds2  << " " << ds3  << " "
                                 <<  ux  << " " <<  uy  << " " <<  un  << " "
                                 << eps  << " " <<  T   << " " <<  P   << " "
                                 << pixx << " " << pixy << " " << pixn << " " << piyy << " " << piyn << " "
                                 <<  Pi  << endl;
          }
          //add the fo cell to fo surface
          //#pragma omp critical
          //if (SAVE_FO_SURF_VECTOR) fo_surf.push_back(fo_cell);

        } //for (int i = 0; i < cor.get_Nelements(); i++)
      } // for (int iy = 0; iy < ny-1; iy++)
    } // for (int ix = 0; ix < nx-1; ix++)
  }

//returns the number of cells with T > T_c
int checkForCellsAboveTc(int nx, int ny, int nz, float freezeoutEnergyDensity, PRECISION *e)
{
  int accumulator = 0;

  #pragma omp parallel for collapse(3)
  for (int ix = 2; ix < nx+2; ix++)
  {
    for (int iy = 2; iy < ny+2; iy++)
    {
      for (int iz = 2; iz < nz+2; iz++)
      {
        int s = columnMajorLinearIndex(ix, iy, iz, nx+4, ny+4, nz+4);
        if (e[s] > freezeoutEnergyDensity) accumulator += 1;
      }
    }
  }
  return accumulator;
}
