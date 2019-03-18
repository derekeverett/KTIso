#pragma once
struct parameters
{
  int OUTPUTFORMAT;
  int IC_ENERGY;
  int DIM_X;
  int DIM_Y;
  int DIM_PHIP;
  int DIM_VZ;
  int DIM_T;
  float DX;
  float DY;
  float DT;
  float T0;
  int EOS_TYPE;
  float E_FREEZE;
  float ALPHA;
  //these are computed based on the chosen parameters above; they are constrained
  int DIM;
  //float TAU;
};
