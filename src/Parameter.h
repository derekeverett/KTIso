#pragma once
struct parameters
{
  int output_format;
  int ic_energy;
  int ic_flow;
  int nx;
  int ny;
  int nphip;
  int nvz;
  //int nt;
  float dx;
  float dy;
  float dt;
  float t0;
  float tf;
  int eos_type;
  float e_sw;
  float eta_over_s;
  int collisions;
  int sources;
  int adapt_time;
  float angular_acc_factor;
  float fs_acc_factor;
  float coll_acc_factor;
  //these are computed based on the chosen parameters above; they are constrained
  int ntot;
  //float TAU;

  //useful constants
  float hbarc;
};
