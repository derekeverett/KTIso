#include "KTiso.cpp"

int main(void)
{
  //Declare an instance of ITA class
  ITA isotropize;

  //try to initialize the energy density from a vector
  /*
  int npoints = 101 * 101;
  size_t size = npoints * sizeof(double);
  std::vector<double> init_e(size); //initial energy density vector
  std::vector<double> final_e(size); //final energy density vector
  std::vector<double> final_p(size); //final pressure vector
  std::vector<double> final_ut(size); // etc...
  std::vector<double> final_ux(size);
  std::vector<double> final_uy(size);
  std::vector<double> final_un(size);
  std::vector<double> final_pitt(size);
  std::vector<double> final_pitx(size);
  std::vector<double> final_pity(size);
  std::vector<double> final_pitn(size);
  std::vector<double> final_pixx(size);
  std::vector<double> final_pixy(size);
  std::vector<double> final_pixn(size);
  std::vector<double> final_piyy(size);
  std::vector<double> final_piyn(size);
  std::vector<double> final_pinn(size);
  std::vector<double> final_Pi(size);
  */
  
  //run the evolution
  isotropize.run_ita();

  /*
  //grab the final hydro vectors to pass to another module
  fsmilne.output_to_vectors(final_e,
                            final_p,
                            final_ut,
                            final_ux,
                            final_uy,
                            final_un,
                            final_pitt,
                            final_pitx,
                            final_pity,
                            final_pitn,
                            final_pixx,
                            final_pixy,
                            final_pixn,
                            final_piyy,
                            final_piyn,
                            final_pinn,
                            final_Pi);
  */
  printf("run_ita() ran sucessfully \n");

  /*
  //check that the vectors were upated
  for (int i = 0; i < npoints; i++)
  {
    printf("e [ %d ] = %f \n", i, final_e[i]);
  }
  */
}
