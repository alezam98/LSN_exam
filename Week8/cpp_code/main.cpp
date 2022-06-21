#include "main.h"

using namespace std;

int main (int argc, char *argv[])
{
  // loading data
  Input();
  

  // SA annealing
  if(verbose) cout << "Start of the simulation" << endl;
  best_params = annealer->annealing(x, Ts, rnd, gen, 1);
  annealer->save_HTs("./"+dir+"/HTs.dat");
  annealer->save_params("./"+dir+"/params.dat");


  // MC sampling with best parameters
  pdf->set_params(best_params);
  annealer->set_pdf(pdf);
  x = mu;

  for(int blk=0; blk<nblk; blk++) mean.push_back(annealer->metropolis(x, rnd, 0));
  gen.block_mean(mean);
  gen.save("./"+dir+"/best_HT.dat");
  annealer->save_xs("./"+dir+"/xs.dat");


  // print
  cout << endl;
  cout << "Best parameters:" << endl;
  cout << "mu = " << best_params[0] << ",\tsigma = " << best_params[1] << endl << endl;

  cout << "Best <HT> with statistical error:" << endl;
  cout << "<HT> = " << gen.get_last_mean() << " +- " << gen.get_last_err() << endl;
  cout << endl;


  rnd.SaveSeed();
  return 0;
}



void Input(void)
{
  ifstream ReadInput, Primes, Seed;

  //Read input informations
  ReadInput.open("input.in");

  ReadInput >> verbose;

  if(verbose)
  {
    cout << "Monodimensional quantum particle" << endl;
    cout << "Simulated annealing simulation, performing Metropolis moves with uniform translations" << endl << endl;

    cout << "Hamiltonian: H = -1/2 * d^2/dx^2 + V(x), where V(x) = (x)^4 - (5/2)*x^2" << endl;
    cout << "Test wave function: psiT = 1/N * (exp(-(x-mu)^2/(2*sigma^2)) + exp(-(x+mu)^2/(2*sigma^2)))" << endl;
    cout << "The program uses h/(2*pi)=1 and m=1 " << endl << endl;
  }

  ReadInput >> restart;

  //Read seed for random numbers
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  if(restart) Seed.open("seed.out");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  ReadInput >> nblk;

  ReadInput >> mc_steps;

  ReadInput >> x_width;

  ReadInput >> sa_steps;

  ReadInput >> mu;

  ReadInput >> sigma;

  ReadInput >> mu_width;

  ReadInput >> sigma_width;

  ReadInput >> T_steps;

  ReadInput >> T;

  ReadInput >> T_width;

  ReadInput >> dir;

  // posizione di partenza
  x = mu;

  // inizializzazione vector
  best_params.push_back(mu);
  best_params.push_back(sigma);
  Ts.push_back(T);
  for(int step=1; step<T_steps; step++) Ts.push_back(Ts[step-1] - T_width);

  // inizializzazione classi
  pdf->set_params(best_params);
  annealer->set_pdf(pdf);
  annealer->set_mc_steps(mc_steps);
  annealer->set_x_width(x_width);
  annealer->set_sa_steps(sa_steps);
  annealer->set_nblk(nblk);
  annealer->set_mu_width(mu_width);
  annealer->set_sigma_width(sigma_width);


  if(verbose)
  {
    cout << "-- Simulation parameters" << endl;
    cout << "Blocks:                 " << nblk << endl;
    cout << "MC steps per block:     " << mc_steps << endl;
    cout << "SA steps with T fixed:  " << sa_steps << endl;
    cout << "Number of temperatures: " << T_steps << endl << endl;

    cout << "-- Starting configuration" << endl;
    cout << "Starting x:     " << x << endl;
    cout << "Starting mu:    " << mu << endl;
    cout << "Starting sigma: " << sigma << endl;
    cout << "Starting T:     " << T << endl << endl;

    cout << "-- Move parameters" << endl;
    cout << "dx:     " << x_width << endl;
    cout << "dmu:    " << mu_width << endl;
    cout << "dsigma: " << sigma_width << endl;
    cout << "dT:     " << T_width << endl << endl;
  }
  
  ReadInput.close();

  return;
}
