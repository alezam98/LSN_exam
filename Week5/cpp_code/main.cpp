#include <iostream>
#include <string>
#include <vector>

#include "random.h"
#include "utils.h" 
#include "stats.h"
#include "functions.h"
#include "algorithm.h"

using namespace std;


int main (int argc, char *argv[]){

   // ********* Generatore numeri pseudo-casuali *********
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2;
      Primes >> p1 >> p2;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;



   // ********* Exercise 05.1 *********
   int M = 100000, N = 100, steps = int(M/N);
   double prec = 0.1;
   bool save;

   double pars[] = {0., 0.};
   string Ts[] = {"uniform", "gaussian"};

   vector<double> mean_r;
   vector<double> config;

   psi_100* pdf_gs = new psi_100();       // pdf, ground state and first excited state
   psi_210* pdf_1e = new psi_210();

   algorithm* metropolis = new algorithm(pdf_gs, 1.);

   stats gen;

   
   // GROUND STATE
   pars[0] = 1.25; 
   pars[1] = 0.75;

   for(int iT=0; iT<2; iT++)
   {
      cout << "- Ground state (1s), " << Ts[iT] << " sampling:" << endl;

      // settings
      config.assign(3, 20);
      mean_r.clear();
      metropolis->set_T(Ts[iT]);
      metropolis->set_par(pars[iT]);

      // equilibration
      cout << "Starting equilibration of the system..." << endl;
      save = 1;
      config = metropolis->equilibration(config, rnd, prec, save);
      metropolis->save("../data/ground_state/eq_configs_"+Ts[iT]+".dat");

      // metropolis
      cout << "Metropolis algorithm..." << endl;
      save = 0;
      for(int i=0; i<N; i++)
         mean_r.push_back( metropolis->metropolis(config, rnd, steps, save) );

      // saving data
      gen.block_mean(mean_r);
      gen.save("../data/ground_state/r_"+Ts[iT]+".dat");

      // 3D space exploration
      cout << "3D space exploration..." << endl;
      save = 1;
      metropolis->metropolis(config, rnd, 10*steps, save);
      metropolis->save("../data/ground_state/configs_"+Ts[iT]+".dat");

      cout << "... done." << endl << endl;
   }


   // EXCITED STATE
   pars[0] = 3.; 
   pars[1] = 1.75;
   metropolis->set_pdf(pdf_1e);

   for(int iT=0; iT<2; iT++)
   {
      cout << "- Excited state (2p), " << Ts[iT] << " sampling:" << endl;

      // settings
      config.assign(3, 30);
      mean_r.clear();
      metropolis->set_T(Ts[iT]);
      metropolis->set_par(pars[iT]);

      // equilibration
      cout << "Starting equilibration of the system..." << endl;
      save = 1;
      config = metropolis->equilibration(config, rnd, prec, save);
      metropolis->save("../data/excited_state/eq_configs_"+Ts[iT]+".dat");

      // metropolis
      cout << "Metropolis algorithm..." << endl;
      save = 0;
      for(int i=0; i<N; i++)
         mean_r.push_back( metropolis->metropolis(config, rnd, steps, save) );

      // saving data
      gen.block_mean(mean_r);
      gen.save("../data/excited_state/r_"+Ts[iT]+".dat");

      // 3D space exploration
      cout << "3D space exploration..." << endl;
      save = 1;
      metropolis->metropolis(config, rnd, 10*steps, save);
      metropolis->save("../data/excited_state/configs_"+Ts[iT]+".dat");
            
      cout << "... done." << endl << endl;
   }

   
   rnd.SaveSeed();
   return 0;
}
