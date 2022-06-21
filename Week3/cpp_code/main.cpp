#include <iostream>
#include <string>
#include <vector>
#include <string.h>

#include "random.h"
#include "utils.h" 
#include "stats.h"

using namespace std;


int main (int argc, char *argv[]){

   // ********* Generatore numeri pseudo-casuali *********
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
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



   // ********* Exercise 03.1 *********
   int M = 10000, N = 100;
   double T = 1.;

   int nsteps[] = {1, 1000};

   vector<double> aux;
   vector<double> mean_C, mean_P;

   stats gen;

   // ciclo sul numero di step
   for(int i=0; i<2; i++)
   {

      mean_C.clear();
      mean_P.clear();

      // calcolo della media per ciascun blocco
      for(int j=0; j<N; j++)
      {
         aux = black_scholes(T, nsteps[i], int(M/N), rnd);
         mean_C.push_back(aux[0]);
         mean_P.push_back(aux[1]); 
      }

      // media a blocchi con errore
      gen.block_mean(mean_C);
      gen.save("../data/mean_C_"+to_string(nsteps[i])+".dat");

      gen.block_mean(mean_P);
      gen.save("../data/mean_P_"+to_string(nsteps[i])+".dat");

   }
   
   rnd.SaveSeed();
   return 0;
}
