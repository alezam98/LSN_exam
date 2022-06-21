#include <iostream>
#include <string>
#include <vector>

#include "random.h"
#include "utils.h" 
#include "stats.h"
#include "rwalk.h"

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



   // ********* Exercise 02.1 *********
   int M1 = 100000, N1 = 100, n_pts = 1000;

   vector<double> I_media, I_sample, I_ar;
   vector<double> media_mean, sample_mean, ar_mean;
   
   stats gen_media;
   stats gen_sample;
   stats gen_ar;

   for(int i=0; i<N1; i++){
   	
      // calcolo M integrali in blocchi da M/N (N blocchi)
      for(int j=0; j<int(M1/N1); j++){

         I_media.push_back( integrale_media(n_pts, rnd) );
         I_sample.push_back( integrale_sample(n_pts, rnd) );
         I_ar.push_back( integrale_ar(n_pts, rnd) );

      }

      media_mean.push_back( get_mean(I_media) );
      sample_mean.push_back( get_mean(I_sample) );
      ar_mean.push_back( get_mean(I_ar) );

   }
   
   // media a blocchi
   gen_media.block_mean( media_mean );
   gen_sample.block_mean( sample_mean );
   gen_ar.block_mean( ar_mean );
   
   
   // salvataggio dati
   gen_media.save("../data/data_02.1/media_data.dat");
   gen_sample.save("../data/data_02.1/sample_data.dat");
   gen_ar.save("../data/data_02.1/ar_data.dat");
   
   
   
   // ********* Exercise 02.2 *********
   int M2 = 10000, N2 = 100, n_step = 100, it = 0;
   bool continous = 0;
   
   vector<double> aux;
   vector<double> mean, mean2, error;
   
   rwalk random_walk(3.);

   
   // calcolo di sigma per rw discreto (continous = 0) e continuo (continous = 1)
   while(it <= 1){

      // inizializzazione
      mean.clear();
      mean2.clear();
      for(int i=0; i<n_step+1; i++){
		   mean.push_back(0);
		   mean2.push_back(0);
      }

      // N2 blocchi, ciascuno da M2/N2 random walk da n_step=100 passi
      for(int i=0; i<N2; i++){
         
         aux = random_walk.sigma2(rnd, n_step, int(M2/N2), continous);
         mean = mean + aux;
         mean2 = mean2 + (aux^2);
      
      }
      
      mean = mean/N2;
      mean2 = mean2/N2;
      
      // errore di sigma2
      error = ( (mean2 - (mean^2))/N2 )^0.5;
      
      
      // da sigma2 a sigma
      mean.erase( mean.begin() );
      error.erase( error.begin() );
      
      mean = mean^0.5;
      error = (error/mean)/2.;

      mean.insert(mean.begin(), 0);
      error.insert(error.begin(), 0);
      
      // salvataggio dati
      if(continous == 0)
         save(mean, error, "../data/data_02.2/discrete.dat", "Mean", "Error");
      else
         save(mean, error, "../data/data_02.2/continous.dat", "Mean", "Error");

      // da discreto a continuo
      it ++;
      continous = it;

   }
   

   rnd.SaveSeed();
   return 0;
}