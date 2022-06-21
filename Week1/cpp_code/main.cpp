#include <iostream>
#include <string>
#include <vector>

#include "random.h"
#include "stats.h"
#include "utils.h"

using namespace std;
 

int main (int argc, char *argv[]){

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

   
   
   // ********* Exercise 01.1 *********
   int M1 = 1000000, N1 = 100, n_int = 100;
   double mean_teo = 0.5;
   
   stats gen;
   
   vector<double> extr;
   vector<double> mean, var;
   vector<double> chi_sq;
   
   
   for(int i=0; i<N1; i++){
   
   	// azzero le estrazioni
   	extr.clear();
   	
   	// genero M numeri casuali in blocchi da M/N (N blocchi)
        extr = rnd.Rannyu( int(M1/N1) );
        
        // salvo media e varianza del singolo blocco
        mean.push_back( get_mean(extr) );
        var.push_back( get_var(extr, mean_teo) );
        	
        // calcolo il chi quadro del singolo blocco
        chi_sq.push_back( get_chi_sq( extr, M1/N1/double(n_int), 0., 1., n_int ) );
        
   }
   
   // calcolo (e salvo su file) la media a blocchi per "mean"...
   gen.block_mean( mean );
   gen.save("../data/data_01.1/mean_data.txt");
   
   // ... e per "dev"
   gen.block_mean( var );
   gen.save("../data/data_01.1/var_data.txt");
   
   // salvo su file il chi quadro
   save(chi_sq, "../data/data_01.1/chi_data.txt", "Chi squared");
   
   
   
   // ********* Exercise 01.2 *********
   double lambda = 1;
   double gamma = 1, mu = 0;
   
   int M2 = 10000;
   int Ns[] = {1, 2, 10, 100};
   
   vector<double> v_unif, v_expo, v_cauc;
   
   ofstream fout;
   
   
   for(int i=0; i<4; i++){
   	
   	// apertura file
   	if(Ns[i] == 1)
		fout.open("../data/data_01.2/Sampling_1.txt");
	else if(Ns[i] == 2)
		fout.open("../data/data_01.2/Sampling_2.txt");
	else if(Ns[i] == 10)
		fout.open("../data/data_01.2/Sampling_10.txt");
	else
		fout.open("../data/data_01.2/Sampling_100.txt");
		
	// azzeramento
	v_unif.clear();
	v_expo.clear();
	v_cauc.clear();
		
	// intestazione
	fout << "Uniform" << "\t" << "Exponential" << "\t" << "Cauchy" << endl;  
   	
   	for(int j=0; j<M2; j++){
   	
   		// genero Ns[i] numeri casuali dalle tre distribuzioni
   		v_unif = rnd.Rannyu(1., 6., Ns[i]);
   		v_expo = rnd.Expo(lambda, Ns[i]);
   		v_cauc = rnd.Cauchy(gamma, mu, Ns[i]);
   		
   		// salvo il valor medio dei numeri generati delle tre distribuzioni (n volte)
   		fout << get_mean( v_unif ) << "\t" << get_mean( v_expo ) << "\t" << get_mean( v_cauc ) << endl;
   		
   	}
   	
   	// chiusura file
   	fout.close();
   	
   }
   
   
   
   // ********* Exercise 01.3 *********
   int M3 = 10000, N3 = 100, N_lines = 50;
   
   double throw_ = 0, x = 0, y = 0, cos_theta = 0;	
   double d = 1., D = d/2., const_ = 2.*D/d;		// in questo modo le linee sono i numeri naturali [0, N_lines)
   
   // riutilizzo "extr" e "mean" per il calcolo di pi greco
   mean.clear();
   
   for(int i=0; i<N3; i++){
   	
   	// azzeramento delle hit
   	extr.clear();
   	for(int k=0; k<N3; k++)
   		extr.push_back(0);
   	
   	// calcolo le hit (verificando di aver colpito almeno una volta)
   	do{
   		for(int j=0; j<int(M3/N3); j++){
   			
   			// effettuo un lancio (throw_ == centro dell'asticella)
   			throw_ = rnd.Rannyu(0, N_lines);
   			
   			// estraggo un angolo tra [-pi/2, pi/2] per ruotare l'asticella
   			do{
   				x = rnd.Rannyu();
   				y = rnd.Rannyu(-1, 1);
   			}while( x*x + y*y > 1 );
   			cos_theta = x/sqrt(x*x + y*y);
   			
   			// controllo che l'asticella sia a cavallo di una riga
   			if( int(throw_-D*cos_theta/2.) < int(throw_+D*cos_theta/2.) )
   				extr[j] ++;
   		}
   	}while(get_sum(extr) == 0);
   	
   	// salvo media del singolo blocco
       mean.push_back( const_ * 1./get_mean(extr) );
       
   }
       
   // calcolo (e salvo su file) la media a blocchi per "mean"
   gen.block_mean( mean );
   gen.save("../data/data_01.3/pi_data.txt");
   
   
   
   rnd.SaveSeed();
   return 0;
}
