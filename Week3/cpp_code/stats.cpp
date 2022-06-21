#include "stats.h"


// ------- COSTRUTTORE E DISTRUTTORE -------
stats::stats(){
}

stats::~stats(){
}


// ------- METODI -------
// media a blocchi
void stats::block_mean(vector<double>& v){
	
	double aux = 0, aux2 = 0;
	
	clear();
   	for(int i=0; i<int(v.size()); i++){

      	// somme progressive
      	aux = 0;
      	aux2 = 0;
      
      	for(int k=0; k<=i; k++){
      		aux += v[k];
			aux2 += pow(v[k], 2);
		}
		
		aux /= i+1;
		aux2 /= i+1;

      	// calcolo media con incertezza
      	m_mean_prog.push_back( aux );
      	m_err_prog.push_back( sqrt((aux2 - aux*aux)/(i+1)) );
      		
	}
      		
}

// salva i vector protetti
void stats::save(string filename){

	ofstream fout(filename);
	
	fout << "Mean \t Error" << endl;
	for(int i=0; i<int(m_mean_prog.size()); i++)
		fout << m_mean_prog[i] << "\t" << m_err_prog[i] << endl;
	
	fout.close();
}

// pulizia degli m_prog
void stats::clear(){

	m_mean_prog.clear();
	m_err_prog.clear();
	
}
