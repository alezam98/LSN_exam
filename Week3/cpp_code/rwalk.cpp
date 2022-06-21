#include "rwalk.h"

	//COSTRUTTORI E DISTRUTTORE

// costruttore senza argomenti
rwalk::rwalk() {
	double aux[] = {0, 0, 0};
	m_a = 1.;
	m_rw.push_back(aux);
	m_rN2.push_back(0);
}

// costruttore, lunghezza del passo
rwalk::rwalk(double a){

	if( a <= 0 ){
		cout << "Errore: il passo deve essere positivo!" <<endl;
		exit(4);
	}

	else{
		double aux[] = {0, 0, 0};
		m_a = a;
		m_rw.push_back(aux);
		m_rN2.push_back(0);
	}

}

// distruttore
rwalk::~rwalk() {
}



	// METODI

// singolo passo, random walk discreto (continous = 0) e continuo (continous = 1)
void rwalk::step(Random& rnd, bool continous){

	double* v = new double(3);

	if(continous){
	
		v = rnd.versor();
		
		for(int i=0; i<3; i++){
			v[i] *= m_a;
			v[i] += m_rw[ m_rw.size()-1 ][i];
		}
		
	}
	else{
	
		int dir = int( rnd.Rannyu(0., 3.) );
		int sign = 2 * int( rnd.Rannyu(0., 2.) ) - 1;
		
		for(int i=0; i<3; i++){
			v[i] = 0;
			v[i] += m_rw[ m_rw.size()-1 ][i];
		}
		v[dir] += m_a * sign;
		
	}
	
	m_rw.push_back(v);
	m_rN2.push_back( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
	
	
}

// n passi, random walk discreto (continous = 0) e continuo (continous = 1)
void rwalk::step(Random& rnd, bool continous, int n_step){

	for(int i=0; i<n_step; i++)
		step(rnd, continous);
	
}

// varianza media su n_rw random walk di n_step passi
vector<double> rwalk::sigma2(Random& rnd, int n_step, int n_rw, bool continous){

	vector<double> mean;
	for(int i=0; i<n_step+1; i++)
		mean.push_back(0);
	
	for(int i=0; i<n_rw; i++){
		
		clear();
		step(rnd, continous, n_step);
		mean = mean + m_rN2;
	
	}
	
	return mean/n_rw;
	
}

// azzeramento
void rwalk::clear(){

	m_rw.clear();
	m_rN2.clear();
	
	double aux[] = {0, 0, 0};
	m_rw.push_back(aux);
	m_rN2.push_back(0);
	
}


// stampa a video del rw
void rwalk::print(){

	cout << "Coordinate del random walk in funzione del passo:" << endl;
	for(int i=0; i<int(m_rw.size()); i++)
		cout << i << ": " << m_rw[i][0] << "\t" << m_rw[i][1] << "\t" << m_rw[i][2] << endl;
	cout << endl;
		
}