#include "random.h"

// ------- COSTRUTTORE E DISTRUTTORE -------
Random :: Random(){
}

Random :: ~Random(){
}


// ------- METODI -------
// numero casuale tra 0 ed 1
double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

// numero casuale tra min e max
double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

// estrazione da gaussiana
double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

// estrazione da lorentziana
double Random :: Cauchy(double gamma, double mu) {
   double x = Rannyu();
   return mu + gamma*tan( M_PI*( x-(1./2.) ) );
}

// estrazione da esponenziale
double Random :: Expo(double lambda) {
   double x = Rannyu();
   return (-1./lambda)*log(1-x);
}




// come sopra, ma con n estrazioni
vector<double> Random :: Rannyu(int n){
  vector <double> v;
  for(int i=0; i<n; i++)
  	v.push_back( Rannyu() );
  return v;
}

vector<double> Random :: Rannyu(double min, double max, int n){
  vector <double> v;
  for(int i=0; i<n; i++)
  	v.push_back( Rannyu(min, max) );
  return v;
}

vector<double> Random :: Gauss(double mean, double sigma, int n){
  vector <double> v;
  for(int i=0; i<n; i++)
  	v.push_back( Gauss(mean, sigma) );
  return v;
}

vector<double> Random :: Cauchy(double gamma, double mu, int n){
  vector <double> v;
  for(int i=0; i<n; i++)
  	v.push_back( Cauchy(gamma, mu) );
  return v;
}

vector<double> Random :: Expo(double lambda, int n){
  vector <double> v;
  for(int i=0; i<n; i++)
  	v.push_back( Expo(lambda) );
  return v;
}



// estrazione di un versore 3D (da pdf sferica uniforme)
vector<double> Random :: tridim_versor(){

	double phi = 2*M_PI * Rannyu();
	double theta = acos( 1 - 2*Rannyu() );
	
	vector<double> v;
	v.push_back( cos(phi)*sin(theta) );
	v.push_back( sin(phi)*sin(theta) );
	v.push_back( cos(theta) );
	
	return v;

}

// estrazione di un versore 2D (da pdf sferica uniforme)
vector<double> Random :: bidim_versor(){

	double phi = 2*M_PI * Rannyu();
	
	vector<double> v;
	v.push_back( cos(phi) );
	v.push_back( sin(phi) );
	
	return v;

}
  	



// ??
void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

// salva l'ultimo seme
void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}
