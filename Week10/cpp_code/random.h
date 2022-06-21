#ifndef __Random__
#define __Random__


#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;


class Random {


  public:
  
  	Random();
  	~Random();
  	
  	void SetRandom(int* , int, int);
  	void SaveSeed();
  	
  	double Rannyu(void);
  	double Rannyu(double min, double max);
  	double Gauss(double mean, double sigma);
  	double Expo(double lambda);
  	double Cauchy(double gamma, double mu);
  	
  	vector<double> Rannyu(int n);
  	vector<double> Rannyu(double min, double max, int n);
  	vector<double> Gauss(double mean, double sigma, int n);
  	vector<double> Expo(double lambda, int n);
  	vector<double> Cauchy(double gamma, double mu, int n);
  	
  	vector<double> tridim_versor();
	vector<double> bidim_versor();

  private:
  
	int m1, m2, m3, m4, l1, l2, l3, l4, n1, n2, n3, n4;
  
};


#endif
