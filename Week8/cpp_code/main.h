#ifndef __var_MC__
#define __var_MC__

#include <iostream>
#include <string>
#include <vector>

#include "random.h"
#include "utils.h" 
#include "stats.h"
#include "functions.h"
#include "algorithm.h"

using namespace std;


//Random numbers
#include "random.h"
int seed[4];
Random rnd;

// metropolis
double x, x_width;

// annealing
double T, T_width;
double mu, mu_width;
double sigma, sigma_width;
vector<double> best_params, Ts;

// best_parameters MC sampling
vector<double> mean;

// simulation
int nblk, mc_steps, sa_steps, T_steps, verbose, restart;
string dir;

// classi
stats gen;
psiT2* pdf = new psiT2();
SA_algorithm* annealer = new SA_algorithm(pdf);

//functions
void Input(void);


#endif