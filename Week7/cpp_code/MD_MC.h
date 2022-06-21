/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid__
#define __fluid__

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, iv, ik, it, ie, iw, ip;
int r_bins;
double vtail, ptail, bin_size, nbins, sd;
double walker[m_props];
vector<double> Norm;

// averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, stima_kin, stima_etot, stima_temp, stima_press;
double err_pot, err_kin, err_etot, err_temp, err_press, err_gdir;
vector<double> blk_g, glob_g, glob_g2, stima_g;
double err_g;

//configuration
const int m_part=108;
double x[m_part],    y[m_part],    z[m_part];
double xold[m_part], yold[m_part], zold[m_part];
double vx[m_part],  vy[m_part],   vz[m_part];

// equilibration
int MD_eqsteps, MC_eqsteps;
double prec;

// thermodynamical state
int npart;
double beta,temp,energy,vol,rho,box,rcut;

// simulation
int iNVET, nstep, nblk, restart, save, verbose;
double delta;
string state, dir;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Initial_Conf(void);
void Equilibration(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
void Compute_g(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);
double Force(int, int);

// vector functions
double get_mean(vector<double> v);
vector<double> operator+ (const vector<double>& v1, const vector<double>& v2);
vector<double> operator/ (const vector<double>& v1, const vector<double>& v2);
vector<double> operator/ (const vector<double>& v1, double k);
vector<double> operator^ (const vector<double>& v1, double esp);


#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
