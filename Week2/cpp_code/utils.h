#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>

#include "random.h"

using namespace std;


double get_mod(vector<double>& v);
double get_sum(vector<double>& v);
double get_mean(vector<double>& v);
double get_var(vector<double>& v);
double get_var(vector<double>& v, double m);
double get_dev_std(vector<double>& v);
double get_dev_std(vector<double>& v, double m);
double get_chi_sq(vector<double>& v, double m, double x_min, double x_max, int n_int);
void save(vector<double>& v, string filename, string head);
void save(vector<double>& v1, vector<double>& v2, string filename, string head1, string head2);

double integranda (double x);
double integrale_media (int n, Random & rnd);
double integrale_sample (int n, Random & rnd);
double integrale_ar (int n, Random & rnd);

vector<double> operator+ (const vector<double>& v1, const vector<double>& v2);
vector<double> operator+ (const vector<double>& v1, double k);
vector<double> operator- (const vector<double>& v1, const vector<double>& v2);
vector<double> operator- (const vector<double>& v1, double k);
double operator* (const vector<double>& v1, const vector<double>& v2);
vector<double> operator* (double k, const vector<double>& v1);
vector<double> operator/ (const vector<double>& v1, const vector<double>& v2);
vector<double> operator/ (const vector<double>& v1, double k);
vector<double> operator^ (const vector<double>& v1, double k);
