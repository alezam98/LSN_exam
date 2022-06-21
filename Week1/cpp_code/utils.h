#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;


double get_sum(vector<double>& v);
double get_mean(vector<double>& v);
double get_var(vector<double>& v);
double get_var(vector<double>& v, double m);
double get_dev_std(vector<double>& v);
double get_dev_std(vector<double>& v, double m);
double get_chi_sq(vector<double>& v, double m, double x_min, double x_max, int n_int);
void save(vector<double>& v, string filename, string head);

vector<double> operator+ (const vector<double>& v1, const vector<double>& v2);
vector<double> operator* (double k, const vector<double>& v1);
double operator* (const vector<double>& v1, const vector<double>& v2);
