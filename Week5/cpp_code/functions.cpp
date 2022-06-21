#include "functions.h"


// ************ classi figlie: scalar_function ************
double psi_100::eval (vector<double> x) const {

    double pi = acos(-1);
    double r = get_mod(x);
    return pow(1./sqrt(pi) * exp(-r), 2);

}


double psi_210::eval (vector<double> x) const {

    double pi = acos(-1);
    double r = get_mod(x);
    return pow(1./8. * sqrt(2./pi) * x[2] * exp(-r/2.), 2);

}
