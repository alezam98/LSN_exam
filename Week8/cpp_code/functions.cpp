#include "functions.h"


// ************ classi figlie: scalar_function ************
psiT2::psiT2()
{
    m_params.push_back(0.);
    m_params.push_back(1.);
}

psiT2::psiT2(vector<double> params)
{
    if(params.size() != 2)
    {
        cout << "Errore: la dimensione del vettore deve essere 2 (media e deviazione std). Esco!" << endl;
        exit(1);
    }
    else m_params = params;
}

psiT2::psiT2(double mu, double sigma)
{
    m_params.push_back(mu);
    m_params.push_back(sigma);
}

psiT2::~psiT2() {}



double psiT2::eval (double x) const 
{
    double mu = m_params[0], sigma = m_params[1];
    double Norm = 2. * sqrt(M_PI*sigma*sigma) * (1. + exp(-pow(mu/sigma, 2.)));
    double xplus = (x+mu)/sigma, xminus = (x-mu)/sigma;
    return 1./Norm * pow(exp(-xplus*xplus/2.) + exp(-xminus*xminus/2.), 2.);
}