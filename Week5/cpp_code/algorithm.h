#ifndef __algorithm_h__
#define __algorithm_h__

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <vector>

#include "utils.h"
#include "random.h"
#include "functions.h"



class algorithm {

    public:

        algorithm(scalar_function* pdf);
        algorithm(scalar_function* pdf, double par);
        algorithm(scalar_function* pdf, double par, string T);
        ~algorithm();

        vector<double> change(vector<double>& config, Random& rnd);
        vector<double> equilibration(vector<double>& config, Random& rnd, double prec, bool save);
        double metropolis(vector<double>& config, Random& rnd, int steps, bool save);

        string get_T() const {return m_T;};
        double get_par() const {return m_par;};
        scalar_function* get_pdf() const {return m_pdf;};
        vector<vector<double>> get_configs() const {return m_configs;};

        void set_T(string T) {m_T = T;};
        void set_par(double par) {m_par = par;};
        void set_pdf(scalar_function* pdf) {m_pdf = pdf;};

        void save(string filename);


    protected:

        int m_proposed, m_accepted;

        string m_T;
        double m_par;
        scalar_function* m_pdf;
        
        vector<vector<double>> m_configs;
        
};



#endif
