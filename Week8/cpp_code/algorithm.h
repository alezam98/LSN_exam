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
#include "stats.h"


// sampling Monte Carlo
class MC_algorithm 
{

    public:

        MC_algorithm(scalar_function* pdf);
        MC_algorithm(scalar_function* pdf, double x_width, int mc_steps);
        ~MC_algorithm();

        double step(double& x, Random& rnd);
        double get_HT(double& x);
        double metropolis(double& x, Random& rnd, bool clear);

        void set_mc_steps(double mc_steps) {m_mc_steps = mc_steps;};
        void set_x_width(double x_width) {m_x_width = x_width;};
        void set_pdf(scalar_function* pdf) {m_pdf = pdf;};

        int get_mc_steps() const {return m_mc_steps;};
        double get_x_width() const {return m_x_width;};
        scalar_function* get_pdf() const {return m_pdf;};

        vector<double> get_xs() const {return m_xs;};

        void save_xs(string filename);

    protected:

        int m_mc_steps;
        double m_x_width;
        scalar_function* m_pdf;

        vector<double> m_xs;
        
};



// SA annealing
class SA_algorithm : public MC_algorithm
{

    public:

        SA_algorithm(scalar_function* pdf);
        SA_algorithm(scalar_function* pdf, double mu_width, double sigma_width, int sa_steps, int nblk);
        SA_algorithm(scalar_function* pdf, double mu_width, double sigma_width, double x_width, int sa_steps, int mc_steps, int nblk);
        ~SA_algorithm();

        vector<double> annealing(double& x, vector<double> Ts, Random& rnd, stats& gen, bool clear);

        void set_nblk(int nblk) {m_nblk = nblk;};
        void set_sa_steps(int sa_steps) {m_sa_steps = sa_steps;};
        void set_mu_width(double mu_width) {m_mu_width = mu_width;};
        void set_sigma_width(double sigma_width) {m_sigma_width = sigma_width;};

        int get_nblk() const {return m_nblk;};
        int get_sa_steps() const {return m_sa_steps;};
        double get_mu_width() const {return m_mu_width;};
        double get_sigma_width() const {return m_sigma_width;};

        vector<vector<double>> get_HTs() const {return m_HTs;};
        vector<vector<double>> get_params() const {return m_params;};

        void save_HTs(string filename);
        void save_params(string filename);

    protected:

        int m_nblk, m_sa_steps;
        double m_mu_width, m_sigma_width;

        vector<vector<double>> m_HTs;
        vector<vector<double>> m_params;

};



#endif