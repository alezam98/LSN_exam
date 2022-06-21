#ifndef __functions_h__
#define __functions_h__

#include <cmath>
#include <cstdlib>
#include <vector>

#include "utils.h"
#include "random.h"

using namespace std;



// ************ classi madri ************
class scalar_function {

	public:
	
		virtual double eval(double x) const =0;
		vector<double> get_params() const {return m_params;};
		void set_params(vector<double> params) {m_params = params;};

	protected:

		vector<double> m_params;
		
};



// ************ classi figlie: scalar_function ************
class psiT2: public scalar_function {

	public:

		psiT2();
		psiT2(vector<double> params);
		psiT2(double mu, double sigma);
		~psiT2();

		double eval(double x) const;

		double get_mu() const {return m_params[0];};
		double get_sigma() const {return m_params[1];};

		void set_mu(double mu) {m_params[0] = mu;};
		void set_sigma(double sigma) {m_params[1] = sigma;};

};



#endif