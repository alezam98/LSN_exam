#ifndef __functions_h__
#define __functions_h__

#include <cmath>
#include <cstdlib>
#include <vector>

#include "utils.h"
#include "random.h"

using namespace std;



// ************ classe madre ************
class scalar_function {

	public:
	
		virtual double eval(vector<double>) const =0;
		
};



// ************ classi figlie ************
class psi_100: public scalar_function {

	public:
	
		psi_100() {};
		~psi_100() {};
		
		virtual double eval(vector<double> x) const;
		
};


class psi_210: public scalar_function {

	public:
	
		psi_210() {};
		~psi_210() {};
		
		virtual double eval(vector<double> x) const;
		
};



#endif
