#ifndef __rwalk_h__
#define __rwalk_h__

#include <iostream>
#include <cmath>
#include <vector>

#include "random.h"
#include "utils.h"

using namespace std;


class rwalk {

  public:

	rwalk();
	rwalk(double a);
	~rwalk();
	
	void set_a(double a) {m_a = a;};
	
	double get_a() {return m_a;};
	vector<double*> get_rw() {return m_rw;};
	vector<double> get_rN2() {return m_rN2;};

	void step(Random& rnd, bool continous);
	void step(Random& rnd, bool continous, int n_step);
	vector<double> sigma2(Random& rnd, int n_step, int n_rw, bool continous);
	
	void clear();
	void print();

  protected:

	double m_a;
	vector<double*> m_rw;
	vector<double> m_rN2;

};


#endif
