#ifndef __stats_h__
#define __stats_h__

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>

using namespace std;


class stats {

  public:
 	
	stats();
	~stats();
	
	void block_mean(vector<double>& v);
	
	vector<double> get_mean() const {return m_mean_prog;};
	vector<double> get_err() const {return m_err_prog;};
	
	double get_last_mean() const {return m_mean_prog[m_mean_prog.size()-1];};
	double get_last_err() const {return m_err_prog[m_err_prog.size()-1];};
	
	void save(string filename);
	
	void clear();
	

  protected:
  
  	vector<double> m_mean_prog;
	vector<double> m_err_prog;
 	
};


#endif
