#include "algorithm.h"


// COSTRUTTORI E DISTRUTTORE
algorithm::algorithm(scalar_function* pdf){

    m_pdf = pdf;
    m_par = 1.;
    m_T = "uniform";
	m_proposed = 0;
	m_accepted = 0;

}

algorithm::algorithm(scalar_function* pdf, double par){

    m_pdf = pdf;
    m_par = par;
    m_T = "uniform";
	m_proposed = 0;
	m_accepted = 0;

}

algorithm::algorithm(scalar_function* pdf, double par, string T){

    m_pdf = pdf;
    m_par = par;
    m_T = T;
	m_proposed = 0;
	m_accepted = 0;

}

algorithm::~algorithm(){
}


// METODI
vector<double> algorithm::change(vector<double>& config, Random& rnd){

	vector<double> new_config;
	if(m_T == "uniform")
	{
		for(auto& el : config) new_config.push_back( rnd.Rannyu(el-m_par, el+m_par) );
	}
	else if(m_T == "gaussian")
	{
		for(auto& el : config) new_config.push_back( rnd.Gauss(el, m_par) );
	}
	
	double ratio = m_pdf->eval(new_config)/m_pdf->eval(config);
	double random = rnd.Rannyu();
	if(random <= ratio) 
	{
		m_proposed ++;
		m_accepted ++;
		return new_config;
	}
	else 
	{
		m_proposed ++;
		return config;
	}

}


vector<double> algorithm::equilibration(vector<double>& config, Random& rnd, double prec, bool save)
{

	vector<double> new_config;
	vector<double> measure;
	double mean, new_mean;
	int check = 0;

	m_configs.clear();

	for(int i=0; i<1000; i++)
	{
		new_config = change(config, rnd);
		measure.push_back( get_mod(new_config) );
		config = new_config;
		if(save) m_configs.push_back(config);
	}

	mean = get_mean(measure);

	do
	{
		for(int i=0; i<100; i++){
			new_config = change(config, rnd);
			measure.push_back( get_mod(new_config) );
			measure.erase(measure.begin());
			config = new_config;
			if(save) m_configs.push_back(config);
		}

		new_mean = get_mean(measure);
			
		if( (new_mean-mean)/mean < prec ) check ++;
		else check = 0;

		mean = new_mean;
	}while(check < 10);

	return new_config;

}


double algorithm::metropolis(vector<double>& config, Random& rnd, int steps, bool save)
{
    
    vector<double> new_config;
	double measure = 0;

	m_configs.clear();
	m_proposed = 0;
	m_accepted = 0;

	if(save)
	{
		m_configs.push_back(config);

		for(int i=0; i<steps; i++)
		{
			new_config = change(config, rnd);
			measure = measure + get_mod(new_config);
			config = new_config;
			m_configs.push_back(config);
		}
	}

	else
	{
		for(int i=0; i<steps; i++)
		{
			new_config = change(config, rnd);
			measure = measure + get_mod(new_config);
			config = new_config;
		}
	}

	cout << "\tAcceptance rate: " << double(m_accepted)/double(m_proposed) << endl;

	return measure/steps;

}


void algorithm::save(string filename){

	if(m_configs.size() == NULL){
		cout << "Error (algorithm::save): no configuration saved." << endl;
	}

	else{

		ofstream fout(filename);
	
		fout << "x\ty\tz" << endl;
		for(auto& config : m_configs)
			fout << config[0] << "\t" << config[1] << "\t" << config[2] << endl;
		
		fout.close();
	}

}
