#include "algorithm.h"

// --------------------------------------------------------------
// -------------------- MC_algorithm METHODS --------------------
// --------------------------------------------------------------

// COSTRUTTORI E DISTRUTTORE
MC_algorithm::MC_algorithm(scalar_function* pdf)
{
    m_pdf = pdf;
    m_x_width = 1.;
	m_mc_steps = 1;
}

MC_algorithm::MC_algorithm(scalar_function* pdf, double x_width, int mc_steps)
{
    m_pdf = pdf;
    m_x_width = x_width;
	m_mc_steps = mc_steps;
}

MC_algorithm::~MC_algorithm() {}


// METODI
double MC_algorithm::step(double& x, Random& rnd)
{
	double new_x = rnd.Rannyu(x-m_x_width, x+m_x_width);
	double ratio = m_pdf->eval(new_x)/m_pdf->eval(x);
	double random = rnd.Rannyu();

	if(random <= ratio) return new_x;
	else return x;
}


double MC_algorithm::get_HT(double& x)
{
    double mu = m_pdf->get_params()[0], sigma = m_pdf->get_params()[1];
    double Norm = 2. * sqrt(M_PI*sigma*sigma) * (1. + exp(-pow(mu/sigma, 2.)));
    double xplus = (x+mu)/sigma, xminus = (x-mu)/sigma;
	double psiT = 1./sqrt(Norm) * ( exp(-xplus*xplus/2.) + exp(-xminus*xminus/2.) );
	
	double K_psiT = 1./(2.*sqrt(Norm)*sigma*sigma) * ( (1. - xplus*xplus)*exp(-xplus*xplus/2.) + (1. - xminus*xminus)*exp(-xminus*xminus/2.) );
	double V_psiT = ( pow(x, 4.) - 5./2. * pow(x, 2.) ) * psiT;
	
	return (K_psiT + V_psiT)/psiT;
}


double MC_algorithm::metropolis(double& x, Random& rnd, bool clear)
{
    double new_x;
	double meas = 0;

	if(clear) m_xs.clear();

	m_xs.push_back(x);

	for(int i=0; i<m_mc_steps; i++)
	{
		new_x = step(x, rnd);
		meas = meas + get_HT(new_x);
		x = new_x;
		m_xs.push_back(x);
	}

	return meas/m_mc_steps;
}


void MC_algorithm::save_xs(string filename)
{
	if(m_xs.size() == NULL)
		cout << "Errore: nessuna posizione salvata." << endl;

	else
	{
		int step = 0;
		ofstream fout(filename);
	
		fout << "step\tx" << endl;
		for(auto& x : m_xs)
		{
			fout <<	step << "\t" << x << endl;
			step ++;
		}
		
		fout.close();
	}
}











// --------------------------------------------------------------
// -------------------- SA_algorithm METHODS --------------------
// --------------------------------------------------------------

// COSTRUTTORI E DISTRUTTORE
SA_algorithm::SA_algorithm(scalar_function* pdf) : 
MC_algorithm(pdf) 
{
	m_mu_width = 1.;
    	m_sigma_width = 1.;
	m_sa_steps = 1;
	m_nblk = 1;
}

SA_algorithm::SA_algorithm(scalar_function* pdf, double mu_width, double sigma_width, int sa_steps, int nblk) : 
MC_algorithm(pdf)
{
    	m_mu_width = mu_width;
    	m_sigma_width = sigma_width;
	m_sa_steps = sa_steps;
	m_nblk = nblk;
}

SA_algorithm::SA_algorithm(scalar_function* pdf, double mu_width, double sigma_width, double x_width, int sa_steps, int mc_steps, int nblk) : 
MC_algorithm(pdf, x_width, mc_steps)
{
    	m_mu_width = mu_width;
    	m_sigma_width = sigma_width;
	m_sa_steps = sa_steps;
	m_nblk = nblk;
}

SA_algorithm::~SA_algorithm() {}


// METODI
vector<double> SA_algorithm::annealing(double& x, vector<double> Ts, Random& rnd, stats& gen, bool clear)
{

	double ratio, random;
	vector<double> blk_HT;
	vector<double> old_HT, new_HT, old_params, new_params;
	vector<double> best_params, best_HT;

	if(clear)
	{
		m_HTs.clear();
		m_params.clear();
	}

	// inizializzazione parametri
	old_params = m_pdf->get_params();
	best_params = old_params;
	m_params.push_back(old_params);

	// inizializzazione HT
	for(int blk=0; blk<m_nblk; blk++)  blk_HT.push_back( metropolis(x, rnd, 1) );
	gen.block_mean(blk_HT);
	old_HT = {gen.get_last_mean(), gen.get_last_err()};
	best_HT = old_HT;
	m_HTs.push_back(old_HT);

	// ciclo su T1, ... TN
	for(auto& T : Ts)
	{

		// # passi (m_sa_steps) uguali per ogni Ti
		for(int sa_step=0; sa_step<m_sa_steps; sa_step++)
		{

			// estrazione nuovi parametri
			new_params.push_back( abs(rnd.Rannyu(old_params[0]-m_mu_width, old_params[0]+m_mu_width)) );	// new_mu
			new_params.push_back( rnd.Rannyu(old_params[1]-m_sigma_width, old_params[1]+m_sigma_width) );	// new_sigma
			m_pdf->set_params(new_params);

			// calcolo HT con nuovi parametri
			blk_HT.clear();
			for(int blk=0; blk<m_nblk; blk++)  blk_HT.push_back( metropolis(x, rnd, 1) );
			gen.block_mean(blk_HT);
			new_HT = {gen.get_last_mean(), gen.get_last_err()};

			// aggiornamento best_params
			if(new_HT[0] <= best_HT[0])
			{
				best_HT = new_HT;
				best_params = new_params;
			}

			// shift parametri
			random = rnd.Rannyu();
			ratio = exp(-(new_HT[0]-old_HT[0])/T);
			if(random <= ratio)
			{
				old_params = new_params;
				old_HT = new_HT;
			}
			
			m_params.push_back(old_params);
			new_params.clear();
				
			m_HTs.push_back(old_HT);
			new_HT.clear();


			// aggiornamento a schermo
			if(sa_step%(m_sa_steps/5) == 0)
			{
				if(sa_step == 0 && get_index(Ts, T) == 0) cout << "----------------------------------------" << endl;;
				cout << "T = " << T << ", step " << sa_step << "/" << m_sa_steps << endl;
				cout << "mu = " << best_params[0] << ", 	sigma = " << best_params[1] << endl;
				cout << "<HT> = " << best_HT[0] << " +- " << best_HT[1] << endl;
				cout << "----------------------------------------" << endl;
			}


		}	

	}

	return best_params;

}



void SA_algorithm::save_HTs(string filename)
{
	if(m_HTs.size() == NULL)
		cout << "Errore: nessun valore di energia salvato." << endl;

	else
	{
		int step = 0;
		ofstream fout(filename);
	
		fout << "step\tHT\terr" << endl;
		for(auto& vec : m_HTs)
		{
			fout <<	step << "\t" << vec[0] << "\t" << vec[1] << endl;
			step ++;
		}
		
		fout.close();
	}
}



void SA_algorithm::save_params(string filename)
{
	if(m_params.size() == NULL)
		cout << "Errore: nessuna sequenza di parametri salvata." << endl;

	else
	{
		int step = 0;
		ofstream fout(filename);
	
		fout << "step\tmu\tsigma" << endl;
		for(auto& vec : m_params)
		{
			fout <<	step << "\t" << vec[0] << "\t" << vec[1] << endl;
			step ++;
		}
		
		fout.close();
	}
}
