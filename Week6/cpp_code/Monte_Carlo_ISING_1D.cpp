/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
  
  Input(); //Inizialization
  Equilibration(); //Initial equilibration

  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput, Primes, Seed;

  //Read input informations
  ReadInput.open("input.dat");

  ReadInput >> verbose;

  if(verbose)
  {
    cout << "Classic 1D Ising model             " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl;
  }

  ReadInput >> restart;
  if(restart) cout << "Restarting from the last saved configuration." << endl;

  ReadInput >> save;

    //Read seed for random numbers
    int p1, p2;
    Primes.open("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    if(restart) Seed.open("seed.out");
    else Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    Seed.close();


  ReadInput >> temp;
  beta = 1.0/temp;
  if(verbose) cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  if(verbose) cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  if(verbose) cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  if(verbose) cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> prec;

  ReadInput >> eqsteps;

  if(verbose)
  {
    if(metro==1) cout << "The program perform Metropolis moves" << endl;
    else cout << "The program perform Gibbs moves" << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
  } 
  
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
  n_props = 4; //Number of observables

  //Restarting from the last saved configuration
  if(restart)
  {
    ifstream ReadConf;
    ReadConf.open("config.out");
    for (int i=0; i<nspin; ++i) ReadConf >> s[i];
    ReadConf.close();
  }
  //Setting up the initial configuration
  else
  {
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and magnetization
  if(verbose)
  {
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
    cout << "Initial magnetization = " << walker[im] << endl;
  }
}


void Equilibration(void)
{

  int check = 0, nchecks = 0;
  vector<double> ms;
  double mean = 0., new_mean = 0.;

  if(verbose)
  {
    cout << endl;
    cout << "Initiating equilibration of the system" << endl;
  }

  for(int i=0; i<eqsteps; i++)
  {
    Move(metro);
    Measure();
		ms.push_back(walker[im]);
	}

  mean = get_mean(ms);
  nchecks = eqsteps;

  do
  {
		for(int i=0; i<int(eqsteps/4); i++)
    {
      Move(metro);
      Measure();
      ms.push_back(walker[im]);
      ms.erase(ms.begin());
		}

    new_mean = get_mean(ms);	
	  if( new_mean-mean < prec*mean ) check ++;
		else check = 0;

    mean = new_mean;
		nchecks += int(eqsteps/4);
	}while(check < 10);

  if(verbose)
  {
    cout << "The system is now equilibrated. Number of steps = " << nchecks << endl;
    cout << "Post-equilibration energy = " << walker[iu]/(double)nspin << endl;
    cout << "Post-equilibration magnetization = " << walker[im] << endl;
  }

}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new;
  double energy_up, energy_down, p_marg;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    p = rnd.Rannyu();

    if(metro==1) //Metropolis
    {
      energy_old = Boltzmann(s[o], o);
      energy_new = Boltzmann(-s[o], o);
      //cout << o << "  " << energy_new - energy_old << "   " << p << "  " << exp(-beta*(energy_new-energy_old)) << endl;
      if( p < exp(-beta*(energy_new-energy_old)) ) s[o] *= -1;
    }
    else //Gibbs sampling
    {
      energy_up = Boltzmann(+1, o);
      energy_down = Boltzmann(-1, o);
      p_marg = 1./(1. + exp(-beta*(energy_up-energy_down)));
      if( p <= p_marg) s[o] = -1;
      else s[o] = +1;
    }
  }
}


double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}


void Measure()
{
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    m += s[i];
  }
  walker[iu] = u;
  walker[ic] = walker[iu]*walker[iu];    // calcolo incompleto, salvo solo il valore di H^2, dopodichè effettuo la stima di C(N, T) in Averages()
  walker[im] = m/(double)nspin;
  walker[ix] = walker[im]*walker[im];    // come sopra, calcolo incompleto
}


double get_mean(vector<double>& v)
{
  double mean = 0;
  for(auto& el : v) mean += el;
  return mean/v.size();
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
  ofstream Ene, Heat, Mag, Chi;
  //const int wd=12;

    if(iblk == 1 && verbose) cout << endl << "----------------------------" << endl << endl;
    
    if(verbose) cout << "Block number " << iblk << endl;
    //cout << "Acceptance rate " << accepted/attempted << endl << endl;

    if(save == 0 && iblk == 1)
    {
      Ene.open("output.ene.dat", ios::out);
      Heat.open("output.heat.dat", ios::out);
      Mag.open("output.mag.dat", ios::out);
      Chi.open("output.chi.dat", ios::out);

      Ene << "Block" << "\t" << "Mean" << "\t" << "Prog_mean" << "\t" << "Error" << endl;
      Heat << "Block" << "\t" << "Mean" << "\t" << "Prog_mean" << "\t" << "Error" << endl;
      Mag << "Block" << "\t" << "Mean" << "\t" << "Prog_mean" << "\t" << "Error" << endl;
      Chi << "Block" << "\t" << "Mean" << "\t" << "Prog_mean" << "\t" << "Error" << endl;
    }
    else
    {
      Ene.open("output.ene.dat", ios::app);
      Heat.open("output.heat.dat", ios::app);
      Mag.open("output.mag.dat", ios::app);
      Chi.open("output.chi.dat", ios::app);
    }
    
    // Internal energy
    stima_u = blk_av[iu]/blk_norm/(double)nspin;
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << iblk <<  "\t" << stima_u << "\t" << glob_av[iu]/(double)iblk << "\t" << err_u << endl;
    Ene.close();
    
    // Heat capacity
    stima_u *= (double)nspin;
    stima_c = beta*beta * (blk_av[ic]/blk_norm - stima_u*stima_u) / nspin;
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << iblk <<  "\t" << stima_c << "\t" << glob_av[ic]/(double)iblk << "\t" << err_c << endl;
    Heat.close();
    
    // Magnetization
    stima_m = blk_av[im]/blk_norm;
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << iblk <<  "\t" << stima_m << "\t" << glob_av[im]/(double)iblk << "\t" << err_m << endl;
    Mag.close();

    // Magnetic susceptibility
    stima_x = nspin * beta * blk_av[ix]/blk_norm;
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << iblk <<  "\t" << stima_x << "\t" << glob_av[ix]/(double)iblk << "\t" << err_x << endl;
    Chi.close();

    if(verbose) cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  if(verbose) cout << "Print final configuration to file config.out " << endl << endl;
  WriteConf.open("config.out");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
