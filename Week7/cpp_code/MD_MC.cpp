/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "MD_MC.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  Equilibration();

  int nconf = 1;
  for(int iblk=1; iblk <= nblk; iblk++) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();
      Measure();
      Compute_g();
      Accumulate(); //Update block averages
      if(istep%10 == 0){
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
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
  ReadInput.open("input.in");

  ReadInput >> iNVET;

  ReadInput >> verbose;

  ReadInput >> restart;

  if(verbose)
  {
    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;
  }

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

  ReadInput >> npart;
  if(verbose) cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  if(verbose) cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  if(verbose) cout << "Volume of the simulation box = " << vol << endl;
  if(verbose) cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  if(verbose) cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> MD_eqsteps;

  ReadInput >> MC_eqsteps;

  ReadInput >> prec;

  ReadInput >> r_bins;
  
  blk_g.assign(r_bins, 0.);
  glob_g.assign(r_bins, 0.);
  glob_g2.assign(r_bins, 0.);
  stima_g.assign(r_bins, 0.);
  for(int i=0; i<r_bins; i++)
    Norm.push_back(npart * rho * 4.*M_PI/3. * (pow((i+1)*box/2./r_bins, 3.) - pow(i*box/2./r_bins, 3.)));

  ReadInput >> save;

  ReadInput >> state;

  ReadInput >> dir;

  if(verbose)
  {
    cout << "The program perform Metropolis moves with uniform translations" << endl;
    cout << "Moves parameter = " << delta << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
  }
  ReadInput.close();

//Prepare arrays for measurements
  iv = 0; //Potential energy
  it = 1; //Temperature
  ik = 2; //Kinetic energy
  ie = 3; //Total energy
  ip = 4; //Pressure
  n_props = 5; //Number of observables

//Read initial configuration
  if(verbose) cout << "Read initial configuration" << endl << endl;
  Initial_Conf();

//Evaluate properties of the initial configuration
  Measure();

//Print initial values for measured properties
if(verbose)
{
  cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
  cout << "Initial kinetic energy   = " << walker[ik]/(double)npart << endl;
  cout << "Initial total energy     = " << walker[ie]/(double)npart << endl;
  cout << "Initial temperature      = " << walker[it] << endl;
  cout << "Initial pressure         = " << walker[ip] << endl;
}

  return;
}



void Initial_Conf(void)
{

  ifstream ReadConf, ReadVelocity;

  if(restart)
  {
    ReadConf.open("config.out");
    ReadVelocity.open("velocity.out");
    for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
  }
  else 
  {
    ReadConf.open("config.in");
    //if(verbose) cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i)
    {
      vx[i] = rnd.Gauss(0.,sqrt(temp));
      vy[i] = rnd.Gauss(0.,sqrt(temp));
      vz[i] = rnd.Gauss(0.,sqrt(temp));
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i)
    {
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    //if(verbose) cout << "velocity scale factor: " << fs << endl << endl;
    for (int i=0; i<npart; ++i)
    {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }

  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];     // le coordinate sono misurate in unita' di misura del lato
    x[i] = Pbc( x[i] * box );             // Pbc: Periodic Boudary Condition (non dovrebbe servire, scrupolo)
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();

  for (int i=0; i<npart; ++i)
  {
    if(iNVET)
    {
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else
    {
      xold[i] = Pbc(x[i] - vx[i] * delta);    // uso Verlet, non usa velocita', mi servono posizioni vecchie (primo passo, uso Eulero con dt = -delta)
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }

}



void Equilibration(void)
{
    
  if(verbose) cout << "Initiating equilibration of the system" << endl;

  // Equilibration for MC simulation
  if(iNVET)
  {

    int check = 0, nchecks = 0;
    vector<double> vs;
    double mean, new_mean;

    for(int i=0; i<MC_eqsteps; i++)
    {
      Move();
      Measure();
      vs.push_back(walker[iv]);
    }

    mean = get_mean(vs);
    nchecks = MC_eqsteps;

    do
    {
      for(int i=0; i<int(MC_eqsteps/4); i++)
      {
        Move();
        Measure();
        vs.push_back(walker[iv]);
        vs.erase(vs.begin());
      }

      new_mean = get_mean(vs);	
      if( abs(new_mean-mean) < abs(mean*prec) ) check ++;
      else check = 0;

      mean = new_mean;
      nchecks += int(MC_eqsteps/4);
    }while(check < 10);

    if(verbose)
    {
      cout << "The system is now equilibrated. Number of steps = " << nchecks << endl;
      cout << "Post-equilibration potential energy = " << walker[iv]/(double)npart << endl;
    }
  
  }

  // Equilibration for MD simulation
  else
  {

    int check = 0;
    double temp_goal = temp;
    double temp_f = 0.0;
    double kin = 0.0;

    do
    {

      for(int i=0; i<MD_eqsteps; i++) Move();

      temp_f = 0.;
      for(int i=0; i<MD_eqsteps/4.; i++){
        Move();
        kin = 0.;
        for (int i=0; i<npart; i++) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
        temp_f += (2.0 / 3.0) * kin/(double)npart;
      }
      temp_f /= MD_eqsteps/4.;

      if(abs(temp_f-temp_goal) < prec)
      {
        check = 1;
        temp = temp_goal;
      }
      else
      {
        temp += (temp_goal-temp_f)/2.;
        beta = 1/temp;
        Initial_Conf();
      }

    }while(check == 0);

    if(verbose)
    {
      cout << "The system is now equilibrated." << endl;
      cout << "Post-equilibration temperature = " << temp_f << endl;
    }
    
  }

}



void Move()
{
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) // Monte Carlo (NVT) move
  {
    for(int i=0; i<npart; ++i)
    {
    //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
      o = (int)(rnd.Rannyu()*npart);

    //Old
      energy_old = Boltzmann(x[o],y[o],z[o],o);

    //New
      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

      energy_new = Boltzmann(x[o],y[o],z[o],o);

    //Metropolis test
      p = exp(beta*(energy_old-energy_new));
      if(p >= rnd.Rannyu())  
      {
      //Update
        xold[o] = x[o];
        yold[o] = y[o];
        zold[o] = z[o];
        accepted = accepted + 1.0;
      } else {
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted = attempted + 1.0;
    }
  } else // Molecular Dynamics (NVE) move
  {
    double fx[m_part], fy[m_part], fz[m_part];

    for(int i=0; i<npart; ++i){ //Force acting on particle i
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }

    for(int i=0; i<npart; ++i){ //Verlet integration scheme

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;
      y[i] = ynew;
      z[i] = znew;

      accepted = accepted + 1.0;        // variabili sensate solo nel MC (qua accetto tutte le mosse)
      attempted = attempted + 1.0;
    }
  }
  return;
}



void Compute_g(void)
{

  int k;
  double dx, dy, dz, d_ij;

  for(int i=0; i<npart-1; i++)
  {
    for(int j=i+1; j<npart; j++)
    {
      dx = Pbc( x[i]-x[j] );
      dy = Pbc( y[i]-y[j] );
      dz = Pbc( z[i]-z[j] );
      d_ij = sqrt( dx*dx + dy*dy + dz*dz );

      if(d_ij/box < 1)
      {
        k = int( r_bins * d_ij/box );
        blk_g[k] += 2.;
      }
    }
  }

}



double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
 
  return f;
}

void Measure() //Properties measurement
{
  double v = 0.0, kin=0.0, w=0.0;
  double vij, wij;
  double tail_v, tail_w;
  double dx, dy, dz, dr;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {
// distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        v += vij;

        wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
        w += wij;
      }
    }          
  }

  for (int i=0; i<npart; ++i) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  // tail corrections
  tail_v = 8. * M_PI * npart * rho / 9. * (1./pow(rcut, 9.) - 3./pow(rcut, 3.));
  tail_w = 32. * M_PI * npart * rho * (1./(3. * pow(rcut, 9.)) - 1./(2. * pow(rcut, 3.)));

  walker[iv] = 4.0 * v + tail_v; // Potential energy
  walker[ik] = kin; // Kinetic energy
  walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperature
  walker[ie] = 4.0 * v + kin;  // Total energy;
  walker[ip] = rho * walker[it] + (1./3./vol) * 48.0 * (w + tail_w); // Pressure 

  return;
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
   
  blk_g.clear();
  blk_g.assign(r_bins, 0.);
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
    
   ofstream Epot, Ekin, Etot, Temp, Press, Rdf;
   //const int wd=12; usata per set(wd)

    if(verbose)
    {
      if(iblk == 1) cout << endl << "----------------------------" << endl << endl;
    
      cout << "Block number " << iblk << endl;
      cout << "Acceptance rate " << accepted/attempted << endl << endl;
    }
    
    
    if(save == 0 && iblk == 1)
    {
      Epot.open("../data/"+dir+"/"+state+"/output_epot.dat",ios::out);
      Ekin.open("../data/"+dir+"/"+state+"/output_ekin.dat",ios::out);
      Temp.open("../data/"+dir+"/"+state+"/output_temp.dat",ios::out);
      Press.open("../data/"+dir+"/"+state+"/output_press.dat",ios::out);
      Etot.open("../data/"+dir+"/"+state+"/output_etot.dat",ios::out);
      Rdf.open("../data/"+dir+"/"+state+"/output_rdf.dat",ios::out);

      Epot << "Block" << "\t" << "Mean" << "\t" << "Prog_mean" << "\t" << "Error" << endl;
      Ekin << "Block" << "\t" << "Mean" << "\t" << "Prog_mean" << "\t" << "Error" << endl;
      Temp << "Block" << "\t" << "Mean" << "\t" << "Prog_mean" << "\t" << "Error" << endl;
      Press << "Block" << "\t" << "Mean" << "\t" << "Prog_mean" << "\t" << "Error" << endl;
      Etot << "Block" << "\t" << "Mean" << "\t" << "Prog_mean" << "\t" << "Error" << endl;
      Rdf << "r" << "\t" << "<g>" << "\t" << "err_g" << endl;
    } else 
    {
      Epot.open("../data/"+dir+"/"+state+"/output_epot.dat",ios::app);
      Ekin.open("../data/"+dir+"/"+state+"/output_ekin.dat",ios::app);
      Temp.open("../data/"+dir+"/"+state+"/output_temp.dat",ios::app);
      Press.open("../data/"+dir+"/"+state+"/output_press.dat",ios::app);
      Etot.open("../data/"+dir+"/"+state+"/output_etot.dat",ios::app);
      if(iblk == nblk) Rdf.open("../data/"+dir+"/"+state+"/output_rdf.dat",ios::app);
    }
    
    //Potential energy
    stima_pot = blk_av[iv]/blk_norm/(double)npart; 
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    Epot << iblk <<  "\t" << stima_pot << "\t" << glob_av[iv]/(double)iblk << "\t" << err_pot << endl;
    
    //Kinetic energy
    stima_kin = blk_av[ik]/blk_norm/(double)npart; 
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);
    Ekin << iblk <<  "\t" << stima_kin << "\t" << glob_av[ik]/(double)iblk << "\t" << err_kin << endl;

    //Total energy
    stima_etot = blk_av[ie]/blk_norm/(double)npart; 
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);
    Etot << iblk <<  "\t" << stima_etot << "\t" << glob_av[ie]/(double)iblk << "\t" << err_etot << endl;

    //Temperature
    stima_temp = blk_av[it]/blk_norm; 
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);
    Temp << iblk <<  "\t" << stima_temp << "\t" << glob_av[it]/(double)iblk << "\t" << err_temp << endl;

    //Pressure
    stima_press = blk_av[ip]/blk_norm; 
    glob_av[ip] += stima_press;
    glob_av2[ip] += stima_press*stima_press;
    err_press=Error(glob_av[ip],glob_av2[ip],iblk);
    Press << iblk <<  "\t" << stima_press << "\t" << glob_av[ip]/(double)iblk << "\t" << err_press << endl;
    
    //Radial distribution function
    for(int i=0; i<r_bins; i++)
    {
      stima_g[i] = blk_g[i]/Norm[i]/blk_norm; 
      glob_g[i] += stima_g[i];
      glob_g2[i] += stima_g[i]*stima_g[i];
      //cout << Norm[i] << "\t" << stima_g[i] << endl;
    }

    if(iblk == nblk)
    {
      for(int i=0; i<r_bins; i++)
      {
        err_g = Error(glob_g[i],glob_g2[i],iblk);
        Rdf << i*box/2./r_bins <<  "\t" << glob_g[i]/iblk << "\t" << err_g << endl;
      }
    }

    if(verbose) cout << "----------------------------" << endl << endl;

    Epot.close();
    Ekin.close();
    Etot.close();
    Temp.close();
    Press.close();
    if(iblk == nblk) Rdf.close();
}


void ConfFinal(void)
{
  ofstream WriteConf, WriteVelocity, WriteSeed;

  if(verbose) cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open("config.out");
  WriteVelocity.open("velocity.out");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}





// VECTOR FUNCTIONS

// compute average
double get_mean(vector<double> v)
{
  double mean = 0.;
  for(auto& el : v) mean += el;
  return mean/v.size();
}

// sum between vectors
vector<double> operator+ (const vector<double>& v1, const vector<double>& v2){

	if(v1.size() != v2.size()){
		cout << "Errore: somma di due vettori di taglie diverse!" << endl;
		exit(0);
	}

	vector<double> v;
	for(int i=0; i<int(v1.size()); i++)
		v.push_back(v1[i] + v2[i]);

	return v;
	
}

// division between vectors
vector<double> operator/ (const vector<double>& v1, const vector<double>& v2){

	if(v1.size() != v2.size()){
		cout << "Errore: divisione tra due vettori di taglie diverse!" << endl;
		exit(0);
	}

	vector<double> v;
	for(int i=0; i<int(v1.size()); i++){
		
		if(v2[i] == 0){
			cout << "Errore: divisione per elemento nullo! (elemento " << i << "-esimo)" << endl;
			exit(0);
		}
		
		v.push_back( v1[i]/v2[i] );
		
	}

	return v;
	
}

// division by a constant
vector<double> operator/ (const vector<double>& v1, double k){

	vector<double> v;
	
	if(k == 0){
		cout << "Errore: la costante deve essere non nulla!" << endl;
		exit(0);
	}
	
	for(auto& el : v1)
		v.push_back(el / k);

	return v;
}

// exponentiation of a vector
vector<double> operator^ (const vector<double>& v1, double esp){

	vector<double> v;
	for(int i=0; i<int(v1.size()); i++)
		v.push_back(pow(v1[i], esp));

	return v;
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
