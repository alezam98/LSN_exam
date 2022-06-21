#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "mpi.h"
#include "random.h"
#include "utils.h" 
#include "stats.h"
#include "genetic_algorithm.h"

using namespace std;


// Random numbers
#include "random.h"
int seed[4];
Random rnd;

// simulation
int nchromo, ncities, ngen, nmigr;
int square, circle;
int verbose, restart;
double exchange_fract;
string check, dir, shape, city_type;

// cities
vector<vector<double>> cities;

// evolution
int update_index;
double pc, pm_1, pm_2, pm_3, pm_4;
vector<double> pms;
vector<double> mean_L;

// classi
stats gen;
genetic_algorithm* algorithm = new genetic_algorithm();

//functions
void Input(int rank);
void Migration(int size, int rank);
void Send_Receive(chromosome& chromo, int rank1, int rank2, int rank);
void Save(int rank);