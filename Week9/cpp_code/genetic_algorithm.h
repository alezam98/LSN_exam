#ifndef __genetic_algorithm_h__
#define __genetic_algorithm_h__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

#include "random.h"
#include "utils.h"

using namespace std;



// function for chromosome sorting
int compare(const void *aa, const void *bb);

// struct "chromosome"
typedef struct _chromosome
{
    vector<double> sequence;
    double fitness;
} chromosome;


// class "genetic_algorithm"
class genetic_algorithm
{

    public:

        genetic_algorithm();
        ~genetic_algorithm();

        void initialize(Random& rnd);
        int* selection(Random& rnd);
        void crossover(int update_index, Random& rnd);
        void mutation(int update_index, Random& rnd);
        void update();

        double fitness(vector<double> sequence) const;
        void check() const;
        void clear() const;

        void set_nchromo(int nchromo);
        void set_cities(vector<vector<double>> cities) {m_cities = cities;};
        void set_pc(double pc) {m_pc = pc;};
        void set_pms(vector<double> pms) {m_pms = pms;};

        vector<double> get_L() const {return m_L;};
        
        void save_chromosome(string filename) const;
        void save_cities(string filename) const;
        void save_fitness(string filename) const;

        void print_generation(int igen) const;
        void print_best(int igen) const;

    protected:

        int m_nchromo;
        chromosome* m_population;
        chromosome* m_new_population;
        vector<vector<double>> m_cities;

        vector<double> m_L;
        vector<double> m_best_fitness;
        chromosome m_best_chromosome;

        double m_pc;
        vector<double> m_pms;

};



# endif