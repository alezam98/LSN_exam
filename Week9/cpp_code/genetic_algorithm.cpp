#include "genetic_algorithm.h"

// constructors/destructor
genetic_algorithm::genetic_algorithm() 
{
    m_nchromo = -1;
    m_pc = 0.;
}

genetic_algorithm::~genetic_algorithm() {}



// methods
void genetic_algorithm::initialize(Random& rnd)
{
    int ncities = m_cities.size();
    int city_pos;
    vector<int> basic_sequence;
    vector<int> copy;


    for(int icity=2; icity<=ncities; icity++) basic_sequence.push_back(icity);
    
    for(int ichromo=0; ichromo<m_nchromo; ichromo++)
    {
        m_population[ichromo].sequence.push_back(1);    // first city fixed (sequence[0] = 1)
        copy = basic_sequence;

        // initializing sequence
        for(int icity=2; icity<=ncities; icity++)
        {
            city_pos = int( rnd.Rannyu(0, copy.size()) );
            m_population[ichromo].sequence.push_back(copy[city_pos]);
            copy.erase(copy.begin() + city_pos);
        }

        // setting fitness
        m_population[ichromo].fitness = fitness(m_population[ichromo].sequence);
    }

    // check of the starting position
    check();

    // sorting by fitness, from worst to best
    qsort(m_population, m_nchromo, sizeof(chromosome), compare);
    
    // updating best_chromosome/fitness
    m_best_chromosome.sequence = m_population[0].sequence;
    m_best_chromosome.fitness = m_population[0].fitness;
    m_best_fitness.push_back(m_best_chromosome.fitness);
}



void genetic_algorithm::crossover(int update_index, Random& rnd)
{
    int ipar1, ipar2;
    int ncities = m_cities.size();
    chromosome parent1, parent2, offspring1, offspring2;
    vector<double> copy1, copy2;

    // selection with replacement
    do
    {
        ipar1 = int( m_nchromo * pow(rnd.Rannyu(), 2) );
        ipar2 = int( m_nchromo * pow(rnd.Rannyu(), 2) );
    } while(ipar1 == ipar2);

    // crossover
    parent1.sequence = m_population[ipar1].sequence;
    parent2.sequence = m_population[ipar2].sequence;

    offspring1.sequence = parent1.sequence;
    offspring2.sequence = parent2.sequence;

    if( rnd.Rannyu() <= m_pc )
    {
        int cross_cut = int(rnd.Rannyu(1, ncities));

        // offspring sequences
        for(int i=0; i<ncities; i++)
        {
            for(int j=cross_cut; j<ncities; j++)
            {
                if(parent2.sequence[i] == parent1.sequence[j])
                    copy1.push_back(parent2.sequence[i]);

                if(parent1.sequence[i] == parent2.sequence[j])
                    copy2.push_back(parent1.sequence[i]);  
            }
            if( int(copy1.size()) == ncities-cross_cut && int(copy2.size()) == ncities-cross_cut ) break;
        }

        // fitness
        offspring1.fitness = fitness(offspring1.sequence);
        offspring2.fitness = fitness(offspring2.sequence);
    }
    else
    {
        offspring1.fitness = m_population[ipar1].fitness;
        offspring2.fitness = m_population[ipar2].fitness;
    }

    m_new_population[update_index].sequence = offspring1.sequence;
    m_new_population[update_index].fitness = offspring1.fitness;

    m_new_population[update_index+1].sequence = offspring2.sequence;
    m_new_population[update_index+1].fitness = offspring2.fitness;
}



void genetic_algorithm::mutation(int update_index, Random& rnd)
{
    int ncities = m_cities.size();
    chromosome mutated;

    // for each offspring...
    for(int ioff=0; ioff<2; ioff++)
    {
        // old chromosome
        mutated.sequence = m_new_population[update_index+ioff].sequence;

        // mutation
        if(rnd.Rannyu() <= m_pms[0])               // pair swap
        {
            int pos1, pos2;
            do
            {
                pos1 = int( rnd.Rannyu(1, ncities-1) );
                pos2 = int( rnd.Rannyu(1, ncities-1) );
            } while(pos1 == pos2);
            swap(mutated.sequence, pos1, pos2);
        }

        else if(rnd.Rannyu() <= m_pms[1])          // shift
        {
            int shift;
            vector<double> copy = mutated.sequence;
            shift = int( rnd.Rannyu(1, ncities-1) );
            for(int i=1; i<ncities; i++)
            {
                if(i+shift < ncities) mutated.sequence[i+shift] = copy[i];
                else mutated.sequence[i+shift-ncities+1] = copy[i];
            }
        }

        else if(rnd.Rannyu() <= m_pms[2])          // contiguos swap
        {
            int pos, m, shift;
            vector<double> copy = mutated.sequence;
            pos = int( rnd.Rannyu(1, ncities-1) );
            m = int( rnd.Rannyu(1, int((ncities-pos)/2)) );
            do{
                shift = int( rnd.Rannyu(1-pos, ncities-pos-m+1) );
            } while(shift == 0 || abs(shift) < m);
            for(int i=0; i<m; i++) swap(mutated.sequence, pos+i, pos+i+shift);
        }

        else if(rnd.Rannyu() <= m_pms[3])          // inversion
        {
            int pos1, pos2;
            do
            {
                pos1 = int( rnd.Rannyu(1, ncities-1) );
                pos2 = int( rnd.Rannyu(1, ncities-1) );
            } while(pos1 == pos2);
            if(pos2 > pos1)
                for(int i=0; i<int((pos2-pos1+1)/2); i++) swap(mutated.sequence, pos1+i, pos2-i);
            else
                for(int i=0; i<int((pos1-pos2+1)/2); i++) swap(mutated.sequence, pos2+i, pos1-i);
        }

        // calculating new fitness
        mutated.fitness = fitness(mutated.sequence);

        // updating
        m_new_population[update_index+ioff].sequence = mutated.sequence;
        m_new_population[update_index+ioff].fitness = mutated.fitness;
    }
}



void genetic_algorithm::update()
{
    double L = 0.;

    // update population
    for(int ichromo=0; ichromo<m_nchromo; ichromo++)
    {
        m_population[ichromo].sequence = m_new_population[ichromo].sequence;
        m_population[ichromo].fitness = m_new_population[ichromo].fitness;
    }

    // check
    check();

    // sorting by fitness
    qsort(m_population, m_nchromo, sizeof(chromosome), compare);

    // updating best_chromosome/fitness
    m_best_chromosome.sequence = m_population[0].sequence;
    m_best_chromosome.fitness = m_population[0].fitness;
    m_best_fitness.push_back(m_best_chromosome.fitness);
    
    // calculating best-half <L> 
    for(int ichromo=0; ichromo<int(m_nchromo/2); ichromo++) L += m_population[ichromo].fitness;
    m_L.push_back( L/int(m_nchromo/2) );
}



void genetic_algorithm::clear() const
{
    delete []m_population;
    delete []m_new_population;
}



void genetic_algorithm::check() const
{
    int ncities = m_cities.size();

    for(int ichromo=0; ichromo<m_nchromo; ichromo++)
    {
        if( m_population[ichromo].sequence[0] != 1)
        {
            cout << "Error (genetic_algorithm::check): first city must always be equal to 1. Exit.";
            exit(1);
        }

        if( get_sum(m_population[ichromo].sequence) != int(ncities*(ncities+1)/2) )
        {
            cout << "Error (genetic_algorithm::check): the same city cannot be visited twice. Exit.";
            exit(1);
        }
    }
}



double genetic_algorithm::fitness(vector<double> sequence) const
{
    double distance = 0;
    vector<double> appo;

    for(int i=0; i<int(sequence.size()-1); i++)
    {
        appo = m_cities[sequence[i+1]-1] - m_cities[sequence[i]-1];
        distance += get_mod(appo);
    }
    
    appo = m_cities[sequence[sequence.size()-1]-1] - m_cities[sequence[0]-1];
    distance += get_mod(appo);

    return distance;
}



void genetic_algorithm::set_nchromo(int nchromo)
{
    delete []m_population;
    delete []m_new_population;
    m_nchromo = nchromo;
    m_population = new chromosome[nchromo];
    m_new_population = new chromosome[nchromo];
}



void genetic_algorithm::print_generation(int igen) const
{
    cout << endl;
    cout << igen << "th GENERATION" << endl;
    for(int ichromo=0; ichromo<m_nchromo; ichromo++)
    {
        cout << "- Chromosome: " << ichromo+1 << endl;
        
        cout << "\t sequence: ";
        for(auto& city:m_population[ichromo].sequence) cout << city << " - ";
        cout << "1" << endl;
        
        cout << "\t fitness:  ";
        cout << m_population[ichromo].fitness << endl << endl;
    }
}



void genetic_algorithm::print_best(int igen) const
{
    if(igen == 0) cout << endl << "-------------------------------------" << endl;

    cout << endl << igen << "th GENERATION" << endl;
    
    cout << "- Best chromosome: " << endl;
    
    cout << "\t sequence: ";
    for(auto& city:m_best_chromosome.sequence) cout << city << " - ";
    cout << "1" << endl;
    
    cout << "\t fitness:  ";
    cout << m_best_chromosome.fitness << endl;
    
    cout << endl << "-------------------------------------" << endl;
}



void genetic_algorithm::save_chromosome(string filename) const
{
    vector<double> coordinates;
	ofstream fout(filename);
	
	fout << "step\tcity_sequence\tcity_x\tcity_y\tfitness" << endl;
	for(int i=0; i<int(m_best_chromosome.sequence.size()); i++)
    {
        coordinates = m_cities[int(m_best_chromosome.sequence[i] - 1)];

        if(i == 0) fout << i << "\t" << m_best_chromosome.sequence[i] << "\t" << coordinates[0] << "\t" << coordinates[1] << "\t" << m_best_chromosome.fitness << endl;
        else fout << i << "\t" << m_best_chromosome.sequence[i] << "\t" << coordinates[0] << "\t" << coordinates[1] << "\t" << "" << endl;
    }

    coordinates = m_cities[0];
    fout << int(m_best_chromosome.sequence.size()) << "\t" << 1 << "\t" << coordinates[0] << "\t" << coordinates[1] << "\t" << "" << endl;   
	
	fout.close();
}



void genetic_algorithm::save_fitness(string filename) const
{
    vector<double> coordinates;
	ofstream fout(filename);
	
	fout << "gen\tfitness" << endl;
	for(int i=0; i<int(m_best_fitness.size()); i++)
        fout << i << "\t" << m_best_fitness[i]<< endl;

    fout.close();
}



void genetic_algorithm::save_cities(string filename) const
{
	ofstream fout(filename);
	
	fout << "city_num\tcity_x\tcity_y" << endl;
	for(int i=0; i<int(m_cities.size()); i++)
        fout << i+1 << "\t" << m_cities[i][0] << "\t" << m_cities[i][1] << endl;   
	
	fout.close();
}






// sorting function
int compare(const void *aa, const void *bb)
{
    chromosome *a, *b;
    a = (chromosome *) aa;
    b = (chromosome *) bb;

    if (a->fitness > b->fitness) return 1;
    else if (a->fitness < b->fitness) return -1;
    else return 0;
}