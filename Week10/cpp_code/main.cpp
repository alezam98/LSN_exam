#include "main.h"

int main(int argc, char* argv[])
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(size%2 != 0)
    {
        if(rank == 0) cout << "Error: number of processes must be even. Exit." << endl;
        exit(-1);
    }

    // initialization
    Input(rank);
    if(verbose && rank == 0) cout << "- Initializing the populations...";
    algorithm->initialize(rnd);
    if(verbose && rank == 0) cout << " end." << endl;

    // simulation
    if(verbose && rank == 0) cout << "- Start of the simulation." << endl;
    for(int igen=1; igen<=ngen; igen++)
    {
        if(rank == 0 && igen%int(ngen/10) == 0) cout << "Generation: " << igen << endl;
        update_index = 0;
        do
        {
            algorithm->crossover(update_index, rnd);
            algorithm->mutation(update_index, rnd);
            update_index += 2;
        } while(update_index < nchromo);
        
        algorithm->update();
        if(nmigr%igen == 0 && size > 1) Migration(size, rank);
    }
    algorithm->clear();

    mean_L = algorithm->get_L();
    gen.block_mean(mean_L);
    if(verbose && rank == 0) cout << "- End of the simulation." << endl;
    
    // saving data
    if(verbose && rank == 0) cout << "- Saving data...";
    Save(rank);
    if(verbose && rank == 0) cout << " saved." << endl << endl;
    
    rnd.SaveSeed();
    MPI_Finalize();
    return 0;
}



void Input(int rank)
{
    ifstream ReadInput, Primes, Seed;

    //Read input informations
    ReadInput.open("input.in");

    ReadInput >> verbose;

    ReadInput >> restart;

    //Read seed for random numbers
    int p1, p2, irank = -1;
    Primes.open("Primes");
    do
    {
        Primes >> p1 >> p2;
        irank ++;
    } while(irank != rank);
    Primes.close();

    if(restart) Seed.open("seed.out");
    else Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    Seed.close();

    ReadInput >> nchromo;

    if(nchromo%2 == 1)
    {
        cout << "The number of chromosomes must be even. Substituting " << nchromo << " with " << nchromo+1 << "." << endl;
        cout << "Do you want to proceed? (y/n) ";
        
        do{
            cin >> check;
            if(check == "y") nchromo ++;
            else if(check == "n") exit(0);
        } while(check != "y" && check != "n");
    }

    ReadInput >> ncities;

    ReadInput >> ngen;

    ReadInput >> nmigr;

    ReadInput >> exchange_fract;

    ReadInput >> pc;

    ReadInput >> pm_1;

    ReadInput >> pm_2;

    ReadInput >> pm_3;

    ReadInput >> pm_4;

    ReadInput >> city_type;

    if(city_type != "geometrical" && city_type != "capitals")
    {
        cout << "Error (Input): 'city_type' variable must be 'geometrical' or 'capitals'. Exit." << endl;
        exit(0);
    }

    ReadInput >> circle;

    ReadInput >> square;

    if(circle*square == 1)
    {
        cout << "Error (Input): city positions must be defined or on a circle or inside a square. Exit." << endl;
        exit(0);
    }

    ReadInput >> dir;

    // initialize vectors
    pms = {pm_1, pm_2, pm_3, pm_4};
    
    if(city_type == "geometrical")
    {
        if(circle)
        {
            for(int icity=0; icity<ncities; icity++) cities.push_back( rnd.bidim_versor() );
            shape = "circle";
        }
        else if(square)
        {
            for(int icity=0; icity<ncities; icity++) cities.push_back( rnd.Rannyu(2) );
            shape = "square";
        }
        
        // broadcasting city positions
        for(int icity=0; icity<ncities; icity++)                                    
            for(int ipos=0; ipos<2; ipos++) MPI_Bcast(&cities[icity][ipos], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    else if(city_type == "capitals")
    {
        ifstream ReadCity;
        ReadCity.open("capitals.dat");

        double lat, lon;
        while(!ReadCity.eof())
        {
            ReadCity >> lat >> lon;
            cities.push_back({lat, lon});
        }
        cities.erase(cities.end());
        ReadCity.close();

        ncities = cities.size();
    }


    // initialize classes
    algorithm->set_nchromo(nchromo);
    algorithm->set_cities(cities);
    algorithm->set_pc(pc);
    algorithm->set_pms(pms);


    if(verbose && rank == 0)
    {
        cout << endl;
        cout << "Genetic algorithm, the traveling salesman problem." << endl;
        cout << "The chromosome population evolves through crossover and mutations." << endl;
        cout << "Each chromosome is defined by a sequence (i.e. visited city) and a fitness (i.e. goodness of the sequence)." << endl;
        cout << "Every sequence starts from the first city." << endl;
        if(city_type == "geometrical")
        {
            cout << "Each city is randomly placed ";
            if(circle) cout << "on a circumference of radius r=1." << endl;
            else if(square) cout << "inside a square of length l=1." << endl;
        }
        else if(city_type == "capitals")
        {
            cout << "Each city has the coordinates (latitude and longitude) of an american capital." << endl;
        }

        cout << endl;
        cout << "The algorithm uses MPI library to search fot the optimal path." << endl;
        cout << "The different processes (Continents) communicate exchanging a fraction of their best sequences (Migration) every fixed number of steps." << endl;

        cout << endl;
        cout << "Parameters" << endl;
        cout << "Number of cities: " << ncities << endl;
        cout << "Number of chromosomes per generation: " << nchromo << endl;
        cout << "Number of generations: " << ngen << endl;
        cout << "Number of generations between migrations: " << nmigr << endl;
        cout << "Fraction of exchanged chromosomes: " << exchange_fract << endl;
        cout << "Crossover probability: " << pc << endl;
        cout << "Mutation probability: " << pm_1 << " (pair swap)" << endl;
        cout << "                      " << pm_2 << " (shift)" << endl;
        cout << "                      " << pm_3 << " (contiguos swap)" << endl;
        cout << "                      " << pm_4 << " (inversion)" << endl << endl;
    }
    
    ReadInput.close();
}


void Migration(int size, int rank)
{
    chromosome* population = algorithm->get_population();
    vector<double> exchange_sequence;

    for(int irank=0; irank<size; irank++) exchange_sequence.push_back(irank);
    random_shuffle(exchange_sequence.begin(), exchange_sequence.end());

    for(int index=0; index<size; index+=2)
    {
        for(int ichromo=0; ichromo<int(nchromo*exchange_fract); ichromo++)
            Send_Receive(population[ichromo], exchange_sequence[index], exchange_sequence[index+1], rank);
    }
    algorithm->set_population(population, rank);
}



void Send_Receive(chromosome& chromo, int rank1, int rank2, int rank)
{
    MPI_Status stat1, stat2;
    int itag1 = 1, itag2 = 2;
    
    // exchanging sequences
    if(rank == rank1)
    {
        MPI_Send(&chromo.sequence[0], int(chromo.sequence.size()), MPI_DOUBLE, rank2, itag1, MPI_COMM_WORLD);
        MPI_Recv(&chromo.sequence[0], int(chromo.sequence.size()), MPI_DOUBLE, rank2, itag2, MPI_COMM_WORLD, &stat2);
    }
    else if(rank == rank2)
    {
        MPI_Send(&chromo.sequence[0], int(chromo.sequence.size()), MPI_DOUBLE, rank1, itag2, MPI_COMM_WORLD);
        MPI_Recv(&chromo.sequence[0], int(chromo.sequence.size()), MPI_DOUBLE, rank1, itag1, MPI_COMM_WORLD, &stat1);
    }

    // exchanging fitnesses
    if(rank == rank1)
    {       
        MPI_Send(&chromo.fitness, 1, MPI_DOUBLE, rank2, itag1, MPI_COMM_WORLD);
        MPI_Recv(&chromo.fitness, 1, MPI_DOUBLE, rank2, itag2, MPI_COMM_WORLD, &stat2);
    }
    else if(rank == rank2)
    {
        MPI_Send(&chromo.fitness, 1, MPI_DOUBLE, rank1, itag2, MPI_COMM_WORLD);
        MPI_Recv(&chromo.fitness, 1, MPI_DOUBLE, rank1, itag1, MPI_COMM_WORLD, &stat1);
    }
}


void Save(int rank)
{
    string rank_num;
    if(rank == 0) rank_num = "r0";
    else if(rank == 1) rank_num = "r1";
    else if(rank == 2) rank_num = "r2";
    else if(rank == 3) rank_num = "r3";

    if(city_type == "geometrical")
    {
        gen.save("./"+dir+"/"+shape+"/"+rank_num+"_mean_L.dat");
        algorithm->save_fitness("./"+dir+"/"+shape+"/"+rank_num+"_fitness.dat");
        algorithm->save_chromosome("./"+dir+"/"+shape+"/"+rank_num+"_chromosome.dat");
        algorithm->save_cities("./"+dir+"/"+shape+"/"+rank_num+"_cities.dat");
    }
    else if(city_type == "capitals")
    {
        gen.save("./"+dir+"/"+city_type+"/"+rank_num+"_mean_L.dat");
        algorithm->save_fitness("./"+dir+"/"+city_type+"/"+rank_num+"_fitness.dat");
        algorithm->save_chromosome("./"+dir+"/"+city_type+"/"+rank_num+"_chromosome.dat");
        algorithm->save_cities("./"+dir+"/"+city_type+"/"+rank_num+"_cities.dat");
    }
}
