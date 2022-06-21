#include "main.h"

int main()
{
    // initialization
    Input();
    if(verbose) cout << "Initializing the population...";
    algorithm->initialize(rnd);
    if(verbose) cout << " end." << endl;

    // simulation
    if(verbose) cout << "Start of the simulation.";
    for(int igen=1; igen<=ngen; igen++)
    {
        update_index = 0;
        do
        {
            algorithm->crossover(update_index, rnd);
            algorithm->mutation(update_index, rnd);
            update_index += 2;
        } while(update_index < nchromo);
        
        algorithm->update();
        if(verbose) algorithm->print_best(igen);
    }
    algorithm->clear();

    mean_L = algorithm->get_L();
    gen.block_mean(mean_L);
    if(verbose) cout << "End of the simulation." << endl;
    
    // saving data
    if(verbose) cout << "Saving data...";
    gen.save("../data/"+dir+"/"+shape+"/mean_L.dat");
    algorithm->save_fitness("../data/"+dir+"/"+shape+"/fitness.dat");
    algorithm->save_chromosome("../data/"+dir+"/"+shape+"/chromosome.dat");
    algorithm->save_cities("../data/"+dir+"/"+shape+"/cities.dat");
    if(verbose) cout << " saved." << endl << endl;
    
    rnd.SaveSeed();
    return 0;
}



void Input(void)
{
    ifstream ReadInput, Primes, Seed;

    //Read input informations
    ReadInput.open("input.in");

    ReadInput >> verbose;

    ReadInput >> restart;

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

    ReadInput >> pc;

    ReadInput >> pm_1;

    ReadInput >> pm_2;

    ReadInput >> pm_3;

    ReadInput >> pm_4;

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

    // initialize classes
    algorithm->set_nchromo(nchromo);
    algorithm->set_cities(cities);
    algorithm->set_pc(pc);
    algorithm->set_pms(pms);


    if(verbose)
    {
        cout << endl;
        cout << "Genetic algorithm, the traveling salesman problem." << endl;
        cout << "The chromosome population evolves through crossover and mutations." << endl;
        cout << "Each chromosome is defined by a sequence (i.e. visited city) and a fitness (i.e. goodness of the sequence)." << endl;
        cout << "Every sequence starts from the first city." << endl;
        cout << "Each city is randomly placed ";
        if(circle) cout << "on a circumference of radius r=1." << endl;
        else if(square) cout << "inside a square of length l=1." << endl;

        cout << endl;
        cout << "Parameters" << endl;
        cout << "Number of cities: " << ncities << endl;
        cout << "Number of chromosomes per generation: " << nchromo << endl;
        cout << "Number of generations: " << ngen << endl;
        cout << "Crossover probability: " << pc << endl;
        cout << "Mutation probability: " << pm_1 << " (pair swap)" << endl;
        cout << "                      " << pm_2 << " (shift)" << endl;
        cout << "                      " << pm_3 << " (contiguos swap)" << endl;
        cout << "                      " << pm_4 << " (inversion)" << endl << endl;
    }
    
    ReadInput.close();
}
