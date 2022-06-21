#!/bin/bash

# presentazione
echo "--------------------------------------------------------------"
echo ""

echo "Classic Lennard-Jones fluid"
echo "MD(NVE) simulation"
echo "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]"
echo "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T "
echo "The program uses Lennard-Jones units "
echo ""

echo "Simulations in three phases of matter: solid, liquid, gas."

# parametri generali
sed -i "1s/.*/0/" input.in             # simul
sed -i "17s/.*/MD/" input.in           # dir
sed -i "9s/.*/100/" input.in           # nblk
sed -i "10s/.*/10000/" input.in        # nstep
sed -i "2s/.*/0/" input.in             # verbose
sed -i "13s/.*/0.01/" input.in         # prec

for state in solid liquid gas
do
    # parametri stato fisico
    sed -i "16s/.*/$state/" input.in        # state
    if [[ $state == "solid" ]]
    then
        sed -i "4s/.*/0.8/" input.in        # T
        sed -i "6s/.*/1.1/" input.in        # rho
        sed -i "7s/.*/2.2/" input.in        # rcut
        sed -i "8s/.*/0.0005/" input.in     # delta
        echo "-- $state: T=0.8, rho=1.1, rcut=2.2, delta=0.0005"
    elif [[ $state == "liquid" ]]
    then
        sed -i "4s/.*/1.1/" input.in        # T
        sed -i "6s/.*/0.8/" input.in        # rho
        sed -i "7s/.*/2.5/" input.in        # rcut
        sed -i "8s/.*/0.0005/" input.in     # delta
        echo "-- $state: T=1.1, rho=0.8, rcut=2.5, delta=0.0005"
    else
        sed -i "4s/.*/1.2/" input.in        # T
        sed -i "6s/.*/0.05/" input.in       # rho
        sed -i "7s/.*/5.0/" input.in        # rcut
        sed -i "8s/.*/0.0005/" input.in     # delta
        echo "-- $state: T=1.2, rho=0.05, rcut=5.0, delta=0.0005"
    fi

    # simulazione
    echo -ne "   start... "
    ./MD_MC.exe
    echo -ne "   ...end\n"
done

echo ""
echo "--------------------------------------------------------------"