#!/bin/bash

# Presentazione
echo "Classic 1D Ising model"
echo "Monte Carlo simulations"
echo "Nearest neighbour interaction"
echo "Boltzmann weight exp(- beta * H ), beta = 1/T"
echo "The program uses k_B=1 and mu_B=1 units"
echo ""


# verbose = 0
sed -i "1s/.*/0/" input.dat


# ciclo di simulazioni: algoritmo di Metropolis e di Gibbs
for algo in Metro Gibbs
do

    # fisso l'algoritmo
    if [[ $algo == Metro ]]
    then
        echo ""
        echo "**********************************"
        echo "Metropolis algorithm"
        sed -i "8s/.*/1/" input.dat
    else
        echo ""
        echo "**********************************"
        echo "Gibbs sampling"
        sed -i "8s/.*/0/" input.dat
    fi

    # ciclo di simulazioni con h variabile, h = (0.0, 0.02)
    for h in 0.0 0.02
    do

        # fisso h
        sed -i "7s/.*/$h/" input.dat

        # Singola simulazione a temperatura fissata, T=0.5
        echo "-- Single simulation: T=0.5 - h=$h"

        sed -i '4s/.*/0.5/' input.dat
        ./Monte_Carlo_ISING_1D.exe

        if [[ $h == 0.0 ]]
        then
            rm output.mag.dat
            mv output* ./$algo/Single_sim
        else
            rm output.ene.dat output.heat.dat output.chi.dat
            mv output* ./$algo/Single_sim
        fi

        echo "----------------------------------"
        echo ""

        # Serie di simulazioni a temperatura variabile, T=[0.5, 2.0]
        echo "-- Multi simulation: T=[0.5, 2.0] - h=$h"

        for T in 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0
        do

            # fisso T
            echo "  -T=$T"

            sed -i "4s/.*/$T/" input.dat
            ./Monte_Carlo_ISING_1D.exe


            if [[ $T == 0.5 ]]
            then

                if [[ $h == 0.0 ]]
                then
                    rm output.mag.dat

                    tail -1 output.ene.dat > multi.ene.dat
                    sed -i "s/$/\t$T/" multi.ene.dat >> multi.ene.dat

                    tail -1 output.heat.dat > multi.heat.dat
                    sed -i "s/$/\t$T/" multi.heat.dat >> multi.heat.dat

                    tail -1 output.chi.dat > multi.chi.dat
                    sed -i "s/$/\t$T/" multi.chi.dat >> multi.chi.dat
                else
                    rm output.ene.dat output.heat.dat output.chi.dat

                    tail -1 output.mag.dat > multi.mag.dat
                    sed -i "s/$/\t$T/" multi.mag.dat >> multi.mag.dat
                fi
            
            else
                if [[ $h == 0.0 ]]
                then
                    rm output.mag.dat

                    tail -1 output.ene.dat > multi.ene0.dat
                    sed -i "s/$/\t$T/" multi.ene0.dat >> multi.ene0.dat
                    cat multi.ene0.dat >> multi.ene.dat

                    tail -1 output.heat.dat > multi.heat0.dat
                    sed -i "s/$/\t$T/" multi.heat0.dat >> multi.heat0.dat
                    cat multi.heat0.dat >> multi.heat.dat

                    tail -1 output.chi.dat > multi.chi0.dat
                    sed -i "s/$/\t$T/" multi.chi0.dat >> multi.chi0.dat
                    cat multi.chi0.dat >> multi.chi.dat
                else
                    rm output.ene.dat output.heat.dat output.chi.dat

                    tail -1 output.mag.dat > multi.mag0.dat
                    sed -i "s/$/\t$T/" multi.mag0.dat >> multi.mag0.dat
                    cat multi.mag0.dat >> multi.mag.dat
                fi

            fi

        done

        echo "----------------------------------"
        echo ""

        rm output* 
        mv multi* ../data/$algo/Multi_sim
        cd ../data/$algo/Multi_sim
        rm multi.ene0.dat multi.heat0.dat multi.chi0.dat multi.mag0.dat
        cd ..
        cd ..

    done

    echo "**********************************"
    echo ""

done
