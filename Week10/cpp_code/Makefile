MC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o utils.o genetic_algorithm.o stats.o
	$(MC) -no-pie random.o main.o stats.o utils.o genetic_algorithm.o -o main.exe
main.o : main.cpp
	$(MC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(MC) -c random.cpp -o random.o $(CFLAGS)
utils.o : utils.cpp utils.h random.h
	$(MC) -c utils.cpp -o utils.o $(CFLAGS)
stats.o : stats.cpp stats.h
	$(MC) -c stats.cpp -o stats.o $(CFLAGS)
genetic_algorithm.o : genetic_algorithm.cpp genetic_algorithm.h random.h utils.h
	$(MC) -c genetic_algorithm.cpp -o genetic_algorithm.o $(CFLAGS)
clean :
	rm *.o seed.out
