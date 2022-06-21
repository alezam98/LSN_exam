CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o stats.o utils.o
	$(CC) -no-pie random.o main.o stats.o utils.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
stats.o : stats.cpp stats.h
	$(CC) -c stats.cpp -o stats.o $(CFLAGS)
utils.o : utils.cpp utils.h random.h
	$(CC) -c utils.cpp -o utils.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
