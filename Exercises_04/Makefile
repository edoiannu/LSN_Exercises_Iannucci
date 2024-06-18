CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
CPATH = /home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/NSL_SIMULATOR/SOURCE/
RPATH = /home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/

all: main.exe
main.exe : system.o main.o particle.o random.o
	$(CC) system.o particle.o random.o main.o -o main.exe $(AFLAGS)
main.o : main.cpp $(CPATH)system.h $(CPATH)particle.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
system.o : $(CPATH)system.cpp $(CPATH)system.h $(RPATH)random.h $(CPATH)particle.h
	$(CC) -c $(CPATH)system.cpp -o system.o $(CFLAGS)
particle.o : $(CPATH)particle.cpp $(CPATH)particle.h
	$(CC) -c $(CPATH)particle.cpp -o particle.o $(CFLAGS)
random.o : $(RPATH)random.cpp $(RPATH)random.h
	$(CC) -c $(RPATH)random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe
remove:
	rm results/SOLID/*.* results/LIQUID/*.* results/GAS/*.* results/SOLID/CONFIG/*.* results/LIQUID/CONFIG/*.* results/GAS/CONFIG/*.*
