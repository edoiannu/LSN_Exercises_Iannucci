CC = g++
CFLAGS = -Wall -O3 --std=c++17
CPATH = /home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
	
random.o : $(CPATH)random.cpp $(CPATH)random.h
	$(CC) -c $(CPATH)random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o *.exe results/*.out
