CC = g++
CFLAGS = -Wall -O3 --std=c++17
CPATH = /home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/

main.exe : main.o random.o mylib.o randomwalk.o
	$(CC) random.o mylib.o randomwalk.o main.o -o main.exe

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
	
random.o : $(CPATH)random.cpp $(CPATH)random.h
	$(CC) -c $(CPATH)random.cpp -o random.o $(CFLAGS)

mylib.o : $(CPATH)mylib.cpp $(CPATH)mylib.h
	$(CC) -c $(CPATH)mylib.cpp -o mylib.o $(CFLAGS)

randomwalk.o : $(CPATH)randomwalk.cpp $(CPATH)randomwalk.h
	$(CC) -c $(CPATH)randomwalk.cpp -o randomwalk.o $(CFLAGS)

clean :
	rm *.o *.exe
