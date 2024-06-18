CC = g++
CFLAGS = -Wall -O3 --std=c++17
CPATH = ../../../libraries/

main.exe : main.o random.o mylib.o
	$(CC) random.o mylib.o main.o -o main.exe

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
	
random.o : $(CPATH)random.cpp $(CPATH)random.h
	$(CC) -c $(CPATH)random.cpp -o random.o $(CFLAGS)

mylib.o : $(CPATH)mylib.cpp $(CPATH)mylib.h
	$(CC) -c $(CPATH)mylib.cpp -o mylib.o $(CFLAGS)

clean :
	rm *.o *.exe
