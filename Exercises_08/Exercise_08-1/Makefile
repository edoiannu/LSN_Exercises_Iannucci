CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
LPATH = ../../libraries/
MPATH = ../

main.exe : metropolis.o main.o random.o mylib.o
	$(CC) metropolis.o random.o mylib.o main.o -o main.exe $(AFLAGS)
main.o : main.cpp $(MPATH)metropolis.h $(LPATH)random.h $(LPATH)mylib.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
metropolis.o : $(MPATH)metropolis.cpp $(MPATH)metropolis.h $(LPATH)random.h $(LPATH)mylib.h
	$(CC) -c $(MPATH)metropolis.cpp -o metropolis.o $(CFLAGS)
random.o : $(LPATH)random.cpp $(LPATH)random.h
	$(CC) -c $(LPATH)random.cpp -o random.o $(CFLAGS)
mylib.o : $(LPATH)mylib.cpp $(LPATH)mylib.h
	$(CC) -c $(LPATH)mylib.cpp -o mylib.o $(CFLAGS)
clean :
	rm *.o main.exe
remove :
	rm results/*.*
