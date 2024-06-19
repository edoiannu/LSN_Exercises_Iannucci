CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
LPATH = ../../libraries/

main.exe : myGA.o random.o cities.o main.o mylib.o
	$(CC) myGA.o random.o cities.o main.o mylib.o -o main.exe $(AFLAGS)
main.o : main.cpp myGA.h $(LPATH)random.h cities.h $(LPATH)mylib.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : $(LPATH)random.cpp $(LPATH)random.h
	$(CC) -c $(LPATH)random.cpp -o random.o $(CFLAGS)
cities.o : cities.cpp cities.h
	$(CC) -c cities.cpp -o cities.o $(CFLAGS)
myGA.o : myGA.cpp myGA.h
	$(CC) -c myGA.cpp -o myGA.o $(CFLAGS)
mylib.o : $(LPATH)mylib.cpp $(LPATH)mylib.h
	$(CC) -c $(LPATH)mylib.cpp -o mylib.o $(CFLAGS)
clean :
	rm *.o *.exe
remove :
	rm results/circle/* results/square/*
