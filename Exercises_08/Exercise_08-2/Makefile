CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
LPATH = /home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/

main.exe : metropolis.o main.o random.o mylib.o
	$(CC) metropolis.o random.o mylib.o main.o -o main.exe $(AFLAGS)
main.o : main.cpp ../metropolis.h $(LPATH)random.h $(LPATH)mylib.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
metropolis.o : ../metropolis.cpp ../metropolis.h $(LPATH)random.h $(LPATH)mylib.h
	$(CC) -c ../metropolis.cpp -o metropolis.o $(CFLAGS)
random.o : $(LPATH)random.cpp $(LPATH)random.h
	$(CC) -c $(LPATH)random.cpp -o random.o $(CFLAGS)
mylib.o : $(LPATH)mylib.cpp $(LPATH)mylib.h
	$(CC) -c $(LPATH)mylib.cpp -o mylib.o $(CFLAGS)
clean :
	rm *.o main.exe
remove :
	rm results/*.*
