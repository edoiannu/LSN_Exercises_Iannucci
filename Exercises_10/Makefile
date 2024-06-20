CFLAGS = -O3 --std=c++11
LPATH = /home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/

main.exe : random.o PT_fun.o mylib.o main.cpp
	mpicxx -o main.exe main.cpp random.o PT_fun.o mylib.o

execute : main.exe
	mpiexec -np 8 ./main.exe

random.o : $(LPATH)random.cpp $(LPATH)random.h
	$(CC) -c $(LPATH)random.cpp -o random.o $(CFLAGS)

PT_fun.o : PT_fun.cpp PT_fun.h
	$(CC) -c PT_fun.cpp -o PT_fun.o $(CFLAGS)

mylib.o : $(LPATH)mylib.cpp $(LPATH)mylib.h
	$(CC) -c $(LPATH)mylib.cpp -o mylib.o $(CFLAGS)

clean:
	rm *.exe *.o