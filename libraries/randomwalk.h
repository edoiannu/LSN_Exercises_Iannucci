#ifndef __RandomWalk__
#define __RandomWalk__

#include <cmath>
#include "/home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/random.h"

class RandomWalk {

private:
double x,y,z;
int steps;

protected:

public:

RandomWalk();

void SetCoords(double , double , double );

double GetX();

double GetY();

double GetZ();

void Evolvex_lattice(int , double );

void Evolvey_lattice(int , double );

void Evolvez_lattice(int , double );

void Evolve_polar(double , double , double );

void Evolve_coords(double , double , double );

double Getr2();

};

#endif