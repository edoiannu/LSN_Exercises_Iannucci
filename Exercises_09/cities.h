#ifndef _Cities_
#define _Cities_

#include <iostream>
#include <armadillo>
#include <cmath>
#include <iomanip>

#include "../libraries/random.h"
#include "../libraries/mylib.h"

using namespace std;
using namespace arma;

const unsigned int N_CITIES = 34;

class Cities {

private:
    Random _rnd;
    mat pos;    // vector containing the positions of the cities

public:
    Cities();
    void create_config(int );
    double distance(int , int );   // calculate the distance between two cities
    int pbc(int );
    double get_coord(int , int );

};

#endif