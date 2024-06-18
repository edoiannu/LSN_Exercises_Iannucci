#ifndef __Metropolis__
#define __Metropolis__

#include <iostream>
#include <fstream>
#include <armadillo>
#include <iomanip>
#include "../libraries/random.h"
#include "../libraries/mylib.h"

using namespace std;
using namespace arma;

class Metropolis {

private:
    int _ndim;  // dimensionality of the system
    int _dis_type;
    Random _rnd;    // object of type Random
    vec _coords;    // coordinates vector
    int _step_index = 0;
    int _blk_index = 0;
    int _acc_step = 0;
    vec _acceptance; // vector filled with the acceptance of each block
    int _nsteps;
    int _nblocks;
    int _blksteps;

    // properties
    string property;

public:
    void initialize(string inputf);
    void setcoord(int i, double x);
    bool step(double (*fun)(double , double , double));
    void increase_step();
    double get_coord(int i);
    // void print_acceptance(string outputd);
    int get_nsteps();
    int get_nblocks();
    int get_blksteps();
    void print_results(vec ave, string outputf);
    int get_distype();
};

using namespace std;

#endif