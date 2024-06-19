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
    int _nsteps_sa;
    int _nblocks;
    int _blksteps;
    double _temp;
    double _beta;
    double _enermin;
    bool _sa;
    vec _param;    // simulated annealing vector parameters
    int _acc_step_sa = 0;
    int _step_index_sa = 0;
    ofstream WriteResults;
    ofstream WriteParams;

    // properties
    string property;

public:
    void initialize(string inputf, string outputd);
    void setcoord(int i, double x);
    void setparam(int i, double x);
    void settemp(double x);
    void setenermin(double x);
    void restart();
    void step(double (*fun)(vec ,vec ));
    void increase_step();
    double get_coord(int i);
    vec get_vcoord();
    vec get_vparam();
    // void print_acceptance(string outputd);
    int get_nsteps();
    int get_nblocks();
    int get_blksteps();
    int get_nsasteps();
    void print_results(vec ave);
    int get_distype();
    double get_temp();
    void sa_step(double (*fun)(vec ,vec ), double (*fun2)(vec ,vec ));
    void evolve_param(int i, double delta);
    void close_stream();
};

#endif