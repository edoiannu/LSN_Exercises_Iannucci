#ifndef _MyGA_
#define _MyGA_

#include "cities.h"

const unsigned int POP_DIM = 1000;
const unsigned int SURV_CUT = POP_DIM/10;

class MyGA {

private:
    Random _rnd;
    Cities _cit;
    mat popul;

public:
    MyGA();
    void set_config(string );
    int get_gene(int , int );
    vec mutation(vec , double );
    vec pair_perm(vec );
    vec shift(vec );
    vec contig_perm(vec );
    vec inversion(vec );
    mat crossover(int , int , double );
    void substitute(vec , int );
    bool check();
    double loss_function(int );
    int selection();
    double fitness();
    vec get_chromosom(int );
    void print_chromosom(int );
    void print_popul();
    void print_coords(int , string );
};

#endif