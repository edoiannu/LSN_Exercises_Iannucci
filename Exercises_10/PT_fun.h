#ifndef PT_FUN_H
#define PT_FUN_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "../libraries/random.h"

using namespace std;

const unsigned int N_PROV = 110;
const unsigned int N_IND = 100;  // n individual
const unsigned int N_GENE = N_PROV*N_IND;
const unsigned int SURV_CUT = N_IND / 10;

void substitute(int *v, int * vnew, int j);

// function that performs a pair permutation of two provs for the i-th element of the population
int * pair_perm(int *v, int i, int j, int k);

// function that performs a shift of n position of m contiguous provs for the i-th element of the population
int * shift(int *v, int i, int start, int m, int n);

// function that performs a permutation among m contiguous provs with other different m contiguous provs (m < N/2) for the i-th element of the population
int * contig_perm(int *v, int i, int m, int start1, int start2);

// inversion of the order in which they appear in the path of m provs (m <= M)
int * inversion(int *v, int i, int m, int start);

int * crossover(int * v, int i, int j, int cut);

int pbc(int i);

double distance(double x1, double y1, double x2, double y);

double loss_function(double x[], double y[], int * v, int i);

bool check(int *v);

double fitness(double x[], double y[], int * v);

void print_genome(int * v);

#endif