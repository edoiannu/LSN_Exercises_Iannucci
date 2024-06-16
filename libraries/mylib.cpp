// library of useful functions

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "mylib.h"

using namespace std;

double var (double A, double A2, int n){
    // this function returns the statistical uncertainty (i.e. the standard deviation of the mean) of an event given its mean A the mean of the squares A2 and the number of samples n
   if ( n == 0 ) { return 0; }
   return sqrt((A2-A*A)/n);}

vector<double> MeanBlock (vector<double> v, int M, int N){

    double sum;             // helper variable
    vector<double> ave(N);  // array of the N means of each block

    int L = M / N;  // number of throws in each block

    for(int i=0; i<N; i++){    // for loop that iterates over the blocks

            sum = 0;

            for(int j=i*L; j<(i+1)*L; j++)  // for loop that iterates over the throws of the (j+1)-th block
                    sum += v[j];

            ave[i] = sum / L;      // <A_i>
        }

    return ave;
}

double Theta (double x, double y){

    return atan(y/x);
}

// Function to find the minimum value between x_min and x_max
double min (double x_min, double x_max){
    if (x_min < x_max) { return x_min; }
    else if (x_min > x_max) { return x_max; }
    else { return 0; }
}