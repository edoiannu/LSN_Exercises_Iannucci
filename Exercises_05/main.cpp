#include <iostream>
#include <vector>
#include <cstdlib> // Per la funzione exit
#include <string>  // Per std::stoi
#include "metropolis.h"

using namespace std;

const string ResultsDirectory = "results/";

const int n = 2;
const int l = 1;
const int m = 0;

double wavef_prob(double x, double y, double z){

    if (n == 1 and l == 0 and m == 0){
        return (1. / M_PI) * exp(-2*sqrt(x * x + y * y + z * z));
    }
    else if (n == 2 and l == 1 and m == 0){
        return (1./32.) * (1. / M_PI) * exp(-sqrt(x * x + y * y + z * z)) * z * z;
    }
    else { return 0; }
}

double r(double x, double y, double z){

    return sqrt(x * x + y * y + z * z);
}

int main(){

    Metropolis metro;
    metro.initialize("metropolis.in");

    const int M = metro.get_nsteps();
    const int N = metro.get_nblocks();
    const int L = metro.get_blksteps();

    string resultsf;

    if (n == 1 and l == 0 and m == 0 and metro.get_distype() == 0) { resultsf = "1s-uniform.out"; }
    else if (n == 1 and l == 0 and m == 0 and metro.get_distype() == 1) { resultsf = "1s-gaussian.out"; }
    else if (n == 2 and l == 1 and m == 0 and metro.get_distype() == 0) { resultsf = "2p-uniform.out"; }
    else if (n == 2 and l == 1 and m == 0 and metro.get_distype() == 1) { resultsf = "2p-gaussian.out"; }

    vec ave(N);

    double appo;
    bool appo_;

    ofstream print_steps;
    print_steps.open(ResultsDirectory+"steps-"+resultsf);

    for (int i=0; i<N; i++){

        appo = 0.;
        
        for(int j=i*L; j<L*(i+1); j++){

            appo_ = metro.step(wavef_prob);
            appo += r(metro.get_coord(0), metro.get_coord(1), metro.get_coord(2));

            if(appo_){
                print_steps << metro.get_coord(0) << " " << metro.get_coord(1) << " " << metro.get_coord(2) << endl;
            }

            // cout << j << endl;
        }

        ave[i] = appo / L;
    }

    // metro.print_results(ave, ResultsDirectory+resultsf);
    print_steps.close();

    return 0;
}