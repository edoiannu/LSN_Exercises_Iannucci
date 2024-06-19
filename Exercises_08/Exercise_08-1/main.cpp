#include <iostream>
#include "../metropolis.h"

const string ResultsDirectory = "results/";

using namespace std;



double wavef(vec coord, vec param){

    double mu = param[0];
    double sigma = param[1];

    return exp( -0.5*(coord[0]-mu)*(coord[0]-mu) / (sigma*sigma) ) + exp( -0.5*(coord[0]+mu)*(coord[0]+mu) / (sigma*sigma) );
}

double wavef_prob(vec coord, vec param){

    double mu = param[0];
    double sigma = param[1];

    double appo = exp( -0.5*(coord[0]-mu)*(coord[0]-mu) / (sigma*sigma) ) + exp( -0.5*(coord[0]+mu)*(coord[0]+mu) / (sigma*sigma) );

    return appo*appo;
}

double energy(vec coord, vec param){

    double x = coord[0];
    double mu = param[0];
    double sigma = param[1];

    double wavef = exp( -0.5*(x-mu)*(x-mu) / (sigma*sigma) ) + exp( -0.5*(x+mu)*(x+mu) / (sigma*sigma) );

    double kenergy = (-0.5)*(exp( -0.5*(x-mu)*(x-mu) / (sigma*sigma) ) * ((x-mu)*(x-mu)/pow(sigma,4) - 1/(sigma*sigma) ) + exp( -0.5*(x+mu)*(x+mu) / (sigma*sigma) ) * ((x+mu)*(x+mu)/pow(sigma,4) - 1/(sigma*sigma) )) / wavef;

    double penergy = pow(x,4) - (5./2.)*x*x;

    return kenergy + penergy;
}

int main(){

    Metropolis metro;
    metro.initialize("metropolis.in", ResultsDirectory);

    const int M = metro.get_nsteps();
    const int N = metro.get_nblocks();
    const int L = metro.get_blksteps();

    vec ave(N);

    double appo;

    for (int i=0; i<N; i++){

        appo = 0.;

        for(int j=i*L; j<L*(i+1); j++){

            metro.step(wavef_prob);
            appo += energy(metro.get_vcoord(), metro.get_vparam());
        }

        ave[i] = appo / L;
    }

    metro.print_results(ave);

    return 0;
}