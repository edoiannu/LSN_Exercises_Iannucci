#include "../metropolis.h"

const string ResultsDirectory = "results/";

double temp(int step){
    double temp = static_cast<double>(3./(step+1));
    // cout << step << " " << temp << endl;
    return temp;
}

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

    double appo = 0.;
    vec ave(N);

    int M_SA = metro.get_nsasteps();

    for (int i=0; i<M_SA; i++){

        metro.settemp(temp(i));

        for (int j=0; j<N; j++){

            appo = 0.;

            for(int k=j*L; k<L*(j+1); k++){
                metro.step(wavef_prob);
                appo += energy(metro.get_vcoord(), metro.get_vparam());
            }

            ave[j] = appo / L;
        }

        appo = 0;

        if (i==0){
            for (int j=0; j<N; j++)
                appo += ave[j];

            appo /= N;

            metro.setenermin(appo);
        } 

        metro.print_results(ave);
        metro.sa_step(wavef_prob, energy);

        metro.restart(); 
    }

    metro.close_stream();

    return 0;
}