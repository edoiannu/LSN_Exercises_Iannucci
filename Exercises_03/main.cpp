#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <vector>
#include "/home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/random.h"

using namespace std;

const int M = 10000;        // number of asset prices
const int N = 100;          // number of blocks
const int L = M / N;        // number of asset prices within a block

const double r = 0.1;       // risk free interest rate
const double K = 100;       // strike price
const double sigma = 0.25;  // volatility
const double S0 = 100;      // asset price at t = 0
const double T = 1.;        // delivery time
const double dt = T/100;

const string path = "/home/edoiannu/Documenti/Lab_simulazione_numerica/Exercises_03/results/";

double var (double A, double A2, int n){
   return sqrt((A2-A*A)/n);}

double price (double mu, double sigma, double t, double x, double S_0){
    return S_0 * exp((mu - 0.5*sigma*sigma)*t + sigma*x*sqrt(t));
}

double max (double x_min, double x_max){
    if (x_min < x_max) { return x_max; }
    else if (x_min > x_max) { return x_min; }
    else if (x_max == x_min) { return 0; }
}

int main(){

    Random rnd;
    vector<double> S(M);
    double sum;
    vector<double> ave(N);
    vector<double> ave2(N);
    vector<double> prog_sum(N);
    vector<double> prog_sum2(N);

    ofstream WriteResults;

    string options[2] = {"call", "put"};
    string sample[2] = {"direct", "discretized"};

    for (int a=0; a<2; a++){    // ciclo sui sample
    for (int b=0; b<2; b++){    // ciclo sulle options

    S.clear();

    for (int i=0; i<N; i++)
        ave[i] = ave2[i] = prog_sum[i] = prog_sum2[i] = 0;

    for (int i=0; i<M; i++){
        if (sample[a] == "direct") { S[i] = price(r, sigma, T, rnd.Gauss(0, T), S0);
        // cout << S[i] << endl;
        }
        else if (sample[a] == "discretized") {
            for (int t=0; t<100; t++){
                if (t==0) { S[i]=S0; }
                S[i] = price(r, sigma, dt, rnd.Gauss(0,1), S[i]);
            }
        }
    }
    
    for (int i=0; i<N; i++){

        sum = 0;

        for (int j=i*L; j<(i+1)*L; j++){
            if (options[b] == "call") { sum += exp(-r*T) * max(0,S[j]-K); }
            else if (options[b] == "put") { sum += exp(-r*T) * max(0,K-S[j]); }
        }

        ave[i] = sum / L;          
        ave2[i] = ave[i]*ave[i];
    }

    WriteResults.open(path+options[b]+"_"+sample[a]+".out");

    for (int i=1; i<N; i++){

        for(int j=0; j<=i; j++){

            prog_sum[i] += ave[j];
            prog_sum2[i] += ave2[j];

        }

        prog_sum[i] /= (i+1);
        prog_sum2[i] /= (i+1);

        WriteResults << prog_sum[i] << setw(12) << var(prog_sum[i], prog_sum2[i],i) << endl;

    }

    WriteResults.close();
    
    }
    }

    return 0;
}