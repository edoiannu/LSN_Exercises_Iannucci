#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "/home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/random.h"
#include "/home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/mylib.h"

using namespace std;

const int M = 10000;    // number of throws
const int N = 100;      // number of blocks

const string ResultsDirectory = "/home/edoiannu/Documenti/Lab_simulazione_numerica/LSN_Exercises_Iannucci/Exercises_02/results/";

double f (double x) { return 1.-sqrt(1.-x); };

int main(){

    Random rnd;     // definition of an object of type Random
    vector<double> r(M);
    vector<double> I(M);
    vector<double> ave(N);
    vector<double> ave2(N);

    for (int i=0; i<M; i++)
        r[i] = rnd.Distributionf(f);

    for (int i=0; i<M; i++){
        I[i] = (M_PI / 2) * cos(M_PI * r[i] / 2) / (2 * (1 - r[i]));
        // cout << rnd.Distributionf(f) << endl;
        }

    ave = MeanBlock(I, M, N);

    for (int i=0; i<N; i++)         
        ave2[i] = ave[i] * ave[i];

    vector<double> prog_sum(N);
    vector<double> prog_sum2(N);

    ofstream WriteResults;
    WriteResults.open(ResultsDirectory+"02-1-2.out");

    WriteResults << "# STEPS " << M << endl;
    WriteResults << "# BLOCKS " << N << endl;

    for (int i=0; i<N; i++){

        for(int j=0; j<=i; j++){

            prog_sum[i] += ave[j];
            prog_sum2[i] += ave2[j];

        }

        prog_sum[i] /= (i+1);
        prog_sum2[i] /= (i+1);

        WriteResults << prog_sum[i] << setw(12) << var(prog_sum[i], prog_sum2[i],i) << endl;

    }

    WriteResults.close();

    return 0;
}