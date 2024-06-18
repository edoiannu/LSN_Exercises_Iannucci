#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "../../../libraries/random.h"
#include "../../../libraries/mylib.h"

using namespace std;

const int M = 1000000;    // number of throws
const int N = 100;      // number of blocks

const string ResultsDirectory = "../../results/";

int main(){

    Random rnd;     // definition of an object of type Random
    vector<double> I(M);
    vector<double> ave(N);
    vector<double> ave2(N);

    for (int i=0; i<M; i++)
        I[i] = (M_PI / 2) * cos(M_PI * rnd.Rannyu() / 2);

    ave = MeanBlock(I, M, N);

    for (int i=0; i<N; i++)         
        ave2[i] = ave[i] * ave[i];

    vector<double> prog_sum(N);
    vector<double> prog_sum2(N);

    ofstream WriteResults;
    WriteResults.open(ResultsDirectory+"02-1-1.out");

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