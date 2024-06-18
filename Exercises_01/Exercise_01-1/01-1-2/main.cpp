#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <vector>
#include "../../../libraries/random.h"
#include "../../../libraries/mylib.h"

using namespace std;

const int M = 100000;   // number of total throws
const int N = 100;      // number of blocks

// const string ResultsDirectory = "/home/edoiannu/Documenti/Lab_simulazione_numerica/Exercises_01/Exercise_01-1/01-1-2/";

int main(){

    Random rnd;     // definition of an object of type Random
    vector<double> r(M);    // array containing the M numbers randomly generated

    for (int i=0; i<M; i++) // for loop that loads an array with M numebers randmoly generated
        r[i] = rnd.Rannyu();

    vector<double> sigma2(M);

    for (int i=0; i<M; i++){
        sigma2[i] = (r[i] - 0.5) * (r[i] - 0.5);
    }

    vector<double> ave(N);  // <A_i>
    vector<double> ave2(N); // <A_i>^2

    ave = MeanBlock(sigma2, M, N);

    for (int i=0; i<N; i++)
        ave2[i] = ave[i] * ave[i];  

    vector<double> prog_sum(N);
    vector<double> prog_sum2(N);

    ofstream WriteResults;
    WriteResults.open("01-1-2.out");

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

    rnd.SaveSeed();
    return 0;
}