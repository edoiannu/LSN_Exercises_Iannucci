#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <vector>
#include "/home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/random.h"
#include "/home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/mylib.h"

using namespace std;

const string ResultsDirectory = "/home/edoiannu/Documenti/Lab_simulazione_numerica/Exercises_01/Exercise_01-1/01-1-3/results/";

const int M = 100;      // number of subintervals of [0,1]
const int N = 100;      // number of groups of throws
const int n = 10000;    // number of throws in each group

int main(){

    Random rnd;     // definition of an object of type Random

    // ofstream prova;
    // prova.open("prova.txt");

    vector<vector<double>> r(N, vector<double>(n));

    for (int i=0; i<N; i++)
       for (int j=0; j<n; j++)
            r[i][j] = rnd.Rannyu();
    
    //                      i               j
    vector<vector<int>> n_i(N,  vector<int>(M));     // number of throws of the (i+1)-th generation fallen in the (j+1)-th sub-interval
    vector<double> chi2(N);

    for (int i=0; i<N; i++){            // iteration over the generations
        for (int j=0; j<n; j++){        // throws in the i-th generation
            for (int k=0; k<M; k++){    // iteration over the M subintervals
                if (r[i][j] >= k * (1./M) && r[i][j] < (k + 1) * (1./M)){
                    n_i[i][k]++;}
                }
            }
        }

    for (int i=0; i<N; i++){
        for (int j=0; j<M; j++){
            chi2[i] += pow((n_i[i][j] - n / M),2) / (n / M);
            }
        }

    ofstream WriteResults;
    WriteResults.open(ResultsDirectory+"01-1-3.out");

    for (int i=0; i<N; i++)
        WriteResults << chi2[i] << endl;
    
    WriteResults.close();

    rnd.SaveSeed();
    return 0;
}