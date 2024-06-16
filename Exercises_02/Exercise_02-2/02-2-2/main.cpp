#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "/home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/random.h"
#include "/home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/mylib.h"
#include "/home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/randomwalk.h"

using namespace std;

const int M = 10000;        // number of RWs
const int steps = 100;      // number of steps within one RW
const int N = 100;          // number of blocks
const int L = M / N;        // number of RWs within one block
const int a = 1;            // lattice spacing

const string ResultsDirectory = "/home/edoiannu/Documenti/Lab_simulazione_numerica/Exercises_02/results/";

// facciamo una classe "random walk"
// membri privati: le coordinate spaziali
// membri publici: funzioni che incrementino le coordinate spaziali

int main(){

    Random rnd;

    vector<vector<double>> phi(M, vector<double>(steps));
    vector<vector<double>> theta(M, vector<double>(steps));
    vector<vector<double>> ave(N, vector<double>(steps));
    vector<vector<double>> ave2(N, vector<double>(steps));

    for (int i=0; i<M; i++){
        for(int j=0; j<steps; j++){
            phi[i][j] = rnd.Rannyu(0,2*M_PI);
            theta[i][j] = rnd.Theta();
        }
    }

    // crea un array di tipo "random walk" di M (numbero di rws) elementi

    vector<RandomWalk> rw(M);   // M rws
    double sum;
    double sum2;

    ofstream WriteResults;
    WriteResults.open(ResultsDirectory+"02-2-2.out");

    for (int i=0; i<steps; i++){ // ciclo sui passi
    
        for (int j=0; j<N; j++){    // ciclo sui blocchi

            sum = 0;

            for(int k=j*L; k<(j+1)*L; k++){
                rw[k].Evolve_polar(phi[k][i], theta[k][i], a);
                sum += rw[k].Getr2();}

            ave[j][i] = sum / L;                    // <|r_i|^2>    dove i indicizza il blocco
            ave2[j][i] = ave[j][i] * ave[j][i];     // <|r_i|^2>^2
        }

        sum = 0;
        sum2 = 0;

        for (int j=0; j<N; j++){
            sum += ave[j][i];
            sum2 += ave2[j][i];
        }

        sum /= N;
        sum2 /= N;

        WriteResults << sqrt(sum) << setw(12) << 0.5 * (1 / sqrt(sum)) * var(sum, sum2, N) << endl;
    }

    WriteResults.close();

    return 0;
}