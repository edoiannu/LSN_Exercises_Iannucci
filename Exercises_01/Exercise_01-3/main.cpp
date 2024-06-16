#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <vector>
#include "/home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/random.h"
#include "/home/edoiannu/Documenti/Lab_simulazione_numerica/libraries/mylib.h"

using namespace std;

const string ResultsDirectory = "/home/edoiannu/Documenti/Lab_simulazione_numerica/Exercises_01/results/";

const double l = 0.8;   // neddle length
const double d = 1.;    // distance between plane lines
const int M = 100000;    // number of throws
const int N = 100;      // number of blocks
const int L = M / N;    // number of throws in each block

int main(){

    Random rnd;
    vector<double> ave(N);
    vector<double> ave2(N);

    double ycm, x, y;
    double tan, y1, y2;
    int Nhit;

    for (int i=0; i<N; i++){

        Nhit = 0;

        for (int j=0; j<L; j++){

            ycm = rnd.Rannyu(-0.3,20.7);

            do{
                x = rnd.Rannyu(-1,1);
                y = rnd.Rannyu(-1,1);
            }while(x*x+y*y < 1);

            tan = atan(x/y);

            y1 = ycm + (l/2) * sin(tan);
            y2 = ycm - (l/2) * sin(tan);

            if (floor(y1) != floor(y2)) { Nhit++; }

            // if (i==0) { cout << y1 << " " << y2 << " " << Nhit << endl; }
        }

        ave[i] = (2*l*L) / (Nhit*d);
        ave2[i] = ave[i] * ave[i];
    }

    vector<double> prog_sum(N);
    vector<double> prog_sum2(N);

    ofstream WriteResults;
    WriteResults.open(ResultsDirectory+"buffon.out");

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

    rnd.SaveSeed();
    return 0;
}