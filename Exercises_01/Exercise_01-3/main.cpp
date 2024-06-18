#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <vector>
#include "../../libraries/random.h"
#include "../../libraries/mylib.h"

using namespace std;

const double l = 1.;   // neddle length
const double d = 2.;    // distance between plane lines
const int M = 100000000;    // number of total throws
const int N = 100;      // number of block
const int L = M / N;    // number of throws in each block

int main(){

    Random rnd;
    vector<double> pi(N);
    vector<double> ave(N);
    vector<double> ave2(N);

    double ycm, x, y;
    double theta;
    int Nhit;

    for (int i=0; i<N; i++){

        Nhit = 0;

        for (int j=0; j<L; j++){

            ycm = rnd.Rannyu(-d*0.5,d*0.5);

            do{
                x = rnd.Rannyu(-1,1);
                y = rnd.Rannyu(-1,1);
            }while(x*x+y*y > 1);

            theta = atan(x/y);

            if ( abs(0.5*l*sin(theta))-abs(ycm)>0 ) {
               Nhit += 1;
            }
        }

        pi[i] = 2* l/d * L/Nhit;
        cout << pi[i] << endl;
    }

    for (int i=0; i<N; i++)
        ave2[i] = pi[i] * pi[i];

    vector<double> prog_sum(N);
    vector<double> prog_sum2(N);

    ofstream WriteResults;
    WriteResults.open("buffon.out");

    WriteResults << "# STEPS " << M << endl;
    WriteResults << "# BLOCKS " << N << endl;

    for (int i=0; i<N; i++){

        for(int j=0; j<=i; j++){
            prog_sum[i] += pi[j];
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