#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <vector>
#include "../../../libraries/random.h"
#include "../../../libraries/mylib.h"

using namespace std;

const string ResultsDirectory = "/results/";

const int M = 10000; // number of total throws

int main(){

    Random rnd; // definition of an object of type Random

    int N[4] = {1, 2, 10, 100};
    string Distributions[3] = {"uniform.out", "exponential.out", "Cauchy-Lorentz.out"}; // array with the name of the three distributions
    ofstream WriteResults[3];
    double sum; // helper variable

    for (int i=0; i<3; i++){

        WriteResults[i].open(ResultsDirectory+Distributions[i]);

        for (int j=0; j<4; j++){

            for (int k=0; k<M; k++){

                sum = 0;

                for(int l=0; l<N[j]; l++){
                    if (i==0) { sum += rnd.Rannyu(); }
                    else if (i==1) { sum += rnd.Exp(1); }
                    else if (i==2) { sum += rnd.Cauchy_Lorentz(0, 1); }
                    }

                WriteResults[i] << sum / N[j] << endl;
            }
        }

        WriteResults[i].close();
    }

    rnd.SaveSeed();
    return 0;
}