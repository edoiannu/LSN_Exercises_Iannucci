#include "myGA.h"

using namespace std;

const int N = 100;  // number of generation
const int M = 100;
const int L = M / N;

int main(int argc, char *argv[]){

    if (argc != 2 and (argv[1] != "square" or argv[1] != "circle")){
        cout << "Usage: " << argv[0] << " <circle/square>" << endl;
        exit(-1);
    }

    string config = argv[1];
    Cities city;
    MyGA TSP;
    TSP.set_config(config);

    bool check = 0;

    vec candidate(N_CITIES);

    for (int i=0; i<1000; i++){
        cout << "Mutation cycle # " << i+1 <<" started" << endl;
        for (int j=0; j<POP_DIM; j++){
            candidate = TSP.mutation(TSP.get_chromosom(j), 1.);
            TSP.substitute(candidate, j);
            }

        check = TSP.check();
        if ( check == false ) { break; }
        else if ( check == true ) { cout << "Mutation cycle # " << i+1 <<" completed" << endl; }
    }

    cout << "Your starting population is ready!" << endl;

    mat sons(2,N_CITIES);
    int j1, j2;
    vec bestL(N);
    vec Lmean(N);
    vec Lmean2(N);
    vec Lmean_err(N);
    double appo = TSP.fitness();
    double appo_ = 0.0;

    int k;

    for (int i=0; i<N; i++){
        appo = 0;
        for (int j=i*L; j<(i+1)*L; j++){
            k = SURV_CUT;
            do{
                j1 = TSP.selection();  
                do{
                    j2= TSP.selection();
                }while( j1 == j2 );
                sons = TSP.crossover(j1, j2, 0.9);
                candidate = TSP.mutation(sons.row(0).t(),0.1);
                TSP.substitute(candidate,k);
                candidate = TSP.mutation(sons.row(1).t(),0.1);
                TSP.substitute(candidate,k+1);
                TSP.check();
                k += 2;
            }while( k < POP_DIM );
            appo += TSP.fitness();
            cout << "Generation # " << j+1 << " created!" << endl;
            for (int l=0; l<POP_DIM/2; l++){
                appo_ = TSP.loss_function(l);
                Lmean(i) += appo_;
                Lmean2(i) += appo_ * appo_;
            }
            Lmean(i) /= (POP_DIM/2);
            Lmean2(i) /= (POP_DIM/2);
            Lmean_err(i) = var(Lmean(i), Lmean2(i),i);
        }
        bestL(i) = appo / L;
    }
    /*
    vec ave2(N);
    ave2 = ave%ave;

    vec prog_sum(N);
    vec prog_sum2(N);
    */
    ofstream print_results;
    print_results.open("results/"+config+"/L.out");

    print_results << "GEN #" << setw(15) << "BEST_L" << setw(15) << "L_MEAN" << setw(15) << "ERROR" << endl;

    for (int i=0; i<N; i++){
        print_results << i+1 << setw(15) << bestL(i) << setw(15) << Lmean(i) << setw(15) << Lmean_err(i) << endl;
    }

    cout << "Your best path:" << endl;
    TSP.print_chromosom(0);
    cout << "L2 = " << TSP.loss_function(0) << endl;
    

    TSP.print_coords(0,"results/"+config+"/best_path.out");

    print_results.close();

    return 0;
}