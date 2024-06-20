#include "mpi.h"
#include "PT_fun.h"
#include "../libraries/mylib.h"

string provs[N_PROV];
double longs[N_PROV];
double lats[N_PROV];

unsigned int N_GEN = 500; // number of generations
unsigned int N_MIGR = 10;

int main(int argc, char* argv[]){

    ifstream read_provs;
    ifstream read_cap;

    read_provs.open("prov_ita.txt");
    read_cap.open("cap_prov_ita.dat");

    for (int i=0; i<N_PROV; i++){
        read_provs >> provs[i];
        read_cap >> longs[i] >> lats[i];
    }

    read_cap.close();
    read_provs.close();

    int * gene = new int[N_GENE];

    for(int i=0; i<N_IND; i++){
        for(int j=0; j<N_PROV; j++){
            gene[(i*N_PROV)+j] = j+1;
        }
    }

    int size, rank;
    int p1, p2;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //////////// SETTING DIFFERENT PRIMES FOR EACH NODE ////////////

    ifstream read_primes;
    read_primes.open("../libraries/Primes");

    // Check if file opened successfully
    if (!read_primes.is_open()) {
        cerr << "Error opening file." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Skip lines to reach the correct line for this process
    for (int i = 0; i < rank; ++i) {
        read_primes.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    read_primes >> p1 >> p2;

    // read_primes.close();

    int seeds[4];

    //////////// SETTING DIFFERENT SEEDS FOR EACH NODE ////////////

    ifstream read_seeds;
    read_seeds.open("../libraries/Seeds");

    // Check if file opened successfully
    if (!read_seeds.is_open()) {
        cerr << "Error opening file." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Skip lines to reach the correct line for this process
    for (int i = 0; i < rank; ++i) {
        read_seeds.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    for (int i=0; i<4; i++) { read_seeds >> seeds[i]; }

    read_seeds.close();
    read_primes.close();

    Random _rnd;
    _rnd.SetRandom(seeds, p1, p2);

    int * mutated_gene = new int[N_PROV];

    //////////// GENERATING STARTING POPULATION ////////////

    for (int i=0; i<1000; i++){
        for (int l=0; l<N_IND; l++){

            int mut_type = 0;

            mut_type = int(_rnd.Rannyu(0,4));

            if (mut_type == 0) {
                
                int j = int(_rnd.Rannyu(1,N_PROV-1));
                int k = 0;

                do{
                    k = int(_rnd.Rannyu(1,N_PROV-1));
                }while( k == j );

                mutated_gene = pair_perm(gene, l, j, k);
            }
            else if (mut_type == 1) {
                int start = int(_rnd.Rannyu(1,N_PROV-2));
                int m = int(_rnd.Rannyu(1,N_PROV-start-1));
                int n = int(_rnd.Rannyu(1,N_PROV-start-m));
                mutated_gene = shift(gene, l, start, m, n);
            }
            else if (mut_type == 2) {
                int m = int(_rnd.Rannyu(1,int(N_PROV/2)));
                int start1 = int (_rnd.Rannyu(1,N_PROV-2*m));
                int start2 = int (_rnd.Rannyu(start1+m,N_PROV-m));
                mutated_gene = contig_perm(gene, l, m, start1, start2);
            }
            else if (mut_type == 3) {
                int m = int(_rnd.Rannyu(1,N_PROV));
                int start = int(_rnd.Rannyu(1,N_PROV-m));
                mutated_gene = inversion(gene, l, m, start);
            }
            substitute(gene,mutated_gene,l);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        check(gene);
    }

    int j1, j2;
    int * sons = new int[2*N_PROV];
    int * son1 = new int[N_PROV];
    int * son2 = new int[N_PROV];

    double * fit = new double[N_GEN];

    for (int i=0; i<N_GEN; i++){

        int k = SURV_CUT;

        do{
            j1 = int(SURV_CUT*_rnd.Rannyu());  
            do{
                j2 = int(SURV_CUT*_rnd.Rannyu());
            }while( j1 == j2 );

            int cross_trial = _rnd.Rannyu();
            int mut_trial = _rnd.Rannyu();

            if ( cross_trial < 0.9 ){

                int cut = int(_rnd.Rannyu(2,N_PROV));

                sons = crossover(gene, j1, j2, cut);

                for(int j=0; j<N_PROV; j++){
                    son1[j] = sons[j];
                    son2[j] = sons[j+N_PROV];
                }

                substitute(gene, son1, k);
                substitute(gene, son2, k+1);

                if ( mut_trial < 0.1 ){
                    int mut_type = 0;

                    mut_type = int(_rnd.Rannyu(0,4));

                    if (mut_type == 0) {
                        
                        int l = int(_rnd.Rannyu(1,N_PROV-1));
                        int m = 0;

                        do{
                            m = int(_rnd.Rannyu(1,N_PROV-1));
                        }while( m == l );

                        son1 = pair_perm(gene, j1, l, m);
                        son2 = pair_perm(gene, j2, l, m);
                    }
                    else if (mut_type == 1) {
                        int start = int(_rnd.Rannyu(1,N_PROV-2));
                        int m = int(_rnd.Rannyu(1,N_PROV-start-1));
                        int n = int(_rnd.Rannyu(1,N_PROV-start-m));

                        son1 = shift(gene, j1, start, m, n);
                        son2= shift(gene, j2, start, m, n);
                    }
                    else if (mut_type == 2) {
                        int m = int(_rnd.Rannyu(1,int(N_PROV/2)));
                        int start1 = int (_rnd.Rannyu(1,N_PROV-2*m));
                        int start2 = int (_rnd.Rannyu(start1+m,N_PROV-m));

                        son1 = contig_perm(gene, j1, m, start1, start2);
                        son2 = contig_perm(gene, j2, m, start1, start2);
                    }
                    else if (mut_type == 3) {
                        int m = int(_rnd.Rannyu(1,N_PROV));
                        int start = int(_rnd.Rannyu(1,N_PROV-m));

                        son1 = inversion(gene, j1, m, start);
                        son2 = inversion(gene, j2, m, start);
                    }
                    substitute(gene, son1, k);
                    substitute(gene, son2, k+1);
                }          
            }
            k += 2;
        }while( k < N_IND );
        MPI_Barrier(MPI_COMM_WORLD);
        check(gene);

        fit[i] = fitness(longs, lats, gene);

        /*if(rank==7 and i==99){
            for (int l=0; l<N_IND; l++){
                for (int m=0; m<N_PROV; m++)
                    cout << gene[m+l*N_PROV] << " ";
                cout << endl << "loss = " << loss_function(longs, lats, gene, l) << endl;
            }
        }*/

        int ind1, ind2;

        for (int l=0; l<4; l++){
            if (l==3) seeds[l] = i+1;
            else seeds[l] = 0;
        }

        p1 = 2892;
        p2 = 2587;

        _rnd.SetRandom(seeds,p1,p2);

        // if (rank==0) cout << _rnd.Rannyu() << endl;

        if ( (i+1) % 10 == 0 ){
            for (int j=0; j<(size/2)*SURV_CUT; j++){

                // sorteggiamo i processi coinvolti nello scambio
                j1 = int(size*_rnd.Rannyu());  
                do{
                    j2 = int(size*_rnd.Rannyu());
                }while( j1 == j2 );

                // if (rank==2) cout << j1 << " " << j2 << endl;

                // sorteggiamo i due invididui da scambiare fra i processi j1 e j2
                ind1 = int(SURV_CUT*_rnd.Rannyu());
                ind2 = SURV_CUT + int((N_IND-SURV_CUT)*_rnd.Rannyu());

               // if (rank==0) cout << "Scambiamo l'individuo " << ind1 << " del processo " << j1 << " con l'individuo " << ind2 << " del processo " << j2 << endl;

                
                ////// PROVA DI STAMPA ///////

                /*
                cout << "STAMPA DELLE DUE POPOLAZIONI PRIMA DELLO SCAMBIO" << endl;}

                MPI_Barrier(MPI_COMM_WORLD);

                // cout << j1 << endl;

                if (j1==rank){
                    cout << endl << j1 << endl;
                    for (int l=0; l<N_IND; l++){
                        cout << "Ind " << l << ": ";
                        for (int m=0; m<N_PROV; m++)
                            cout << gene[m+l*N_PROV] << " ";
                        cout << endl << "loss = " << loss_function(longs, lats, gene, l) << endl;
                    }
                }
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank==j2){
                    cout << endl << j2 << endl;
                    for (int l=0; l<N_IND; l++){
                        cout << "Ind " << l << ": ";
                        for (int m=0; m<N_PROV; m++)
                            cout << gene[m+l*N_PROV] << " ";
                        cout << endl << "loss = " << loss_function(longs, lats, gene, l) << endl;
                    }
                }*/

                // Sincronizzazione prima dello scambio
                MPI_Barrier(MPI_COMM_WORLD);
                
                int * temp = new int[N_PROV];

                    if (rank == j1) {
                        MPI_Sendrecv(&gene[ind1*N_PROV], N_PROV, MPI_INT, j2, 0, temp, N_PROV, MPI_INT, j2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        for (int l=0; l<N_PROV; l++)
                            gene[l+ind1*N_PROV] = temp[l];
                    } else if (rank == j2) {
                        MPI_Sendrecv(&gene[ind2*N_PROV], N_PROV, MPI_INT, j1, 0, temp, N_PROV, MPI_INT, j1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        for (int l=0; l<N_PROV; l++)
                            gene[l+ind2*N_PROV] = temp[l];
                    }
                delete[] temp;
                MPI_Barrier(MPI_COMM_WORLD);
                // if (rank==0) cout << endl << "SCAMBIO EFFETTUATO" << endl;
                // dopo ogni scambio tutti i processi devono essere sincronizzati
                /*
                if (rank==0){cout << "STAMPA DELLE DUE POPOLAZIONI DOPO LO SCAMBIO" << endl;
                // cout << "POPOLAZIONE " << j1 << endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank==j1){
                    cout << endl << j1 << endl;
                    for (int l=0; l<N_IND; l++){
                        cout << "Ind " << l << ": ";
                        for (int m=0; m<N_PROV; m++)
                            cout << gene[m+l*N_PROV] << " ";
                        cout << endl << "loss = " << loss_function(longs, lats, gene, l) << endl;
                    }
                }
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank==j2){
                    cout << endl << j2 << endl;
                    for (int l=0; l<N_IND; l++){
                        cout << "Ind " << l << ": ";
                        for (int m=0; m<N_PROV; m++)
                            cout << gene[m+l*N_PROV] << " ";
                        cout << endl << "loss = " << loss_function(longs, lats, gene, l) << endl;
                    }
                }*/
            }    
        }

        fit[i] = fitness(longs, lats, gene);

        check(gene);
    }

    /////// MEAN BLOCK ///////

    vector<int> prog_sum(N_GEN,0);
    vector<int> prog_sum2(N_GEN,0);
    vector<int> fit2(N_GEN,0);

    for (int i=0; i<N_GEN; i++)
        fit2[i] = fit[i] * fit[i];

    ofstream print_results[size];

    for (int i=0; i<size; i++){
        if (i==rank){
            print_results[i].open("results/Lmean_proc"+to_string(i)+".out");
            print_results[i] << "BLK #" << setw(15) << "BLK_AVE" << setw(15) << "CURRENT_AVE" << setw(15) << "ERROR" << endl;

            for (int j=0; j<N_GEN; j++){

                if(i==0 and j==0){
                    cout << prog_sum[j] << " " << prog_sum2[j] << endl;
                }

                for(int k=0; k<=j; k++){
                    prog_sum[j] += fit[k];
                    prog_sum2[j] += fit2[k];
                    if(i==0 and j==0){
                        cout << fit[k] << " " << fit2[k] << endl;
                    }
                }

                if(i==0 and j==0){
                    cout << prog_sum[j] << " " << prog_sum2[j] << endl;
                }

            prog_sum[j] /= (j+1);
            prog_sum2[j] /= (j+1);

            print_results[i] << j+1 << setw(15) << fit[j] << setw(15) << prog_sum[j] << setw(15) << var(prog_sum[j], prog_sum2[j],j) << endl;
            }

        }
        MPI_Barrier(MPI_COMM_WORLD);
        print_results[i].close();
    }

    MPI_Finalize();

    delete[] gene;
    delete[] mutated_gene;
    delete[] sons;
    delete[] son1;
    delete[] son2;

    return 0;
}