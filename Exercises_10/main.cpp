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

    //////////// PRINT PARAMS //////////////////////////////////////

    ofstream print_par;
    print_par.open("params.out");

    print_par << "N_IND " << N_IND << endl;
    print_par << "N_IND " << N_GEN << endl;
    print_par << "N_IND " << N_MIGR << endl;

    print_par.close();

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
        if (rank==0) cout << "Starting population generation: cycle " << i+1 << " completed" << endl;
    }

    int j1, j2;
    int * sons = new int[2*N_PROV];
    int * son1 = new int[N_PROV];
    int * son2 = new int[N_PROV];

    double * fit = new double[N_GEN];
    double * Lmean = new double[N_GEN];
    int * best_path = new int[N_PROV]; // array containing the order of the provinces of the best path
    double * best_fits = new double[N_GEN]; // array containing the best fit for each generation

    for (int i=0; i<N_GEN; i++){

        if (rank==0) cout << "Generation " << i+1 << endl;

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

        // NIENTE SCAMBI
        /*
        if ( (i+1) % 10 == 0 ){

            int ind1, ind2;

            for (int l=0; l<4; l++){
                if (l==3) seeds[l] = i+1;
                else seeds[l] = 0;
            }

            p1 = 2892;
            p2 = 2587;

            _rnd.SetRandom(seeds,p1,p2);

            for (int j=0; j<(size/2)*SURV_CUT; j++){

                // sorteggiamo i processi coinvolti nello scambio
                j1 = int(size*_rnd.Rannyu());  
                do{
                    j2 = int(size*_rnd.Rannyu());
                }while( j1 == j2 );

                // sorteggiamo i due invididui da scambiare fra i processi j1 e j2
                ind1 = int(SURV_CUT*_rnd.Rannyu());
                ind2 = SURV_CUT + int((N_IND-SURV_CUT)*_rnd.Rannyu());

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
            }

            // RISETTIAMO I SEMI DIVERSI PER CIASCUN NODO

            if (rank==0){
                for (int l=0; l<4; l++){
                    if (l==0) seeds[l] = i+1;
                    else seeds[l] = 0;
                }
            }
            else if (rank==1){
                for (int l=0; l<4; l++){
                    if (l==1) seeds[l] = i+1;
                    else seeds[l] = 0;
                }
            }
            else if (rank==2){
                for (int l=0; l<4; l++){
                    if (l==0) seeds[l] = i+1;
                    else if (l==1) seeds[l] = i+1;
                    else seeds[l] = 0;
                }
            }
            else if (rank==3){
                for (int l=0; l<4; l++){
                    if (l==2) seeds[l] = i+1;
                    else seeds[l] = 0;
                }
            }
            else if (rank==4){
                for (int l=0; l<4; l++){
                    if (l==2) seeds[l] = i+1;
                    else seeds[l] = 0;
                }
            }
            else if (rank==5){
                for (int l=0; l<4; l++){
                    if (l==1) seeds[l] = i+1;
                    else if (l==2) seeds[l] = i+1;
                    else seeds[l] = 0;
                }
            }
            else if (rank==6){
                for (int l=0; l<4; l++){
                    if (l!=3) seeds[l] = i+1;
                    else seeds[l] = 0;
                }
            }
            else if (rank==7){
                for (int l=0; l<4; l++){
                    if (l==1) seeds[l] = i+1;
                    else if (l==3) seeds[l] = i+1;
                    else seeds[l] = 0;
                }
            }

            _rnd.SetRandom(seeds,p1,p2);
        }
        */
        fit[i] = fitness(longs, lats, gene);

        double appo = 0.0;

        for (int j=0; j<(N_IND)/2; j++){
            appo += loss_function(longs, lats, gene, j);
        }

        Lmean[i] = appo / ((N_IND)/2);

        check(gene);

        MPI_Barrier(MPI_COMM_WORLD);
        
        double * temp = new double[size];

        if (rank == 0){
            temp[0] = fit[i];
        }

        appo = 0.;
        
        for (int j=1; j<size; j++){
            if (j == rank){
                appo = fit[i];
                MPI_Send(&appo, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
            else if (rank == 0){
                MPI_Recv(&temp[j],1,MPI_DOUBLE,j,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        
        int best = 0;

        if (rank==0){

            best = 0;
            double best_el = temp[best];

            for (int j=1; j<size; j++){

                if (temp[j] < best_el) {
                    best_el = temp[j];
                    best = j;
                }
            }

        }    

        for (int j=1; j<size; j++){
            if (rank == 0){
                MPI_Send(&best, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
            }
            else if (rank == j){
                MPI_Recv(&best,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

        appo = 0.;
        
        if (best != 0){
            if (rank == best){
                appo = fit[i];
                MPI_Send(&appo, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
            else if (rank == 0){
                MPI_Recv(&best_fits[i],1,MPI_DOUBLE,best,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        else{
            if (rank == 0) { best_fits[i] = fit[i]; }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if ( i == N_GEN-1 ){
            if (best != 0){
            if (rank == best){
                appo = fit[i];
                MPI_Send(&gene[0], N_PROV, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
            else if (rank == 0){
                MPI_Recv(&best_path[0],N_PROV,MPI_INT,best,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            }
            else{
                if (rank == 0) {
                    for (int j=0; j<N_PROV; j++){
                        best_path[j] = gene[j];
                    }
                }
            }
        }
        delete[] temp;
    }
    /////// PRINT RESULTS ///////
    
    ofstream print_results[size];
    // ofstream print_bestfits;
    // ofstream print_bestpath;

    // print_bestfits.open("results/Lbest.out");
    // print_bestpath.open("results/best_path.out");

    for (int i=0; i<size; i++){
        if (i==rank){
            print_results[i].open("results/L_proc"+to_string(i)+"_nomigr.out");
            print_results[i] << "GEN #" << setw(15) << "BEST_L" << setw(15) << "L_MEAN"  << endl;

            for (int j=0; j<N_GEN; j++){

                print_results[i] << j+1 << setw(15) << fit[j] << setw(15) << Lmean[j] << endl;
            }

        }
        MPI_Barrier(MPI_COMM_WORLD);
        print_results[i].close();
    }

    /*
    if (rank==0){
        for (int i=0; i<N_GEN; i++) print_bestfits << best_fits[i] << endl;
        for (int i=0; i<N_PROV; i++) print_bestpath << longs[best_path[i]-1] << setw(10) << lats[best_path[i]-1] << endl;
        print_bestpath << longs[0] << setw(10) << lats[0] << endl;
    }
    */
    

    MPI_Finalize();

    read_seeds.close();
    read_primes.close();
    // print_bestfits.close();
    // print_bestpath.close();

    delete[] gene;
    delete[] mutated_gene;
    delete[] sons;
    delete[] son1;
    delete[] son2;
    delete[] Lmean;
    delete[] best_fits;
    delete[] best_path;

    return 0;

}