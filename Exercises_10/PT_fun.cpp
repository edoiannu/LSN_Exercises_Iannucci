#include "PT_fun.h"

void substitute(int * v, int * vnew, int j){
    for (int i=0; i<N_PROV; i++){
        v[i+j*N_PROV] = vnew[i];
    }
}

int * pair_perm(int *v, int i, int j, int k){

    int * mut = new int[N_PROV];

    for (int l=0; l<N_PROV; l++)
        mut[l] = v[l+i*N_PROV];

    int appo;

    appo = mut[k];
    mut[k] = mut[j];
    mut[j] = appo;

    return mut;
}

int * shift(int *v, int i, int start, int m, int n){

    int * mut = new int[N_PROV];

    for (int l=0; l<N_PROV; l++)
        mut[l] = v[l+i*N_PROV];

    int * appo_m = new int[m];

    for (int j=0; j<m; j++)
        appo_m[j] = mut[start+j];

    int * appo_n = new int[n];

    for (int j=0; j<n; j++)
        appo_n[j] = mut[start+j+m];

    for (int j=0; j<n; j++)
        mut[start+j]= appo_n[j];

    for (int j=0; j<m; j++)
        mut[start+j+n] = appo_m[j];

    delete[] appo_m;
    delete[] appo_n;

    return mut;
}

int * contig_perm(int *v, int i, int m, int start1, int start2){

    int * appo_m1 = new int[m];
    int * appo_m2 = new int[m];

    int * mut = new int[N_PROV];

    for (int l=0; l<N_PROV; l++)
        mut[l] = v[l+i*N_PROV];
    
    for (int j=0; j<m; j++){
        appo_m1[j] = mut[start1+j];
        appo_m2[j] = mut[start2+j];
    }

    for (int j=0; j<m; j++){
        mut[start1+j] = appo_m2[j];
        mut[start2+j] = appo_m1[j];
    }

    delete[] appo_m1;
    delete[] appo_m2;

    return mut;
}

int * inversion(int *v, int i, int m, int start){

    int * mut = new int[N_PROV];

    for (int l=0; l<N_PROV; l++)
        mut[l] = v[l+i*N_PROV];

    int * appo = new int[m];
    int * appo_inv = new int[m];

    for (int j=0; j<m; j++)
        appo[j] = mut[start+j];

    for (int j=0; j<m; j++)
        appo_inv[j] = appo[m-1-j];

    for (int j=0; j<m; j++)
        mut[start+j] = appo_inv[j];

    delete[] appo;
    delete[] appo_inv;

    return mut;
}

int * crossover(int * v, int i, int j, int cut){

    int * offspring = new int[2*N_PROV];

    int * parents = new int[2*N_PROV];

    for (int k=0; k<N_PROV; k++){
        parents[k] = v[i*N_PROV+k];
        parents[k+N_PROV] = v[j*N_PROV+k];
    }

    for (int k=0; k<cut; k++){
        offspring[k] = parents[k];
        offspring[k+N_PROV] = parents[k+N_PROV];
    }

    int k = cut;
    int appo;
    bool substitute;

    for (int l=1; l<N_PROV; l++){
        substitute = true;
        appo = parents[l];
        for (int m=1; m<cut; m++){
            if ( offspring[m+N_PROV] == appo ) { substitute = false; }
        }
        if ( substitute ) {
            offspring[k+N_PROV] = appo;
            k++;
        }
    }

    k = cut;

    for (int l=1; l<N_PROV; l++){
        substitute = true;
        appo = parents[l+N_PROV];
        for (int m=1; m<cut; m++){
            if ( offspring[m] == appo ) { substitute = false; }
        }
        if ( substitute ) {
            offspring[k] = appo;
            k++;
        }
    }

    return offspring;
}

int pbc(int i){
    if (i < N_PROV) { return i; }
    else if (i >= 34) { return i-N_PROV; }
}

double distance(double x1, double y1, double x2, double y2){
    double r2 = 0.;

    r2 = pow((x1-x2),2) + pow((y1-y2),2);

    return sqrt(r2);
}

double loss_function(double x[], double y[], int * v, int i){

    // cout << "INDIVIDUO " << i << endl;

    double L = 0.;
    double appo = 0.;

    for(int j=i*N_PROV; j<(i+1)*N_PROV; j++){
        appo = pow(distance(x[v[j]-1],y[v[j]-1],x[pbc(v[j+1])-1],y[pbc(v[j+1])-1]),2);
        // cout << "Distanza fra " << v[j] << " e " << pbc(v[j+1]) << ": " << appo << endl;
        // cout << "Distanza fra " << v[j] << " e " << pbc(v[j+1]) << ": " << appo << endl;
        L += appo;
        // cout << "Progressivo: " << L << endl;
    }

    return L;
}

bool check(int *v){

    int appo = 0;

    for (int i=0; i<N_IND; i++){
        if ( v[i*N_PROV] != 1 ) {
            cout << "ERROR IN MEMBER # " << i+1 << " : the first city must be fixed" << endl;
            return false;
        }
        for (int j=i*N_IND; j<i*N_IND+N_PROV; j++){
            appo = v[j];
            for (int k=j+1; k<N_PROV; k++){
                // cout << appo << " " << popul(i,k) << endl;
                if (v[k] == appo) {
                    cout << "ERROR IN MEMBER # " << i+1 << " : your traveling salesman has already been in this city..." << endl;
                    return false;
                }
            }
        }
    }

    return true;
}

double fitness(double x[], double y[], int * v){

    int * appo = new int[N_PROV];

    for (int i=0; i<N_IND-1; i++){
        int max = i;
        for (int j=i+1; j<N_IND; j++){
            if (loss_function(x, y, v, j) < loss_function(x, y, v, max)) { max = j; }
        }

        for (int l=0; l<N_PROV; l++)
            appo[l] = v[l+i*N_PROV];

        for (int l=0; l<N_PROV; l++)
            v[l+i*N_PROV] = v[l+max*N_PROV];
        
        for (int l=0; l<N_PROV; l++)
            v[l+max*N_PROV] = appo[l];
    }

    delete[] appo;

    return loss_function(x, y, v, 0);
}

void print_genome(int * v){
    for(int i=0; i<N_IND; i++){
        for(int j=0; j<N_PROV; j++){
            cout << v[(i*N_PROV)+j] << " ";
        }
        cout << endl;
    }
}