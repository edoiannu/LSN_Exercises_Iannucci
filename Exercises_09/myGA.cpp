#include "myGA.h"

MyGA :: MyGA(){

    popul.set_size(POP_DIM,N_CITIES);

    vec chromosom;

    for (int i=0; i<POP_DIM; i++){
        chromosom = regspace(1,1,34);
        popul.row(i) = chromosom.t();
    }
}

void MyGA :: set_config(string type){
    int config;
    if (type == "circle") { config = 0; }
    else if (type == "square") { config = 1; }
    _cit.create_config(config);
}

int MyGA :: get_gene(int i, int j){
    return popul(i,j);
}

vec MyGA :: mutation(vec mut, double prob){

    vec mut_(N_CITIES);
    mut_ = mut;

    /*for (int j=0; j<N_CITIES; j++)
        mut(j) = popul(i,j);*/
    
    int mut_type = 0;
    double r = 0.;

    r = _rnd.Rannyu();

    mut_type = int(_rnd.Rannyu(0,4));

    if ( r < prob ){
        if (mut_type == 0) { mut_ = pair_perm(mut); }
        else if (mut_type == 1) { mut_ = shift(mut); }
        else if (mut_type == 2) { mut_ = contig_perm(mut); }
        else if (mut_type == 3) { mut_ = inversion(mut); }
    }

    return mut_;
}

vec MyGA :: pair_perm(vec mut){

    /*vec mut(N_CITIES);

    for (int j=0; j<N_CITIES; j++)
        mut(j) = popul(i,j);*/

    int j = int(_rnd.Rannyu(1,N_CITIES-1));
    int k;

    do{
        k = int(_rnd.Rannyu(1,N_CITIES-1));
    }while( k == j );

    int appo;;

    appo = mut(k);
    mut(k) = mut(j);
    mut(j) = appo;

    return mut;
}

vec MyGA :: shift(vec mut){

    /*vec mut(N_CITIES);

    for (int j=0; j<N_CITIES; j++)
        mut(j) = popul(i,j);*/

    int start = int(_rnd.Rannyu(1,N_CITIES-2));

    int m = int(_rnd.Rannyu(1,N_CITIES-start-1));

    vec appo_m(m);

    for (int j=0; j<m; j++)
        appo_m(j) = mut(start+j);

    int n = int(_rnd.Rannyu(1,N_CITIES-start-m));
    vec appo_n(n);

    for (int j=0; j<n; j++)
        appo_n(j) = mut(start+j+m);

    for (int j=0; j<n; j++)
        mut(start+j) = appo_n(j);

    for (int j=0; j<m; j++)
        mut(start+j+n) = appo_m(j);

    return mut;

}

vec MyGA :: contig_perm(vec mut){

    /*vec mut(N_CITIES);

    for (int j=0; j<N_CITIES; j++)
        mut(j) = popul(i,j);*/

    int m = int(_rnd.Rannyu(1,int(N_CITIES/2)));

    vec appo_m1(m);
    vec appo_m2(m);

    int start1 = int (_rnd.Rannyu(1,N_CITIES-2*m));
    int start2 = int (_rnd.Rannyu(start1+m,N_CITIES-m));
    
    for (int j=0; j<m; j++){
        appo_m1(j) = mut(start1+j);
        appo_m2(j) = mut(start2+j);
    }

    for (int j=0; j<m; j++){
        mut(start1+j) = appo_m2(j);
        mut(start2+j) = appo_m1(j);
    }

    return mut;
}

vec MyGA :: inversion(vec mut){

    /*vec mut(N_CITIES);

    for (int j=0; j<N_CITIES; j++)
        mut(j) = popul(i,j);*/

    int m = int(_rnd.Rannyu(1,N_CITIES));
    int start = int(_rnd.Rannyu(1,N_CITIES-m));

    vec appo(m);
    vec appo_inv(m);

    for (int j=0; j<m; j++)
        appo(j) = mut(start+j);

    for (int j=0; j<m; j++)
        appo_inv(j) = appo(m-1-j);

    for (int j=0; j<m; j++)
        mut(start+j) = appo_inv(j);

    return mut;
}

mat MyGA :: crossover(int i, int j, double prob){

    mat offspring(2,N_CITIES);

    double r = _rnd.Rannyu();

    mat parents(2,N_CITIES);

    for (int k=0; k<N_CITIES; k++){
        parents(0,k) = popul(i,k);
        parents(1,k) = popul(j,k);
    }

    int cut = int(_rnd.Rannyu(2,N_CITIES));

    for (int k=0; k<cut; k++){
        offspring(0,k) = parents(0,k);
        offspring(1,k) = parents(1,k);
    }

    if (r < prob){

        int k = cut;
        int appo;
        bool substitute;

        for (int l=1; l<N_CITIES; l++){
            substitute = true;
            appo = parents(0,l);
            for (int m=1; m<cut; m++){
                if ( offspring(1,m) == appo ) { substitute = false; }
            }
            if ( substitute ) {
                offspring(1,k) = appo;
                k++;
            }
        }

        k = cut;

        for (int l=1; l<N_CITIES; l++){
            substitute = true;
            appo = parents(1,l);
            for (int m=1; m<cut; m++){
                if ( offspring(0,m) == appo ) { substitute = false; }
            }
            if ( substitute ) {
                offspring(0,k) = appo;
                k++;
            }
        }
    }
    else{
        for (int k=cut; k<N_CITIES; k++){
            offspring(0,k) = parents(0,k);
            offspring(1,k) = parents(1,k);
        }
    }

    return offspring;    
}

void MyGA :: substitute(vec mut, int i){
    popul.row(i) = mut.t();
}

bool MyGA :: check(){

    int appo = 0;

    for (int i=0; i<POP_DIM; i++){
        if ( popul(i,0) != 1 ) {
            cout << "ERROR IN MEMBER # " << i+1 << " : the first city must be fixed" << endl;
            return false;
        }
        for (int j=0; j<N_CITIES; j++){
            appo = popul(i,j);
            for (int k=j+1; k<N_CITIES; k++){
                // cout << appo << " " << popul(i,k) << endl;
                if (popul(i,k) == appo) {
                    cout << "ERROR IN MEMBER # " << i+1 << " : your traveling salesman has already been in this city..." << endl;
                    this->print_chromosom(i);
                    return false;
                }
            }
        }
    }

    return true;
}

double MyGA :: loss_function(int i){

    double L = 0.;

    for(int j=0; j<N_CITIES; j++){
        // cout << j << " " << popul(i,j) << " " << popul(i,_cit.pbc(j+1))<< " " << _cit.pbc(j+1) << endl;
        L += pow(_cit.distance(popul(i,j)-1,popul(i,_cit.pbc(j+1))-1),2);
    }

    return L;
}

int MyGA :: selection(/*int p*/){
    
    int j = 0;
    j = int(SURV_CUT*_rnd.Rannyu());

    return j;
}

double MyGA :: fitness(){

    rowvec appo(N_CITIES);

    for (int i=0; i<POP_DIM-1; i++){
        int max = i;
        for (int j=i+1; j<POP_DIM; j++){
            if (this->loss_function(j) < this->loss_function(max)) { max = j; }

            appo = popul.row(i);
            popul.row(i) = popul.row(max);
            popul.row(max) = appo;
        }
    }

    return this->loss_function(0);
}

vec MyGA :: get_chromosom(int i){
    return popul.row(i).t();
}

void MyGA :: print_popul(){
    popul.print();
}

void MyGA :: print_chromosom(int i){
    popul.row(i).print();
}

void MyGA :: print_coords(int i, string filename){

    ofstream outf;
    outf.open(filename);

    for (int j=0; j<=N_CITIES; j++){
        outf << _cit.get_coord(popul(i,_cit.pbc(j))-1,0) << " " << _cit.get_coord(popul(i,_cit.pbc(j))-1,1) << endl;
    }

    outf.close();
}

