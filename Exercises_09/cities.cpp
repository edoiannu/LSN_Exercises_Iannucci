#include "cities.h"

Cities :: Cities(){
    pos.set_size(N_CITIES,2);

    /*vec angles(N_CITIES);

    angles(0) = 0.;

    for (int i=1; i<N_CITIES; i++)
        angles(i) = _rnd.Rannyu(0,2*M_PI);

    angles = sort(angles);

    for (int i=0; i<N_CITIES; i++){
        pos(i,0) = cos(angles(i));
        pos(i,1) = sin(angles(i));
    }*/
}

void Cities :: create_config(int i){

    if (i==0){
        vec angles(N_CITIES);

        angles(0) = 0.;

        for (int i=1; i<N_CITIES; i++)
            angles(i) = _rnd.Rannyu(0,2*M_PI);

        angles = sort(angles);

        for (int i=0; i<N_CITIES; i++){
            pos(i,0) = cos(angles(i));
            pos(i,1) = sin(angles(i));
        }
    }
    else if (i==1){
        for (int i=0; i<N_CITIES; i++){
            pos(i,0) = _rnd.Rannyu(-1,1);
            pos(i,1) = _rnd.Rannyu(-1,1);
        }
    }
}

double Cities :: distance(int i, int j){
    double r2 = 0.;

    for (int k = 0; k<pos.n_cols; k++)
        r2 += pow((pos(i,k)-pos(j,k)),2);

    return sqrt(r2);
}

int Cities :: pbc(int i){
    if (i < 34) { return i; }
    else if (i >= 34) { return i-N_CITIES; }
}

double Cities :: get_coord(int i, int j){
    return pos(i,j);
}