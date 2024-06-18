#include "metropolis.h"

using namespace std;

void Metropolis :: setcoord(int i, double x){
    _coords(i) = x;
}

double Metropolis :: get_coord(int i){  
    return _coords(i);
}

void Metropolis :: increase_step(){
    _step_index++;
}

void Metropolis :: initialize(string inputf){

    double coords_temp;

    ifstream input(inputf);
    while ( !input.eof() ){
        input >> property;
        if ( property == "DIMENSION" ){
            input >> _ndim;
            _coords.resize(_ndim);
        }
        else if ( property == "STARTING_POSITION" ){
            for (int i=0; i<_ndim; i++){
                input >> coords_temp;
                this->setcoord(i,coords_temp);
            }
        }
        else if ( property == "DISTRIBUTION"){
            input >> _dis_type;
        }
        else if ( property == "STEPS"){
            input >> _nsteps;
        }
        else if ( property == "BLOCKS"){
            input >> _nblocks;
        }
        else if ( property == "ENDINPUT"){
            cout << "Reading input completed!" << endl;
            break;
        }
    }
    input.close();
    _blksteps = _nsteps / _nblocks;
    _acceptance.resize(_nblocks);
}

bool Metropolis :: step(double (*fun)(double , double , double )){

    double a = 0.;
    double r = 0.;
    vec proposed_coord;
    vec attempted_step;

    proposed_coord.resize(_ndim);
    attempted_step.resize(_ndim);

    if (_dis_type == 0){
        for (int i=0; i<3; i++)
            attempted_step(i) = _rnd.Rannyu(-3,3);
    }
    else if (_dis_type == 1){
        for (int i=0; i<3; i++)
            attempted_step(i) = _rnd.Gauss(0,0.8);
    }

    proposed_coord = _coords + attempted_step;

    a = min(1., fun(proposed_coord(0), proposed_coord(1), proposed_coord(2)) / fun(_coords(0), _coords(1), _coords(2)));

    r = _rnd.Rannyu();

    this->increase_step();

    if (_step_index % _blksteps == 0/**/) {
        // cout << _step_index << endl;
        _acceptance(_blk_index) = static_cast<double>(_acc_step) / _blksteps;
        // cout << _acc_step << " " <<_blksteps << endl;
        _blk_index++;
        _acc_step = 0;
    }

    if (r <= a ) {
        _coords = proposed_coord;
        _acc_step++;
        // if (_blk_index == 2 ) cout << _step_index << " " << _acc_step << endl;

        return true;
    }

    return false;

}
int Metropolis :: get_nsteps(){
    return _nsteps;
}

int Metropolis :: get_nblocks(){
    return _nblocks;
}

int Metropolis :: get_blksteps(){
    return _blksteps;
}

void Metropolis :: print_results(vec ave, string outputf){

    ofstream WriteResults;
    WriteResults.open(outputf);

    vec ave2 = ave % ave;

    vec prog_sum(_nblocks);
    vec prog_sum2(_nblocks);

    WriteResults << "# BLOCK" << setw(10) << "OBS_AVE" << setw(15) << "ERROR" << setw(15) << "ACCEPTANCE" << endl;

    for (int i=0; i<_nblocks; i++){

        for(int j=0; j<=i; j++){

            prog_sum[i] += ave[j];
            prog_sum2[i] += ave2[j];

        }

        prog_sum[i] /= (i+1);
        prog_sum2[i] /= (i+1);

        WriteResults << setw(7) << i << setw(10) << prog_sum[i] << setw(15) << var(prog_sum[i], prog_sum2[i],i) << setw(15) << _acceptance[i] << endl;

    }

    WriteResults.close();
}

int Metropolis :: get_distype(){
    return _dis_type;
}