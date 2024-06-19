#include "metropolis.h"

const int nparams = 2;

using namespace std;

void Metropolis :: setcoord(int i, double x){
    _coords(i) = x;
}

void Metropolis :: setparam(int i, double x){
    _param(i) = x;
}

void Metropolis :: setenermin(double x){
    _enermin = x;
}

double Metropolis :: get_coord(int i){  
    return _coords(i);
}

void Metropolis :: increase_step(){
    _step_index++;
}

vec Metropolis :: get_vcoord(){
    return _coords;
}

vec Metropolis :: get_vparam(){
    return _param;
}

void Metropolis :: restart(){
    _blk_index = 0;
    _acc_step_sa = 0;
    _nsteps_sa = 0;
    _step_index = 0;
    _acc_step = 0;
}

void Metropolis :: initialize(string inputf, string outputd){

    double coords_temp;
    double param_temp;

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
        /*else if ( property == "TEMPERATURE"){
            input >> _temp;
            _beta = 1. / _temp; 
        }*/
        else if ( property == "SA"){
            input >> _sa;
            _param.resize(nparams);
            for (int i=0; i<nparams; i++){
                input >> param_temp;
                this->setparam(i,param_temp);
            }
            input >> _nsteps_sa;
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

    WriteResults.open(outputd+"/energy.out");

    if (_sa == false) { WriteResults << "# BLOCK" << setw(15) << "CURRENT_AVE" << setw(15) << "ERROR" << setw(15) << "ACCEPTANCE" << endl; }
    else if (_sa == true) {
        WriteResults << "SA_STEP" << setw(15) << "SA_TEMP" << setw(15) << "FINAL_AVE" << setw(15) << "ERROR" << setw(15) << "ACCEPTANCE" << endl;
        WriteParams.open(outputd+"/params.out");
    }
}

void Metropolis :: step(double (*fun)(vec ,vec )){

    double a = 0.;
    double r = 0.;
    vec proposed_coord;
    vec attempted_step;

    proposed_coord.resize(_ndim);
    attempted_step.resize(_ndim);

    if (_dis_type == 0){
        for (int i=0; i<_ndim; i++)
            attempted_step(i) = _rnd.Rannyu(-1.25,1.25);
    }
    else if (_dis_type == 1){
        for (int i=0; i<_ndim; i++)
            attempted_step(i) = _rnd.Gauss(0,1.8);
    }

    proposed_coord = _coords + attempted_step;

    a = min(1., fun(proposed_coord, _param) / fun(_coords, _param));

    r = _rnd.Rannyu();

    this->increase_step();

    if (r <= a ) {
        _coords = proposed_coord;
        _acc_step++;
    }

    if (_step_index % _blksteps == 0) {
        _acceptance(_blk_index) = static_cast<double>(_acc_step) / _blksteps;
        _blk_index++;
        _acc_step = 0;
    }

}

void Metropolis :: sa_step(double (*fun)(vec , vec), double (*fun2)(vec ,vec )){

    double a = 0.;
    double r = 0.;
    double appo = 0.;
    double delta = 0.;

    vec attempted_step;
    vec proposed_param;
    vec proposed_coord(_ndim);
    attempted_step.resize(_param.size());
    proposed_param.resize(_param.size());

    delta = static_cast<double>(3./(_step_index_sa+1));

    attempted_step(0) = _rnd.Rannyu(-1,1) * delta;  // mu
    attempted_step(1) = _rnd.Rannyu(-1,1) * delta;  // sigma

    // cout << attempted_step(0) << " " << attempted_step(1) << endl;

    proposed_param(0) = abs(_param(0) + attempted_step(0));
    proposed_param(1) = abs(_param(1) + attempted_step(1));

    for (int i=0; i<_nsteps; i++){

        for (int j=0; j<_ndim; j++){
            proposed_coord(j) = _coords(j) + _rnd.Rannyu(-1.25,1.25);
            // cout << j << endl;
        }

        a = min(1., fun(proposed_coord, proposed_param) / fun(_coords, proposed_param));

        r = _rnd.Rannyu();

        if (r <= a )
            _coords = proposed_coord;

        appo += fun2(this->get_vcoord(), proposed_param);
    }

    appo /= _nsteps;

    a = min(1., exp(-_beta*appo) / exp(-_beta*_enermin));

    r = _rnd.Rannyu();

    if (r <= a ) {
        _param = proposed_param;
        _acc_step_sa++;
        _enermin = appo;
    }

    _step_index_sa++;
}

void Metropolis :: evolve_param(int i, double delta){
    _param(i) += delta;
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

int Metropolis :: get_nsasteps(){
    return _nsteps_sa;
}

void Metropolis :: settemp(double x){
    _temp = x;
    _beta = 1./_temp;
}

void Metropolis :: print_results(vec ave){

    double final_ave = 0.;
    double final_sigma = 0.;
    double final_acc = 0.;

    vec ave2 = ave % ave;

    vec prog_sum(_nblocks);
    vec prog_sum2(_nblocks);
    vec prog_acc(_nblocks);

    for (int i=1; i<_nblocks; i++){

        for(int j=0; j<=i; j++){

            prog_sum[i] += ave[j];
            prog_sum2[i] += ave2[j];
            prog_acc[i] += _acceptance[j];

        }

        prog_sum[i] /= (i+1);
        prog_sum2[i] /= (i+1);
        prog_acc[i] /= (i+1);

        if (_sa == false) { WriteResults << setw(7) << i << setw(15) << prog_sum[i] << setw(15) << var(prog_sum[i], prog_sum2[i],i) << setw(15) << prog_acc[i] << endl; }

        if (i==_nblocks-1){
            final_ave = prog_sum[i];
            final_sigma = var(prog_sum[i], prog_sum2[i],i);
            final_acc = prog_acc[i];

            if (_sa == true ) {
                WriteResults << setw(7) << _step_index_sa+1 <<setw(15) << _temp << setw(15) << final_ave << setw(15) << final_sigma << setw(15) << final_acc << endl;
                WriteParams << _param(0) << setw(15) << _param(1) << endl;
            }
        }
    }
}

int Metropolis :: get_distype(){
    return _dis_type;
}

double Metropolis :: get_temp(){
    return _temp;
}

void Metropolis :: close_stream(){
    WriteResults.close();
    if (_sa == true){ WriteParams.close(); }
}