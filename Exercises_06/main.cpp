#include <iostream>
#include <string>
#include "../libraries/NSL_SIMULATOR/SOURCE/system.h"

const string inputd = "INPUT/INPUT-DATA/";
const string inputf = inputd+"input.ISING";

using namespace std;

/*bool starting_condition(string sim_type, double t0){
  bool _sim_type = true;
  bool _t0 = true;

  for (double i = 0.5; i <= 2.0; i += 0.1){
    if (fabs(t0 - i) < 0.000001){
      _t0 = false;
      break;
    }
  }

  if (sim_type != "METROPOLIS" and sim_type != "GIBBS"){
    _sim_type = false;
  }

  return _sim_type and _t0;
}*/

int main(int argc, char *argv[]){

  if (argc != 3 /*or starting_condition(argv[1], stod(argv[2])) == false*/){
    cout << "Usage: " << argv[0] << " <METROPOLIS/GIBBS> <TEMPERATURE>" << endl;
    exit(-1);
  }

  string outputd = "results/"+string(argv[1])+"/"+string(argv[2])+"/";

  int nconf = 1;
  System SYS;
  SYS.set_temp(stod(argv[2]));
  SYS.initialize(inputf, outputd);
  SYS.initialize_properties(outputd);
  SYS.block_reset(0, outputd);

  for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks

    for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      SYS.step();
      SYS.measure();
      if(j%10 == 0){
  //        SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf++;
      }
      
    }
    SYS.averages(i+1, outputd);
    if (i+1 == SYS.get_nbl()) { SYS.print_exercise6("results/"+string(argv[1])+"/",i+1); }
    SYS.block_reset(i+1, outputd);
  }
  SYS.finalize(outputd);

  return 0;
}