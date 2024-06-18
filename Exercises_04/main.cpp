#include <iostream>
#include "../libraries/NSL_SIMULATOR/SOURCE/system.h"

const string inputd = "INPUT/INPUT-DATA/";

using namespace std;

int main(int argc, char *argv[]){

    if (argc != 2 and (argv[0] != "SOLID" or argv[0] != "LIQUID" or argv[0] != "GAS")){
    cout << "Usage: " << argv[0] << " <SOLID/LIQUID/GAS>" << endl;
    exit(-1);
  }

  string inputf = inputd+string(argv[1])+".dat";
  string outputd = "results/"+string(argv[1])+"/";

  int nconf = 1;
  System SYS;
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
    SYS.block_reset(i+1, outputd);
  }
  SYS.finalize(outputd);

  return 0;
}