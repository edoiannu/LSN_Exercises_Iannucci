/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

const string LibrariesDirectory = "/home/edoiannu/Documenti/Lab_simulazione_numerica/LSN_Exercises_Iannucci/libraries/";

Random :: Random(){

   int seed[4];    // array of dimension four containing the four seeds for the random number generator
   int p1, p2;     // definition of the parameters p1 and p2 of the random number generator

   // the seeds and the parameters p1 and p2 are read from file

   ifstream Primes(LibrariesDirectory+"Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
        // cout << p1 << " " << p2 << endl;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

   ifstream input(LibrariesDirectory+"seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
               input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
               //cout << seed[0] << " " << seed[1] << " " << seed[2] << " " << seed[3] << endl;
               SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

}
// Default constructor, does not perform any action

Random :: ~Random(){}
// Default destructor, does not perform any action

void Random :: SaveSeed(){
   // This function saves the current state of the random number generator to a file "seed.out"
   ofstream WriteSeed;
   WriteSeed.open(LibrariesDirectory+"seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << "RANDOMSEED	" << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   // This function generates a random number from a Gaussian distribution with given mean and sigma
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   // This function generates a random number in the range [min, max)
   return min+(max-min)*Rannyu();
}

double Random :: Exp (double lambda){
   // This function generates a random number from an exponential distribution with a given lambda coefficient
   return -(1/lambda) * log(1-Rannyu());
}

double Random :: Cauchy_Lorentz (double mean, double gamma){
   // This function generates a random number from a Cauchy-Lorentz distribution with a given mean and gamma coefficient
   return gamma * tan(M_PI * ( Rannyu() - 0.5 )) + mean;
}

bool Random :: Succ_fail(){
   // This function randmly generate 1 or 0 (success of failure)
   return (round(Rannyu()));
}

double Random :: Theta(){
   // This function generate a generic theta valure between 0 and pi according to the frequency function p(x) = 0.5 * sin(x)
   return acos(1.-2*Rannyu());
}

double Random :: Rannyu(void){
  // This function generates a random number in the range [0,1)
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

double Random :: Distributionf(double (*fun)(double)){
   return fun(Rannyu());
}

void Random :: SetRandom(int * s, int p1, int p2){
  // This function sets the seed and parameters of the random number generator
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

void Random :: SetPrimes(int p1, int p2){
  // This function sets the Primes of the random number generator
   n3 = p1;
   n4 = p2;

   return;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
