#include<iostream>
#include "parameter.h"
#include "monte.h"
#include <cmath>
#include <vector>
#include <list>
#include <fstream>
int main(){
	//***********************initialize the system********units SI**************************//
    int size=20;
    int N=size*size;
    std::vector<atom> atomall(N);
    int temp;
    double r_verlet=1.2*r_cut;
    double r_shell_initial=r_verlet-r_cut;
    double r_shell=r_shell_initial;
    std::vector<double> totalstress;
    std::vector<double> initial(2,0);
    for(size_t i=0;i<size;i++)
     for(size_t j=0;j<size;j++){
        temp=i*size+j;
        atomall[temp].setx(i*r_min);
        atomall[temp].sety(j*r_min);
        atomall[temp].setr(0.1*r_min);
        atomall[temp].setm(39.948*1e-3/NA);
        atomall[temp].setv(initial);
        atomall[temp].setf(initial);
       }
    //**********************end initialize ****************************************
   montecarlo(atomall,0.1*r_min,5,2000000);
}
