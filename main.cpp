#include "tJSBMF.h"
#include <iostream>

int main(){
    double J = 0.4;
    double x = 0.15; //hole dopping
    double T = 0;
    int N = 100;
    char wave = 'd';
    tJSBMF example(J, x, T, N, wave);
    example.h = 0.5;
    example.B = 1;
    example.Delta = -10;
    example.self_consistent();
    // std::cout<<-N/2<<std::endl;
    return 0;
}