#include "tJSBMF.h"
#include <iostream>

int main(){
    double J = 0.1;
    double x = 0.1; //hole dopping
    double T = 0;
    int N = 10;
    char wave = 'd';
    tJSBMF example(J, x, T, N, wave);
    example.h = 0.1;
    example.B = 0.1;
    example.Delta = 0.5;
    example.step_forward()
}