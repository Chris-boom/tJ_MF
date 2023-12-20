#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include "tJSBMF.h"

tJSBMF::tJSBMF(double J, double x, double T, int N, char wave,
    int maxstep, double atol, double rtol): J(J), x(x), T(T), N(N), wave(wave),
    maxstep(maxstep), atol(atol), rtol(rtol){
    //Initialize
    if(wave == 's'){
        wf = 1;
    }
    else if(wave == 'd'){
        wf = -1;
    }
    else{
        std::cout<<"wrong symmetry!"<<std::endl;
        exit(2333);
    }
    NN = N*N;
    Delta = 0;
    h = 0;
    B = 0;
    lamb = 0;
    mu = 0;
}
double tJSBMF::fF(const double x, const double Temp){
    //Fermi distribution function
    if(Temp == 0){
        if(x > 0){
            return 0;
        }
        else if (x == 0){
            return 0.5;
        }
        else{
            return 1;
        }
    }
    else{
        return 1./(exp(x/Temp)+1);
    }
}
double tJSBMF::fB(const double x, const double Temp){
    //Boson distribution. x must be positive!!!
    if(x <= 0){
        std::cout<<"Error! Bose distribution undefined."<<std::endl;
        exit(2333333);
    }
    else if(Temp == 0){
        return 0;
    }
    else{
        return 1./(exp(x/Temp)-1);
    }
}
double tJSBMF::holon_numb(double lamb){
    //only when T!=0, or lamb must be the lowest energy level.
    double numholon = 0;
    if(T==0){
        std::cout<<"can not be used under nonzero temperature."<<std::endl;
        exit(22333);
    }
    for(int xi=-N/2; xi<N/2; xi++){
        for(int yi=-N/2; yi<N/2; yi++){
            numholon += fB(-2*B*(cos(ki(xi))+cos(ki(yi)))-lamb, T);
        }
    }
    numholon /= NN;
    return numholon;
}
double tJSBMF::spinon_numb(double mu){
    //calculate the spinon number
    double numspinon = 0;
    double ek=0;
    double Ek=0;
    double Delta_k=0;
    for(int xi=-N/2; xi<N/2; xi++){
        for(int yi=-N/2; yi<N/2; yi++){
            Delta_k = Delta*( cos(ki(xi))+wf*cos(ki(yi)) );
            ek = -2*(h+J*B/2)*(cos(ki(xi)) + cos(ki(yi)))-mu;
            Ek = sqrt(pow(ek,2)+pow(Delta_k,2));
            numspinon += 0.5*(1+ek/Ek)*fF(Ek,T)+0.5*(1-ek/Ek)*(1-fF(Ek,T));
        }
    }
    numspinon *= 2/NN;
    return numspinon;
}
void tJSBMF::step_forward(){
    //calculate the new parameters
    //First determine lambda and h
    if(T == 0){
        //BEC
        lamb = -4*B;
        h = x;
    }
    else{
        double lamb_hi = -4*B;
        double lamb_lo = lamb_hi-10;
        double xn = holon_numb(-4*B);
        if (xn <= x){
            //BEC
            lamb = -4*B;
            h = (x - xn)*NN;
            for(int xi=-N/2; xi<N/2; xi++){
                for(int yi=-N/2; yi<N/2; yi++){
                    if( xi!=0 || yi!=0){
                        h += fB(-2*B*(cos(ki(xi))+cos(ki(yi)))-lamb, T)
                            *(cos(ki(xi))+cos(ki(yi)))/2;
                    }
                }
            }
            h /= NN;
        }
        else{
            while(holon_numb(lamb_lo)>=x){
                lamb_lo -= 50;
            }
            const gsl_root_fsolver_type *sol_type=gsl_root_fsolver_brent;
            gsl_root_fsolver *sol=gsl_root_fsolver_alloc(sol_type);
            gsl_function sol_F;
            //using lambdaFunction to establish F
            auto lambdaFunction = [](double lamb, void *params)->double{
                tJSBMF* ptr_tJ = static_cast<tJSBMF*>(params);
                return ptr_tJ->holon_numb(lamb)-ptr_tJ->x;
            };
            sol_F.function = lambdaFunction;
            sol_F.params = this;
            gsl_root_fsolver_set(sol, &sol_F, lamb_lo, lamb_hi);
            
            //the equation solving process
            int status;
            int iter = 0;
            double lamb_guess;
            do
            {
                iter++;
                status = gsl_root_fsolver_iterate(sol);
                lamb_guess = gsl_root_fsolver_root(sol);
                lamb_lo = gsl_root_fsolver_x_lower(sol);
                lamb_hi = gsl_root_fsolver_x_upper(sol);
                status = gsl_root_test_interval(lamb_lo,lamb_hi, 1e-6, 1e-3);
            } while (status == GSL_CONTINUE && iter < 1000);
            //calculate h
            lamb = lamb_guess;
            h = 0;
            for(int xi=-N/2; xi<N/2; xi++){
                for(int yi=-N/2; yi<N/2; yi++){
                        h += fB(-2*B*(cos(ki(xi))+cos(ki(yi)))-lamb, T)
                            *(cos(ki(xi))+cos(ki(yi)))/2;
                }
            }
            h /= NN;
        }
    }
    //Second determine mu and calculate B&Delta
    
}