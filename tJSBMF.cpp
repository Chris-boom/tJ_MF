#include <iostream>
#include <algorithm>
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
    x0 = 0;
    std::cout<<"Parameters: J="<<J<<" x="<<x<<" T="<<T<<" N="<<N<<" wave="<<wave<<std::endl;
}
double tJSBMF::fF(const double x, const double Temp) const{
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
double tJSBMF::fB(const double x, const double Temp) const{
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
double tJSBMF::holon_numb(double lamb) const{
    //only when T!=0, or lamb must be the lowest energy level.
    double numholon = 0;
    if(T==0){
        std::cout<<"can not be used under nonzero temperature."<<std::endl;
        exit(22333);
    }
    for(int xi=-N/2; xi<N/2; xi++){
        for(int yi=-N/2; yi<N/2; yi++){
            if(xi!=0 || yi!=0){
                numholon += fB(-2*B*(cos(ki(xi))+cos(ki(yi)))-lamb, T);
            }
        }
    }
    if(lamb < -4*B){
        numholon += fB(-4*B-lamb, T);
    }
    numholon /= NN;
    return numholon;
}
double tJSBMF::spinon_numb(double mu) const{
    //calculate the spinon number
    double numspinon = 0;
    double ek=0;
    double Ek=0;
    double Delta_k=0;
    for(int xi=-N/2; xi<N/2; xi++){
        for(int yi=-N/2; yi<N/2; yi++){
            Delta_k = Delta*( cos(ki(xi))+wf*cos(ki(yi)) );
            ek = -2*(h+J*B)*(cos(ki(xi)) + cos(ki(yi)))-mu;
            Ek = sqrt(pow(ek,2)+pow(Delta_k,2));
            if(Ek != 0){
                numspinon += 0.5*(1+ek/Ek)*fF(Ek,T)+0.5*(1-ek/Ek)*(1-fF(Ek,T));
            }
            else{
                numspinon += fF(ek,T);
            }
        }
    }
    numspinon = numspinon/NN*2;
    return numspinon;
}
void tJSBMF::step_forward(){
    //calculate the new parameters
    //First determine lambda and h
    const gsl_root_fsolver_type *sol_type=gsl_root_fsolver_brent;
    gsl_root_fsolver *sol=gsl_root_fsolver_alloc(sol_type);
    gsl_function sol_F;
    int status;
    int iter;
    if(T == 0){
        //BEC
        lamb = -4*B;
        h_new = x;
        x0 = x;
    }
    else{
        double lamb_hi = -4*B;
        double lamb_lo = lamb_hi-10;
        double xn = holon_numb(-4*B);
        if (xn <= x){
            //BEC
            lamb = -4*B;
            h_new = (x - xn)*NN;
            for(int xi=-N/2; xi<N/2; xi++){
                for(int yi=-N/2; yi<N/2; yi++){
                    if( xi!=0 || yi!=0){
                        h_new += fB(-2*B*(cos(ki(xi))+cos(ki(yi)))-lamb, T)
                            *(cos(ki(xi))+cos(ki(yi)))/2;
                    }
                }
            }
            h_new /= NN;
            x0 = x - xn;
        }
        else{
            while(holon_numb(lamb_lo)>=x){
                lamb_lo -= 50;
            }
            //using lambdaFunction to establish F
            auto lambdaFunction = [](double lamb, void *params)->double{
                tJSBMF* ptr_tJ = static_cast<tJSBMF*>(params);
                return ptr_tJ->holon_numb(lamb)-ptr_tJ->x;
            };
            sol_F.function = lambdaFunction;
            sol_F.params = this;
            gsl_root_fsolver_set(sol, &sol_F, lamb_lo, lamb_hi);
            
            //the equation solving process
            double lamb_guess;
            iter = 0;
            do
            {
                iter++;
                status = gsl_root_fsolver_iterate(sol);
                lamb_guess = gsl_root_fsolver_root(sol);
                lamb_lo = gsl_root_fsolver_x_lower(sol);
                lamb_hi = gsl_root_fsolver_x_upper(sol);
                status = gsl_root_test_interval(lamb_lo,lamb_hi, atol, rtol);
            } while (status == GSL_CONTINUE && iter < 1000);
            //calculate h
            lamb = lamb_guess;
            h_new = 0;
            for(int xi=-N/2; xi<N/2; xi++){
                for(int yi=-N/2; yi<N/2; yi++){
                        h_new += fB(-2*B*(cos(ki(xi))+cos(ki(yi)))-lamb, T)
                            *(cos(ki(xi))+cos(ki(yi)))/2;
                }
            }
            h_new /= NN;
            x0 = 0;
        }
    }
    //Second determine mu and calculate B&Delta
    auto muFunction = [](double mu, void *params)->double{
        tJSBMF *ptr_tJ = static_cast<tJSBMF*>(params);
        return ptr_tJ->spinon_numb(mu)-(1-ptr_tJ->x);
    };
    sol_F.function = muFunction;
    sol_F.params = this;
    double mu_lo = -10;
    double mu_hi = 10;
    double mu_guess = 0;
    while( (spinon_numb(mu_lo)-1+x)*(spinon_numb(mu_hi)-1+x)>=0){
        mu_lo -= 100;
        mu_hi += 100;
    }
    gsl_root_fsolver_set(sol, &sol_F, mu_lo, mu_hi);
    iter = 0;
    do{
        iter++;
        status = gsl_root_fsolver_iterate(sol);
        mu_guess = gsl_root_fsolver_root(sol);
        mu_lo = gsl_root_fsolver_x_lower(sol);
        mu_hi = gsl_root_fsolver_x_upper(sol);
        status = gsl_root_test_interval(mu_lo, mu_hi, atol, rtol);
    }while(status == GSL_CONTINUE && iter < 1000);
    //calculate B&Delta
    mu = mu_guess;
    B_new = 0;
    Delta_new = 0;
    // double Delta_y = 0;
    // double Delta_x = 0;
    double Delta_k, ek, Ek;
    for(int xi=-N/2; xi<N/2; xi++){
        for(int yi=-N/2; yi<N/2; yi++){
            Delta_k = Delta*( cos(ki(xi))+wf*cos(ki(yi)) );
            ek = -2*(h+J*B)*(cos(ki(xi)) + cos(ki(yi)))-mu;
            Ek = sqrt(pow(ek,2)+pow(Delta_k,2));
            B_new += 0.5*(0.5*(1+ek/Ek)*fF(Ek,T)+0.5*(1-ek/Ek)*(1-fF(Ek,T)))*
                (cos(ki(xi))+cos(ki(yi)));
            Delta_new += 0.5*Delta_k*( 1 - 2*fF(Ek,T) )/Ek*
                (cos(ki(xi)) + wf*cos(ki(yi)));
            // Delta_x += 0.5*Delta_k*( 1 - 2*fF(Ek,T) )/Ek*
            //     2*cos(ki(xi));
            // Delta_y += 0.5*Delta_k*( 1 - 2*fF(Ek,T) )/Ek*
            //     2*cos(ki(yi));
        }
    }
    B_new /= NN;
    Delta_new = Delta_new/NN*J;
    // Delta_y = Delta_y/NN*J;
    // Delta_x = Delta_x/NN*J;
    // std::cout<<"Delta_y="<<Delta_y<<" Delta_x="<<Delta_x<<std::endl;
}
void tJSBMF::self_consistent(int report_freq){
    //the err of every param must be less than its etol.
    //Remember to set the start state!
    int iter = 0;
    int status_MF = 1;
    do
    {
        iter++;
        step_forward();
        if(abs(h_new-h)<abs(h)*rtol+atol){
            if(abs(B_new-B)<abs(B)*rtol+atol){
                if(abs(Delta_new-Delta)<abs(Delta)*rtol+atol){
                    status_MF = 0;
                }
                else status_MF = -1;
            }
            else status_MF = -2;
        }
        else status_MF = -3;
        if(report_freq!=0 && (iter%report_freq == 0 || iter == 1)){
            std::cout<<"step"<<iter<<"\th="<<h<<"\tB="<<B<<"\tDelta="<<Delta<<std::endl;
        }
        h = h_new*sfac + h*(1-sfac);
        B = B_new*sfac + B*(1-sfac);
        Delta = Delta_new*sfac + Delta*(1-sfac);
    } while (status_MF != 0 && iter < maxstep);
}
int tJSBMF::reset(){
    //initialize the parameters
    h = 0.1;
    B = 1;
    Delta = 1;
    return 0;
}
double tJSBMF::DeltaSC(){
    double DeltaSC = Delta*x0;
    return DeltaSC;
}