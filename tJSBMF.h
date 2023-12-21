#include <vector>
#include <cmath>

class tJSBMF
{
public:
    tJSBMF(double J, double x, double T, int N, char wave='d',
        int maxstep=1000, double atol=1e-6, double rtol=1e-3);
    void step_forward();
    void self_consistent();
    double holon_numb(double lamb) const;
    double spinon_numb(double mu) const;
    inline double ki(int i) const{
        return 2*M_PI*i/N;
    }
    double J, x, T;
    int N;
    char wave;
    double wf;
    double Delta, h, B, lamb, mu;
private:
    double fF(const double x, const double Temp) const;
    double fB(const double x, const double Temp) const;
    double Delta_new, h_new, B_new;
    int NN, maxstep;
    double atol, rtol;
    double sfac = 0.7;
};