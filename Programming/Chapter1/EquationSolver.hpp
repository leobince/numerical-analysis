#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER

#include "Function.hpp"
#include <cmath>
#include <stdexcept>

class EquationSolver{
protected:
    const Function & F;
public:
    EquationSolver(const Function& F) : F(F) {}
    virtual double solve() = 0;
};

class Bisection_Method : public EquationSolver {
private:
    double a, b;
    double eps, delta;
    int Maxiter;
public:
    Bisection_Method(const Function &F, double a, double b, 
        double eps = 1e-7, double delta = 1e-6, int Maxiter = 50) :
        EquationSolver(F), a(a), b(b), eps(eps), delta(delta), Maxiter(Maxiter) {}
    
    virtual double solve() {
        double fa = F(a);
        double fb = F(b);
        if (fa * fb > 0) {
            throw std::invalid_argument("Function has the same sign at endpoints. May not have a root.");
        }
        double c, fc;
        for (int iter = 0; iter < Maxiter; ++iter) {
            c = (a + b) / 2;
            fc = F(c);
            if (fabs(fc) < delta || (b - a)/2 < eps) {
                return c;
            }
            if (fa * fc < 0) {
                b = c;
                fb = fc;
            } else {
                a = c;
                fa = fc;
            }
        }
        return c;
    }
};

class Newton_Method : public EquationSolver {
private:
    double x0;
    double eps;
    int Maxiter;
public:
    Newton_Method(const Function &F, double x0, 
        double eps = 1e-7, int Maxiter = 50) :
        EquationSolver(F), x0(x0), eps(eps), Maxiter(Maxiter) {}
    
    virtual double solve() {
        double x = x0;
        for (int iter = 0; iter < Maxiter; ++iter) {
            double fx = F(x);
            double dfx = F.derivative(x);
            if (fabs(dfx) < 1e-12) {
                throw std::runtime_error("Derivative too small, cannot continue iteration.");
            }
            double dx = fx / dfx;
            x -= dx;
            if (fabs(dx) < eps) {
                return x;
            }
        }
        return x;
    }
};

class Secant_Method : public EquationSolver {
private:
    double x0, x1;
    double eps;
    int Maxiter;
public:
    Secant_Method(const Function &F, double x0, double x1,
        double eps = 1e-7, int Maxiter = 50) :
        EquationSolver(F), x0(x0), x1(x1), eps(eps), Maxiter(Maxiter) {}

    virtual double solve() {
        double x_prev = x0;
        double x_curr = x1;
        double f_prev = F(x_prev);
        double f_curr = F(x_curr);
        for (int iter = 0; iter < Maxiter; ++iter) {
            if (fabs(f_curr - f_prev) < 1e-12) {
                // If denominator is too small, adjust f_curr slightly
                f_curr += 1e-12;
            }
            double x_next = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev);
            if (fabs(x_next - x_curr) < eps) {
                return x_next;
            }
            x_prev = x_curr;
            x_curr = x_next;
            f_prev = f_curr;
            f_curr = F(x_curr);
        }
        return x_curr;
    }
};

#endif