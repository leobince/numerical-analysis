#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>
#include <limits> 

const double Pi = acos(-1.0);

class F1 : public Function {
public:
    double operator() (double x) const {
        return 1.0 / x - tan(x);
    }
};

void solve_f1() {
    std::cout << "Solving 1/x - tan(x) on interval [0.1, π/2 - 0.1]" << std::endl;
    Bisection_Method solver_f1(F1(), 0.1, Pi/2 - 0.1);
    double x = solver_f1.solve();
    std::cout << "Root: " << x << std::endl;
}

class F2 : public Function {
public:
    double operator() (double x) const {
        return 1.0 / x - 2 * x;
    }
};

void solve_f2() {
    std::cout << "Solving 1/x - 2x on interval [0.1, 1]" << std::endl;
    Bisection_Method solver_f2(F2(), 0.1, 1);
    double x = solver_f2.solve();
    std::cout << "Root: " << x << std::endl;
}

class F3 : public Function {
public:
    double operator() (double x) const {
        return pow(2, -x) + exp(x) + 2 * cos(x) - 6;
    }
};

void solve_f3() {
    std::cout << "Solving 2^{-x} + e^{x} + 2cos(x) - 6 on interval [1, 3]" << std::endl;
    Bisection_Method solver_f3(F3(), 1, 3);
    double x = solver_f3.solve();
    std::cout << "Root: " << x << std::endl;
}

class F4 : public Function {
public:
    double operator() (double x) const {
        return (pow(x, 3) + 4 * x * x + 3 * x + 5) / (2 * pow(x, 3) - 9 * x * x + 18 * x - 2);
    }
};

void solve_f4() {
    std::cout << "Solving (x^3 + 4x^2 + 3x + 5)/(2x^3 - 9x^2 + 18x - 2) on interval [0, 4]" << std::endl;
    Bisection_Method solver_f4(F4(), 0, 4);
    double x = solver_f4.solve();
    std::cout << "Root: " << x << std::endl;
}

class F5 : public Function {
public:
    double operator() (double x) const {
        return x - tan(x);
    }
};

void solve_f5_newton() {
    std::cout << "Using Newton's method to solve x = tan(x), initial guesses: 4.5 and 7.7" << std::endl;
    Newton_Method solver_f5a(F5(), 4.5);
    double x1 = solver_f5a.solve();
    std::cout << "Root near 4.5: " << x1 << std::endl;

    Newton_Method solver_f5b(F5(), 7.7);
    double x2 = solver_f5b.solve();
    std::cout << "Root near 7.7: " << x2 << std::endl;
}

class F6 : public Function {
public:
    double operator() (double x) const {
        return sin(x / 2) - 1;
    }
};

void solve_f6_secant() {
    std::cout << "Using Secant method to solve sin(x/2) - 1, initial values x0 = 2, x1 = 4" << std::endl;
    Secant_Method solver_f6(F6(), 2, 4);
    double x = solver_f6.solve();
    std::cout << "Root: " << x << std::endl;
    std::cout << "Using Secant method to solve sin(x/2) - 1, initial values x0 = 3, x1 = 5" << std::endl;
    Secant_Method solver_f6_2(F6(), 3, 5);
    double y = solver_f6_2.solve();
    std::cout << "Root: " << y << std::endl;
}

class F7 : public Function {
public:
    double operator() (double x) const {
        return exp(x) - tan(x);
    }
};

void solve_f7_secant() {
    std::cout << "Using Secant method to solve e^{x} - tan(x), initial values x0 = 1, x1 = 1.4" << std::endl;
    Secant_Method solver_f7(F7(), 1, 1.4);
    double x = solver_f7.solve();
    std::cout << "Root: " << x << std::endl;
    std::cout << "Using Secant method to solve e^{x} - tan(x), initial values x0 = 2, x1 = 3" << std::endl;
    Secant_Method solver_f7_2(F7(), 2, 3);
    double y = solver_f7_2.solve();
    std::cout << "Root: " << y << std::endl;
}

class F8 : public Function {
public:
    double operator() (double x) const {
        return x * x * x - 12 * x * x + 3 * x + 1;
    }
};

void solve_f8_secant() {
    std::cout << "Using Secant method to solve x^3 - 12x^2 + 3x + 1, initial values x0 = 0, x1 = -0.5" << std::endl;
    Secant_Method solver_f8(F8(), 0, -0.5);
    double x = solver_f8.solve();
    std::cout << "Root: " << x << std::endl;
    std::cout << "Using Secant method to solve x^3 - 12x^2 + 3x + 1, initial values x0 = 1, x1 = 0.5" << std::endl;
    Secant_Method solver_f8_2(F8(), 1, 0.5);
    double y = solver_f8_2.solve();
    std::cout << "Root: " << y << std::endl;
}

class TroughVolumeFunction : public Function {
private:
    double L;
    double r;
    double V;
public:
    TroughVolumeFunction(double L, double r, double V) : L(L), r(r), V(V) {}

    double operator() (double h) const {
        // Ensure h is within a valid range [0, r]
        if (h < 0 || h > r) {
            return std::numeric_limits<double>::infinity();
        }

        double theta = acos(h / r);
        double area = 0.5 * Pi * r * r - r * r * theta - h * sqrt(r * r - h * h);
        return L * area - V;
    }
};

void solve_trough_volume() {
    std::cout << "Solving for the depth h of water in the trough using three methods" << std::endl;
    double L = 10;
    double r = 1;
    double V = 12.4;
    TroughVolumeFunction F(L, r, V);

    // Bisection method
    Bisection_Method bisection_solver(F, 0.0, r);
    double h_bisect = bisection_solver.solve();
    std::cout << "Bisection method h = " << h_bisect << std::endl;

    // Newton's method
    Newton_Method newton_solver(F, r / 2.0);
    double h_newton = newton_solver.solve();
    std::cout << "Newton's method h = " << h_newton << std::endl;

    // Secant method
    Secant_Method secant_solver(F, 0.1, r);
    double h_secant = secant_solver.solve();
    std::cout << "Secant method h = " << h_secant << std::endl;
}

class NoseInFailureFunction : public Function {
private:
    double A, B, C, E;
public:
    NoseInFailureFunction(double A, double B, double C, double E) : A(A), B(B), C(C), E(E) {}

    double operator() (double alpha_deg) const {
        double alpha = alpha_deg * Pi / 180.0;
        return A * sin(alpha) * cos(alpha) + B * sin(alpha) * sin(alpha) - C * cos(alpha) - E * sin(alpha);
    }
};

void solve_nose_in_failure(double l, double h, double D, double beta1_deg, double initial_guess) {
    double beta1 = beta1_deg * Pi / 180.0;
    double A = l * sin(beta1);
    double B = l * cos(beta1);
    double C = (h + 0.5 * D) * sin(beta1) - 0.5 * D * tan(beta1);
    double E = (h + 0.5 * D) * cos(beta1) - 0.5 * D;

    NoseInFailureFunction F(A, B, C, E);
    Newton_Method newton_solver(F, initial_guess);
    double alpha = newton_solver.solve();
    std::cout <<"When initial_guess = " << initial_guess << "     Alpha (degrees) = " << alpha << std::endl;
}

void solve_nose_in_failure_secant(double l, double h, double D, double beta1_deg, double x0, double x1) {
    double beta1 = beta1_deg * Pi / 180.0;
    double A = l * sin(beta1);
    double B = l * cos(beta1);
    double C = (h + 0.5 * D) * sin(beta1) - 0.5 * D * tan(beta1);
    double E = (h + 0.5 * D) * cos(beta1) - 0.5 * D;

    NoseInFailureFunction F(A, B, C, E);
    Secant_Method secant_solver(F, x0, x1);
    double alpha = secant_solver.solve();
    std::cout <<"When initial_guess = " << x1 << "      Alpha (degrees) = " << alpha << std::endl;
}

int main() {
    solve_f1();
    solve_f2();
    solve_f3();
    solve_f4();
    solve_f5_newton();
    solve_f6_secant();
    solve_f7_secant();
    solve_f8_secant();
    solve_trough_volume();

    // Problem F(a)
    std::cout << "Problem F(a):" << std::endl;
    solve_nose_in_failure(89, 49, 55, 11.5, 33.0);

    // Problem F(b)
    std::cout << "Problem F(b):" << std::endl;
    solve_nose_in_failure(89, 49, 30, 11.5, 33.0);

    // Problem F(c)
    std::cout << "Problem F(c):" << std::endl;
    // Using Secant method with initial guesses far from 33 degrees
    solve_nose_in_failure_secant(89, 49, 30, 11.5, 10.0, 60.0);
    solve_nose_in_failure_secant(89, 49, 30, 11.5, 10.0, 70.0);
    solve_nose_in_failure_secant(89, 49, 30, 11.5, 10.0, 80.0);

    return 0;
}