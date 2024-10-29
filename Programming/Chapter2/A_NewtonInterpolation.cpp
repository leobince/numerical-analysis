#include <iostream>
#include "Polynomial.h"
#include "InterpCondition.h"
#include "Interp.h"

// Example function f(x)
double f(double x) {
    // Define the function you want to interpolate
    return 1.0 / (1.0 + x * x); // Example: f(x) = 1 / (1 + x^2)
}

int main() {
    int n = 8; // Number of data points
    std::vector<double> x(n);
    std::vector<double> y(n);

    // Sample data points
    for (int i = 0; i <= n; ++i) {
        x[i] = -5.0 + 10.0 * i / n; // Equally spaced points in [-5, 5]
        y[i] = f(x[i]);
    }

    // Create interpolation condition
    InterpCondition cond(x, y);

    // Perform interpolation
    Interp interp(cond);

    // Evaluate the polynomial at a given x value
    double x_val = 2.0;
    double p_val = interp.evaluate(x_val);

    std::cout << "Interpolated value at x = " << x_val << " is " << p_val << std::endl;
    std::cout << "Actual value f(x) = " << f(x_val) << std::endl;

    return 0;
}