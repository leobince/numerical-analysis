#include <iostream>
#include <vector>

// Function to compute divided differences coefficients
void computeDividedDifferences(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& coef) {
    int n = x.size();
    coef = y; // Initialize coefficients with y values
    for (int j = 1; j < n; ++j) {
        for (int i = n - 1; i >= j; --i) {
            coef[i] = (coef[i] - coef[i - 1]) / (x[i] - x[i - j]);
        }
    }
}

// Function to evaluate Newton's interpolating polynomial at x_val
double evaluateNewtonPolynomial(const std::vector<double>& x, const std::vector<double>& coef, double x_val) {
    int n = x.size();
    double result = coef[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        result = result * (x_val - x[i]) + coef[i];
    }
    return result;
}

// Example function f(x)
double f(double x) {
    // Define the function you want to interpolate
    return 1.0 / (1.0 + x * x); // Example: f(x) = 1 / (1 + x^2)
}

int main() {
    // Example usage
    int n = 5; // Number of data points
    std::vector<double> x(n);
    std::vector<double> y(n);
    std::vector<double> coef(n);

    // Sample data points
    for (int i = 0; i < n; ++i) {
        x[i] = -5.0 + 10.0 * i / (n - 1); // Equally spaced points in [-5, 5]
        y[i] = f(x[i]);
    }

    // Compute divided differences coefficients
    computeDividedDifferences(x, y, coef);

    // Evaluate the polynomial at a given x value
    double x_val = 2.0;
    double p_val = evaluateNewtonPolynomial(x, coef, x_val);

    std::cout << "Interpolated value at x = " << x_val << " is " << p_val << std::endl;
    std::cout << "Actual value f(x) = " << f(x_val) << std::endl;

    return 0;
}