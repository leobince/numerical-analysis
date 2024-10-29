#include "Interp.h"
#include <iostream>

// Helper function to compute factorial (used in Hermite interpolation)
unsigned long factorial(unsigned long n) {
    unsigned long result = 1;
    for (unsigned long i = 1; i <= n; ++i) {
        result *= i;
    }
    return result;
}

// Constructor
Interp::Interp(const InterpCondition& con) : condition_(con) {
    computeDividedDifferences();
    interPoly_ = getInterpolationPolynomial();
}

// Compute divided differences
void Interp::computeDividedDifferences() {
    const std::vector<double>& x = condition_.getX();
    const std::vector<double>& y = condition_.getY();
    const std::vector<double>& y_prime = condition_.getYPrime();
    size_t n = x.size();

    // Initialize the table
    tableOfDividedDiffs.resize(n);
    for (size_t i = 0; i < n; ++i) {
        tableOfDividedDiffs[i].resize(n, 0.0);
    }

    // Fill the first column with y values
    for (size_t i = 0; i < n; ++i) {
        tableOfDividedDiffs[i][0] = y[i];
    }

    // If Hermite interpolation is needed (derivatives are provided)
    if (condition_.getOrder() == 1) {
        // Handle the first divided differences
        for (size_t i = 1; i < n; ++i) {
            if (x[i] == x[i - 1]) {
                tableOfDividedDiffs[i][1] = y_prime[i]; // Use derivative when x_i = x_{i-1}
            } else {
                tableOfDividedDiffs[i][1] = (tableOfDividedDiffs[i][0] - tableOfDividedDiffs[i - 1][0]) / (x[i] - x[i - 1]);
            }
        }

        // Compute higher-order divided differences
        for (size_t j = 2; j < n; ++j) {
            for (size_t i = j; i < n; ++i) {
                if (x[i] == x[i - j]) {
                    // Avoid division by zero
                    tableOfDividedDiffs[i][j] = y_prime[i] / factorial(j);
                } else {
                    tableOfDividedDiffs[i][j] = (tableOfDividedDiffs[i][j - 1] - tableOfDividedDiffs[i - 1][j - 1]) / (x[i] - x[i - j]);
                }
            }
        }
    } else {
        // Compute the divided differences for Newton interpolation
        for (size_t j = 1; j < n; ++j) {
            for (size_t i = j; i < n; ++i) {
                tableOfDividedDiffs[i][j] = (tableOfDividedDiffs[i][j - 1] - tableOfDividedDiffs[i - 1][j - 1]) / (x[i] - x[i - j]);
            }
        }
    }

    // Extract the coefficients
    coefficients.resize(n);
    for (size_t i = 0; i < n; ++i) {
        coefficients[i] = tableOfDividedDiffs[i][i];
    }
}

// Get the interpolation polynomial
Polynomial Interp::getInterpolationPolynomial() {
    const std::vector<double>& x = condition_.getX();
    Polynomial poly({coefficients[0]}); // Initialize with the first term

    Polynomial term({1.0});
    for (size_t i = 1; i < coefficients.size(); ++i) {
        term *= Polynomial({-x[i - 1], 1.0}); // (x - x_i)
        poly += term * coefficients[i];
    }

    return poly;
}

// Evaluate the interpolation polynomial at x
double Interp::evaluate(double x_val) {
    double result = 0.0;
    double term = 1.0;
    const std::vector<double>& x = condition_.getX();
    for (size_t i = 0; i < coefficients.size(); ++i) {
        result += coefficients[i] * term;
        term *= (x_val - x[i]);
    }
    return result;
}