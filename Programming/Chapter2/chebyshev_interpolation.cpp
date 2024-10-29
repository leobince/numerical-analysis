#include <iostream>
#include <vector>
#include <fstream>

#define _USE_MATH_DEFINES

#include <cmath>

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

// Runge function scaled f(x)
double f(double x) {
    return 1.0 / (1.0 + 25 * x * x);
}

int main() {
    std::vector<int> n_values = {5, 10, 15, 20};

    // Prepare to write data for plotting
    std::ofstream data_file("chebyshev_exact.txt");

    // Generate x values for plotting
    int num_points = 1000;
    std::vector<double> x_plot(num_points);
    double x_start = -1.0;
    double x_end = 1.0;
    double dx = (x_end - x_start) / (num_points - 1);
    for (int i = 0; i < num_points; ++i) {
        x_plot[i] = x_start + i * dx;
        data_file << x_plot[i] << " " << f(x_plot[i]) << std::endl;
    }
    data_file.close();

    // Loop over different n values
    for (int n : n_values) {
        std::vector<double> x_nodes(n);
        std::vector<double> y_nodes(n);
        std::vector<double> coef(n);

        // Generate Chebyshev nodes
        for (int i = 0; i < n; ++i) {
            x_nodes[i] = cos(M_PI * (2 * i + 1) / (2 * n)); // Chebyshev nodes
            y_nodes[i] = f(x_nodes[i]);
        }

        // Compute divided differences coefficients
        computeDividedDifferences(x_nodes, y_nodes, coef);

        // Evaluate the polynomial at plotting points
        std::ofstream poly_file("chebyshev_poly_n" + std::to_string(n) + ".txt");
        for (int i = 0; i < num_points; ++i) {
            double p_val = evaluateNewtonPolynomial(x_nodes, coef, x_plot[i]);
            poly_file << x_plot[i] << " " << p_val << std::endl;
        }
        poly_file.close();
    }

    std::cout << "Data files generated for n = 5, 10, 15, 20." << std::endl;
    std::cout << "Use a plotting tool to visualize the results." << std::endl;

    return 0;
}