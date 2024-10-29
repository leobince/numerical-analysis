#include <iostream>
#include <vector>
#include <fstream>
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

// Runge function f(x)
double f(double x) {
    return 1.0 / (1.0 + x * x);
}

int main() {
    std::vector<int> n_values = {2, 4, 6, 8};

    // Prepare to write data for plotting
    std::ofstream data_file("runge_data.txt");
    std::ofstream exact_file("runge_exact.txt");

    // Generate x values for plotting
    int num_points = 1000;
    std::vector<double> x_plot(num_points);
    double x_start = -5.0;
    double x_end = 5.0;
    double dx = (x_end - x_start) / (num_points - 1);
    for (int i = 0; i < num_points; ++i) {
        x_plot[i] = x_start + i * dx;
        exact_file << x_plot[i] << " " << f(x_plot[i]) << std::endl;
    }
    exact_file.close();

    // Loop over different n values
    for (int n : n_values) {
        int num_nodes = n + 1;
        std::vector<double> x_nodes(num_nodes);
        std::vector<double> y_nodes(num_nodes);
        std::vector<double> coef(num_nodes);

        // Generate interpolation nodes
        for (int i = 0; i < num_nodes; ++i) {
            x_nodes[i] = -5.0 + 10.0 * i / n; // x_i = -5 + (10i)/n
            y_nodes[i] = f(x_nodes[i]);
        }

        // Compute divided differences coefficients
        computeDividedDifferences(x_nodes, y_nodes, coef);

        // Evaluate the polynomial at plotting points
        std::ofstream poly_file("runge_poly_n" + std::to_string(n) + ".txt");
        for (int i = 0; i < num_points; ++i) {
            double p_val = evaluateNewtonPolynomial(x_nodes, coef, x_plot[i]);
            poly_file << x_plot[i] << " " << p_val << std::endl;
        }
        poly_file.close();
    }

    std::cout << "Data files generated for n = 2, 4, 6, 8." << std::endl;
    std::cout << "Use a plotting tool to visualize the results." << std::endl;

    return 0;
}