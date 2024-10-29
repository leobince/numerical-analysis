#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "Polynomial.h"
#include "InterpCondition.h"
#include "Interp.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Runge function scaled f(x)
double f(double x) {
    return 1.0 / (1.0 + 25 * x * x);
}

int main() {
    std::vector<int> n_values = {5, 10, 15, 20};

    // Prepare to write data for plotting
    std::ofstream exact_file("chebyshev_exact.txt");

    // Generate x values for plotting
    int num_points = 1000;
    std::vector<double> x_plot(num_points);
    double x_start = -1.0;
    double x_end = 1.0;
    double dx = (x_end - x_start) / (num_points - 1);
    for (int i = 0; i < num_points; ++i) {
        x_plot[i] = x_start + i * dx;
        exact_file << x_plot[i] << " " << f(x_plot[i]) << std::endl;
    }
    exact_file.close();

    // Loop over different n values
    for (int n : n_values) {
        std::vector<double> x_nodes(n);
        std::vector<double> y_nodes(n);

        // Generate Chebyshev nodes
        for (int i = 0; i < n; ++i) {
            x_nodes[i] = cos(M_PI * (2.0 * i + 1.0) / (2.0 * n)); // Chebyshev nodes
            y_nodes[i] = f(x_nodes[i]);
        }

        // Create interpolation condition
        InterpCondition cond(x_nodes, y_nodes);

        // Perform interpolation
        Interp interp(cond);

        // Evaluate the polynomial at plotting points
        std::ofstream poly_file("chebyshev_poly_n" + std::to_string(n) + ".txt");
        for (int i = 0; i < num_points; ++i) {
            double p_val = interp.evaluate(x_plot[i]);
            poly_file << x_plot[i] << " " << p_val << std::endl;
        }
        poly_file.close();
    }

    std::cout << "Data files generated for n = 5, 10, 15, 20." << std::endl;
    std::cout << "Use a plotting tool to visualize the results." << std::endl;

    return 0;
}