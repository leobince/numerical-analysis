#include <iostream>
#include <vector>
#include <fstream>
#include "Polynomial.h"
#include "InterpCondition.h"
#include "Interp.h"

// Runge function f(x)
double f(double x) {
    return 1.0 / (1.0 + x * x);
}

int main() {
    std::vector<int> n_values = {2, 4, 6, 8};

    // Prepare to write data for plotting
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

        // Generate interpolation nodes
        for (int i = 0; i < num_nodes; ++i) {
            x_nodes[i] = -5.0 + 10.0 * i / n; // x_i = -5 + (10i)/n
            y_nodes[i] = f(x_nodes[i]);
        }

        // Create interpolation condition
        InterpCondition cond(x_nodes, y_nodes);

        // Perform interpolation
        Interp interp(cond);

        // Evaluate the polynomial at plotting points
        std::ofstream poly_file("runge_poly_n" + std::to_string(n) + ".txt");
        for (int i = 0; i < num_points; ++i) {
            double p_val = interp.evaluate(x_plot[i]);
            poly_file << x_plot[i] << " " << p_val << std::endl;
        }
        poly_file.close();
    }

    std::cout << "Data files generated for n = 2, 4, 6, 8." << std::endl;
    std::cout << "Use a plotting tool to visualize the results." << std::endl;

    return 0;
}