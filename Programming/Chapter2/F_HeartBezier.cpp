#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "Polynomial.h"          // Assuming we may need the Polynomial class
#include "InterpCondition.h"     // and other classes if necessary
#include "Interp.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Function to generate points on the heart curve
void generateHeartCurve(std::vector<double>& x_vals, std::vector<double>& y_vals, int m) {
    for (int i = 0; i <= m; ++i) {
        double t = M_PI * i / m;
        double x = 16 * pow(sin(t), 3);
        double y = 13 * cos(t) - 5 * cos(2 * t) - 2 * cos(3 * t) - cos(4 * t);
        x_vals.push_back(x);
        y_vals.push_back(y);
    }
}

// Function to write points to file
void writePointsToFile(const std::vector<double>& x_vals, const std::vector<double>& y_vals, const std::string& filename) {
    std::ofstream file(filename);
    for (size_t i = 0; i < x_vals.size(); ++i) {
        file << x_vals[i] << " " << y_vals[i] << std::endl;
    }
    file.close();
}

int main() {
    std::vector<int> m_values = {10, 40, 160};

    for (int m : m_values) {
        std::vector<double> x_vals;
        std::vector<double> y_vals;

        // Generate heart curve points
        generateHeartCurve(x_vals, y_vals, m);

        // Write to file
        std::string filename = "heart_m" + std::to_string(m) + ".txt";
        writePointsToFile(x_vals, y_vals, filename);
    }

    std::cout << "Heart curve data files generated for m = 10, 40, 160." << std::endl;
    std::cout << "Use a plotting tool to visualize the results." << std::endl;

    return 0;
}